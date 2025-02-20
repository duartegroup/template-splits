"""
Cleans and filters atom-mapped reaction strings for the initialised USPTO and 
Pistachio datasets.
"""

import os
import requests
import tarfile
import pandas as pd
import numpy as np
import re
import tarfile
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts, ReactionToSmiles
from rdkit.Chem import Draw
from rxn.chemutils.reaction_equation import (
    canonicalize_compounds,
    merge_reactants_and_agents,
    sort_compounds,
    cleanup_compounds
)
from rxn.chemutils.reaction_smiles import (
    ReactionFormat,
    parse_extended_reaction_smiles,
    to_reaction_smiles,
)
from rxn.chemutils.utils import remove_atom_mapping
from pandarallel import pandarallel
from rdkit.Chem.rdMolDescriptors import CalcNumHeavyAtoms
import argparse

pandarallel.initialize(progress_bar=True, nb_workers=4)


def untangle_tildes(rxn_map):
    """Untangles molecules with coordinate bonds ~ in reaction maps"""
    rxn = rxn_map.split(">")
    new_rxn = []
    for side in rxn:
        new_side = []
        side = side.split(".")
        for mol in side:
            if "~" in mol:
                try:
                    mols = Chem.MolFromSmiles(mol)
                    tilde_list = []
                    for b in mols.GetBonds():
                        if str(b.GetBondType()) == "UNSPECIFIED":
                            tilde_list.append(b.GetIdx())
                    mols = Chem.FragmentOnBonds(mols, tilde_list)
                    smiles = Chem.MolToSmiles(mols)
                    smiles = re.sub(r'\[\d+\*\]', '', smiles).replace("~","").replace("*","").replace("()","")
                    #smiles = smiles.replace("~",".").split(".")
                    #smiles = ".".join([m for m in smiles if "*" not in m])
                except:
                    smiles = mol.replace("~",".")
            else:
                smiles = mol
            new_side.append(smiles)
        new_rxn.append(".".join(new_side))
    
    full_rxn = new_rxn[0]+">"+new_rxn[1]+">"+new_rxn[2]
    return full_rxn


def join_reactants_reagents(rxn_map):
    """Joins reactants and reagents in a reaction map."""
    if ">>" in rxn_map:
        return rxn_map
    else:
        rxn = rxn_map.split(">")
        full_rxn = rxn[0]+"."+rxn[1]+">>"+rxn[2]
        return full_rxn


def remove_fragment_info(rxn_map):
    """Removes additional info at the end of a reaction string."""
    return rxn_map.split(" ")[0]


def remove_reagents(rxn_map):
    reactants = rxn_map.split(">>")[0].split(".")
    product = rxn_map.split(">>")[1]
    product_maps = set(re.findall('\:([0-9]+)\]', product))
    
    new_reactants = []
    for smiles in reactants:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        used = False
        for a in mol.GetAtoms():
            if a.HasProp('molAtomMapNumber'):
                if a.GetProp('molAtomMapNumber') in product_maps:
                    used = True
                else:
                    a.ClearProp('molAtomMapNumber')
        if used:
            new_reactants.append(Chem.MolToSmiles(mol, True))
    new_rxn_map = ".".join(new_reactants) + ">>" + product
    return new_rxn_map


def reactant_count_filter(rxn_map):
    reactants = rxn_map.split(">>")[0].split(".")
    count = len(reactants)
    if count < 5 and count >= 1:
        return rxn_map
    else:
        return "Error: too many reactants"


def multiproduct_fixer(rxn_map):
    products = rxn_map.split(">>")[1].split(".")
    if len(products) > 1:
        prod_mols = [Chem.MolFromSmiles(prod) for prod in products]
        try:
            prod_atoms = [CalcNumHeavyAtoms(mol) for mol in prod_mols]
        except:
            return "Error: invalid product"
        prod_bools = [atoms > 5 for atoms in prod_atoms]
        new_prods = [products[i] for i, b in enumerate(prod_bools) if b]
        if len(new_prods) == 1:
            new_rxn_map = rxn_map.split(">>")[0] + ">>" + new_prods[0]
            return new_rxn_map
        else:
            return "Error: too many products"
    else:
        if Chem.MolFromSmiles(products[0]):
            return rxn_map
        else:
            return "Error: product smiles error"


def no_carbon(rxn_map):
    """Tags inorganic products."""
    prod = rxn_map.split(">>")[1]
    mol = Chem.MolFromSmiles(prod)
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[#6]")):
        return rxn_map
    else:
        return "Error: no organic product"
    
    
def remove_stereoalchemy(rxn):
    """Removes product stereochemistry if no stereochemistry
    is present in the reactants."""
    r = rxn.split(">>")[0]
    p = rxn.split(">>")[1]
    
    if "@" in p:
        if "@" not in r:
            # Remove chiral centers from product
            p_mol = Chem.MolFromSmiles(p)
            Chem.RemoveStereochemistry(p_mol)
            new_p = Chem.MolToSmiles(p_mol)
            return rxn.replace(p, new_p)
    return rxn


def product_in_reactants(rxn_map):
    product = rxn_map.split(">>")[1]
    reactants = rxn_map.split(">>")[0].split(".")
    if product in reactants:
        return "Error: product in reactants"
    else:
        return rxn_map
    

def canonicalise(rxn_map):
    """Canonicalise reaction CXSMILES. Returns Invalid SMILES if
    SMILES does not contain 2-10 reactants and 1 product. Mixes
    together reactants and reagents.

    Args:
        rxn (str): Reaction CXSMILES

    Returns:
        str: Canonicalised reaction SMILES
    """
    try:
        # Show fragment interactions with ~
        rxn_type = ReactionFormat.STANDARD_WITH_TILDE
        # Parse CXSMILES
        reacteq = parse_extended_reaction_smiles(rxn_map, remove_atom_maps=False)
        # Canonicalise and sort compounds
        std_rxn = sort_compounds(canonicalize_compounds(reacteq))
        # Convert to SMILES
        new_rxn_map = to_reaction_smiles(std_rxn, rxn_type)
        # Check validity with RDKit
        try:
            _ = ReactionFromSmarts(new_rxn_map)
            return new_rxn_map
        except:
            return "Error: smarts"
    except:
        return "Error: during canonicalisation"


def remove_mapping(rxn_map):
    try:
        # Show fragment interactions with ~
        rxn_type = ReactionFormat.STANDARD_WITH_TILDE
        # Parse CXSMILES
        reacteq = parse_extended_reaction_smiles(rxn_map, remove_atom_maps=True)
        # Canonicalise and sort compounds
        reacteq = cleanup_compounds(reacteq)
        std_rxn = sort_compounds(canonicalize_compounds(reacteq))
        # Convert to SMILES
        rxn = to_reaction_smiles(std_rxn, rxn_type)
        # Check validity with RDKit
        try:
            _ = ReactionFromSmarts(rxn)
            return rxn
        except:
            return "Error: smarts"
    except:
        return "Error: mapping"


def size_filter(rxn):
    """Tags reactions longer than 512 tokens."""
    return (len(rxn)>512)
    
    
def clean_dataset(dataset_name):
    df = pd.read_csv(f"data/raw/1_{dataset_name}.csv", index_col="dataset_id")
    df = df.head(1000)
    print("Full dataset size: ", len(df))

    df = df[df["rxn_map"].str.count(">")==2]
    df = df.drop_duplicates("rxn_map")
    print("Rows after removing duplicates: ", len(df))

    print("Untangling and removing tildes...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(untangle_tildes)
    df = df.drop_duplicates("rxn_map")
    print("Rows after untangling: ", len(df))
    
    print("Joining reactants and reagents...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(join_reactants_reagents)
    df = df.drop_duplicates("rxn_map")
    print("Rows after joining reactants and reagents: ", len(df))

    print("Remove extra info...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(remove_fragment_info)
    df = df.drop_duplicates("rxn_map")
    print("Rows after removing fragment info: ", len(df))

    print("Removing reagents")
    df["rxn_map"] = df["rxn_map"].parallel_apply(remove_reagents)
    df = df.drop_duplicates("rxn_map")
    print("Rows after removing reagents: ", len(df))

    print("Removing reactions with >4 reactants")
    df["rxn_map"] = df["rxn_map"].parallel_apply(reactant_count_filter)
    df = df[~df["rxn_map"].str.contains("Error")]
    print("Rows after removing reactions with >4 reactants: ", len(df))
    
    print("Removing multi-product reactions and side products...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(multiproduct_fixer)
    df = df[~df["rxn_map"].str.contains("Error")]
    print("Rows after removing multi-product reactions: ", len(df))

    print("Removing inorganic reactions...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(no_carbon)
    df = df[~df["rxn_map"].str.contains("Error")]
    df = df.drop_duplicates("rxn_map")
    print("Rows after removing inorganic reactions: ", len(df))

    print("Removing stereoalchemy...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(remove_stereoalchemy)
    df = df.drop_duplicates("rxn_map")
    print("Rows after removing stereoalchemy: ", len(df))

    print("Removing reactions with product listed as a reactant...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(product_in_reactants)
    df = df[~df["rxn_map"].str.contains("Error")]
    df = df.drop_duplicates("rxn_map")
    print("Rows after removing product in reactant reactions: ", len(df))

    print("Canonicalising...")
    df["rxn_map"] = df["rxn_map"].parallel_apply(canonicalise)
    df = df[~df["rxn_map"].str.contains("Error")]
    print("Rows after canonicalising: ", len(df))

    print("Removing atom mapping...")
    df["canonic_rxn"] = df["rxn_map"].parallel_apply(remove_mapping)
    df = df[~df["canonic_rxn"].str.contains("Error")]
    print("Rows after removing atom mapping: ", len(df))
    
    print("Removing extra long reactions...")
    df = df[~df["canonic_rxn"].parallel_apply(size_filter)]
    print("Rows after removing extra long reactions: ", len(df))
    
    print("Removing duplicates...")
    df = df.drop_duplicates("canonic_rxn")
    print("Rows after removing duplicates: ", len(df))

    df.to_csv(f"./data/raw/2_{dataset_name}.csv")
    print(f"Saved to ./data/raw/2_{dataset_name}.csv")
    

def main(args):
    print("Loading USPTO...")
    clean_dataset("uspto")
    
    if not args["uspto_only"]:
        print("Loading Pistachio...")
        clean_dataset("pistachio")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--uspto_only", action="store_true", help="Only run USPTO code.")
    args = vars(parser.parse_args())
    main(args)
    print("done!")