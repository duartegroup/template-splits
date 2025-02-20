"""
Initialises USPTO and Pistachio datasets for further preprocessing. The Pistachio dataset
is filtered to remove US patents and remapped using RXNMapper.

Note that the Pistachio dataset must be manually imported into data/raw as pistachio.smi
before running this script.
"""
import os
import requests
import tarfile
import pandas as pd
import tarfile
from rxn.chemutils.utils import remove_atom_mapping

from rxnmapper import RXNMapper
import pandas as pd
from rxn.chemutils.utils import remove_atom_mapping
from tqdm import tqdm
import argparse
from transformers.utils import logging

logging.set_verbosity(40)

def download_uspto(data_dir="data/raw/"):
    """Download USPTO data (after cleaning and atom mapping)."""

    base_url = "https://drive.switch.ch/index.php/s/t9uvLpScEu1esIy/download?path=%2F&files={}"

    fname = "uspto.tar.gz"
    url = base_url.format(fname)

    target_path = data_dir + fname

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    if not os.path.isfile(target_path):
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(target_path, "wb") as f:
                f.write(response.raw.read())

    with tarfile.open(target_path) as f:
        f.extractall(data_dir)


def reformat_uspto():
    """Reformat USPTO dataset from .tsv to .csv."""
    df = pd.read_csv("data/raw/uspto_public.tsv", sep="\t").rename(columns={"rxnmapper_aam":"rxn_map"})
    df.index.name = "dataset_id"
    df.drop(columns=["yield"]).to_csv("data/raw/1_uspto.csv")


def load_pistachio():
    """Load Pistachio dataset from .smi file."""
    filename = "./data/raw/pistachio.smi"
    if not os.path.exists(filename):
        raise FileNotFoundError("No file named './data/raw/pistachio.smi' - must be manually imported!")
    pistachio_full = pd.read_csv(filename, sep="\t", header=None).drop(columns=[2], axis=1)
    pistachio_full = pistachio_full.rename(columns={0:"rxn_map", 1:"patent", 3:"rxn_class", 4:"rxn_name"})
    pistachio_full.index.name = "dataset_id"
    return pistachio_full


def remove_us_patents(df):
    """Filter out all US-based patents"""
    df["region"] = df["patent"].apply(lambda x: x[:2])
    df_no_us = df.loc[df["region"]!="US"]
    df_no_us = df_no_us.drop(columns=["region"])
    return df_no_us


def size_filter(rxn):
    """Tags reactions longer than 512 tokens."""
    rxn = remove_atom_mapping(rxn)
    return (len(rxn)>512)


def remap(rxns):
    """Remaps a list of reactions using RXNMapper."""
    def batch(iterable, n=1):
        l = len(iterable)
        for ndx in range(0, l, n):
            yield iterable[ndx:min(ndx + n, l)]

    rxn_mapper = RXNMapper()
    results = []
    batch_size = 32
    number_batches = (len(rxns) - 1) // batch_size + 1
    for chunk in tqdm(batch(rxns, batch_size), total=number_batches):
        try:
            result = rxn_mapper.get_attention_guided_atom_maps(chunk)
            result = [x["mapped_rxn"] for x in result]
            results += result
        except:
            # Batched RXNMapper function cannot handle errors 
            for rxn in chunk:
                try:
                    result = rxn_mapper.get_attention_guided_atom_maps([rxn])
                except:
                    result = [{"mapped_rxn": ">>"}]
                results.append(result[0]["mapped_rxn"])
    return results
                

def main(args):
    # Initialising USPTO
    download_uspto()
    reformat_uspto()
    print("Initalised USPTO dataset at data/raw/1_uspto.csv")
    
    if not args["uspto_only"]:
        # Initialising Pistachio
        pistachio_full = load_pistachio()
        pistachio_no_us = remove_us_patents(pistachio_full)
        # Filter and drop duplicates
        pistachio_filtered = pistachio_no_us[~pistachio_no_us["rxn_map"].apply(size_filter)]
        pistachio_filtered = pistachio_filtered.drop_duplicates("rxn_map")
        # Remap with RXNMapper
        pistachio_rxns = list(pistachio_filtered["rxn_map"])
        pistachio_filtered["rxn_map"] = remap(pistachio_rxns)
        
        pistachio_filtered.to_csv("./data/raw/1_pistachio.csv")
        print("Initalised Pistachio dataset at data/raw/1_pistachio.csv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--uspto_only", action="store_true", help="Only run USPTO code.")
    args = vars(parser.parse_args())
    main(args)
    print("done!")