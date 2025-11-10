"""
Splits USPTO by template via three strategies:
    Full split: random split into 90:5:5 train:val:test sets, with equivalent products
        grouped together in the same set.
    Narrow template split: generates smaller and less diverse train, val, and test sets
        which are x% of the full split, sampled by template.
    Broad template split: generates smaller train sets by reducing the number of reactions
        but maintaining diversity by keeping at least 1 reaction per template.
        
Splits Pistachio by template:
    ID test: contains 10k reactions with templates also present in USPTO.
    OOD test: contains 10k reactions with templates not present in USPTO.
"""
import pandas as pd
import os
import argparse

SET_NAMES = ["train","val","test"]

def full_dataset_product_split(data, frac):
    """
    Split reactions dataframe in train, val, and test sets keeping reactions
    with equal products in the same set.
    Adapted from https://github.com/schwallergroup/choriso

    Args:
        data: pd.Df, dataset containing the reactions
        frac: float, fraction of original data for val and test sets

    Out:
        [train, val, test]: list, pd.Dfs corresponding to train, val, and test sets
    """

    # Number of reactions in the test set
    tot = round(len(data) * frac)

    # Create groups and shuffle
    groups = data.groupby("products").size().sample(frac=1.0, random_state=42)

    # Create lists of counts and values
    counts = groups.values
    products = groups.index.values

    counter = 0
    test_products = []
    val_products = []

    # Append products in list until having a number of reactions > tot
    test_on = True
    for i, times in enumerate(counts):
        if test_on:
            test_products.append(products[i])
            counter += times
            if counter > tot:
                test_on = False
                counter = 0
        else:
            val_products.append(products[i])
            counter += times
            if counter > tot:
                break

    # Create boolean array to select test reactions
    test_mask = data["products"].isin(test_products)
    val_mask = data["products"].isin(val_products)

    test = data[test_mask]
    val = data[val_mask]
    train = data[~(test_mask+val_mask)]

    return [train, val, test]


def make_dir(data_dir):
    """Make directory if not already present."""
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)


def save_sets(sets, split_name, frac=None):
    """
    Save train, val, and test sets to the allocated directory,
    given by the split_name and frac.
    """
    data_dir = f"data/processed/uspto_retro/splits/{split_name}"
    make_dir(data_dir)
    if frac:
        data_dir = f"data/processed/uspto_retro/splits/{split_name}/{frac}"
        make_dir(data_dir)
        print(f"{split_name} {frac}% split")
    else:
        print(f"{split_name} split")

    for set, set_name in zip(sets, SET_NAMES):
        print(f"   {set_name}: {str(len(set))}")
        set.drop(columns=["products"]).to_csv(f"{data_dir}/{set_name}.csv")


def full_split(data):
    """
    Full dataset split, keeping all reactions with equal products in
    the same split. Splits into 90% train, 5% val, 5% test.
    """
    data["products"] = data["canonic_rxn"].apply(lambda x: x.split(">>")[1])
    sets = full_dataset_product_split(data, 0.05)
    save_sets(sets, "full")
    return sets

        
def broad_template_split(sets):
    """
    Broad dataset split, keeps the same val and test set as the full
    split but reduces the volume of the train set while keeping at least
    1 reaction from each template.
    """
    fracs = [90, 50, 25, 10]
    
    # Retain 1 example per template
    must_train = sets[0].groupby("template").sample(n=1, random_state=42)
    
    # Create splits sequentially to ensure same sampling
    for i in range(3):
        sets[0] = sets[0].groupby("template").sample(frac=fracs[i+1]/fracs[i], random_state=42)
        sets[0] = pd.concat([sets[0], must_train]).drop_duplicates()
        save_sets(sets, "broad", str(fracs[i+1]))


def narrow_template_split(sets):
    """
    Narrow dataset split, reduces the number of templates sampled in the
    train, val, and test splits according to the assigned fraction.
    """
    fracs = [10, 25, 50]    
    groups = sets[0].groupby("template").size().sample(frac=1.0, random_state=40)
    counts = groups.values
    templates = groups.index.values
    for frac in fracs:
        # Number of reactions desired in the train set
        tot = round(len(sets[0]) * frac/90)

        counter = 0
        selected_templates = []
        # Append templates in list until having a number of reactions > tot
        for i, times in enumerate(counts):
            selected_templates.append(templates[i])
            counter += times
            if counter > tot:
                break
        
        # Select reactions in all sets according to templates
        new_sets = [x.loc[x["template"].isin(selected_templates)] for x in sets]
        save_sets(new_sets, "narrow", frac)
        

def template_test(df, df_base):
    """
    Splits external dataset (df) by ID and OOD templates with regards to the original
    dataset (df_base), and samples 10k reactions for each test set.
    """
    df_ood = df[~df["template"].isin(df_base["template"].drop_duplicates())]
    df_id = df[df["template"].isin(df_base["template"].drop_duplicates())]
    
    for split, split_name in zip([df_ood, df_id], ["OOD", "ID"]):
        print("Pistachio", split_name, "reactions:", len(split))
        split_test = split.sample(10000, random_state=42)
        make_dir(f"data/processed/pistachio_retro/pistachio_{split_name}")
        split_test.to_csv(f"data/processed/pistachio_retro/pistachio_{split_name}/test.csv")
    

def main(args):
    # Load USPTO
    uspto = pd.read_csv("data/processed/uspto_retro/uspto_retro.csv", index_col="dataset_id")
    make_dir("data/processed/uspto_retro/splits")
    # USPTO splits
    uspto_sets = full_split(uspto)
    narrow_template_split(uspto_sets)
    broad_template_split(uspto_sets)

    if not args["uspto_only"]:
        # Load Pistachio
        pistachio = pd.read_csv("data/processed/pistachio_retro/pistachio_retro.csv", index_col="dataset_id")
        # Pistachio splits
        template_test(pistachio, uspto)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--uspto_only", action="store_true", help="Only run USPTO code.")
    args = vars(parser.parse_args())
    main(args)
    print("done!")