"""
Extracts templates from reaction datasets (USPTO and Pistachio) using the LocalTemplate method
from https://github.com/kaist-amsg/LocalRetro. Their scripts are recreated in /LocalTemplate for clarity
and reproducibility.
"""
import pandas as pd
import argparse
from LocalTemplate.extract_from_train_data import build_template_extractor, extract_templates


def extract_localtemplates(dataset_name, args):
    extractor = build_template_extractor(args)
    
    df = pd.read_csv(f"./data/raw/2_{dataset_name}.csv")
    rxns = df['rxn_map'].tolist()
    
    template_infos, template_labels = extract_templates(rxns, args, extractor)
    
    df['template'] = template_labels
    print(f"Saved LocalTemplates to ./data/raw/2_{dataset_name}.csv")
    
    template_infos.to_csv(f"./data/raw/{dataset_name}_template_info.csv")
    print(f"Saved template infos to ./data/raw/2_{dataset_name}_template_info.csv")

def main(args):
    print("Extracting LocalTemplates for USPTO dataset...")
    extract_localtemplates("uspto", args)
    
    if not args["uspto_only"]:
        print("Extracting LocalTemplates for Pistachio dataset...")
        extract_localtemplates("pistachio", args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--uspto-only", action="store_true", help="Only run USPTO code.")
    parser.add_argument("--retro", default=True, help='Retrosynthesis or forward synthesis (True for retrosynthesis)')
    parser.add_argument('-v', '--verbose', default=False,  help='Verbose during template extraction')
    parser.add_argument('-stereo', '--use-stereo', default=True,  help='Use stereo info in template extraction')
    parser.add_argument('-min', '--min-template-n', type=int, default=1,  help='Minimum of template frequency')
    args = vars(parser.parse_args())
    main(args)
    print("done!")