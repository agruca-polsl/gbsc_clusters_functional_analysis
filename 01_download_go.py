"""
GBSC Clusters GO Ontology Functional Analysis Pipeline
=======================================================

Major refactoring of original analysis scripts by Joanna Ziemska-Legiecka (2025).
GO download logic, clusters GO enrichment algorithms and s-measure caluclations preserved with fixes.

Author: Aleksandra Gruca (2026)
Original: Joanna Ziemska-Legiecka (2025)
"""

GO_ANNOTATIONS_FILE = "go_annotations.json"
GO_NAMES_FILE ="go_names.csv"
GO_MAX_PATH_FILE ="go_max_path.csv"


import logging
import os
import sys
import time
import typing
import json
from pathlib import Path
from optparse import OptionParser

import requests


def fill_names(go_ids, save_file="/tmp/tmp_go.csv"):
    for e, go_id in enumerate(list(go_ids)):
        URL = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query={go_id}"
        print(URL, e, len(go_ids))
        res = requests.get(URL)
        for result in res.json()["results"]:
            if result["id"] == go_id:
                go_name = result['name']
                aspect = result["aspect"]
                print(go_id, (go_name, aspect))
                with open(save_file, "a") as f:
                    f.write(f"{go_id}\t{go_name}\t{aspect}\n")
                break


def get_GO(
        protein_list: iter,
        exclude: list,
        #save_go_file: str,
        aspect: str,
        #lack_goes: str,
) -> (typing.Dict, set):
    aspect_dict = dict(F="molecular_function",
                       P="biological_process",
                       C="cellular_component")
    result = {}
    all_go = set()
    e = 0
    number_seq = "?"
    
    #protein_list = [i for i in protein_list]
    begining = len(protein_list)
    protein_go_dict = {}
    while protein_list:
        protein_run = protein_list[:100]
        print(protein_run[0])
        protein_list = protein_list[100:]
        tries = 0
        url = f"https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId={','.join(protein_run)}"
        header = dict(Accept='text/tsv')

        logging.info(f"GO info downloaded for {protein_run} from {url} left {e}/{number_seq}")
        req = requests.get(url, headers=header, timeout=10)
        print(url, f"seq_no={e}", f"status_code={req.status_code}", f"tries={tries}")
        e += 1
        if req.status_code != 200:
            print(req, protein_run)
            print("err", req.status_code)
            break
        if not req.text.strip():
            print("lack of content", req, protein_run)
        for line in req.text.split("\n"):
            new_line = line.split("\t")
            if not line.startswith("GENE PRODUCT DB") and line.strip():
                if new_line[1] in protein_run:
                    protein_go = new_line[4]
                    aspect_go = new_line[5]
                    if aspect_go == aspect or aspect_go == aspect_dict.get(aspect):
                        annotation_type_go = new_line[7]
                        if annotation_type_go not in exclude:
                            all_go.add(protein_go)
                            if new_line[1] not in result:
                                result[new_line[1]] = [protein_go]
                            else:
                                result[new_line[1]].append(protein_go)
        for protein_acc in protein_run:
            if protein_acc in result.keys():
                #TO CLEAN
                # logging.info(f"Save go in {save_go_file} data:{str(set(result[protein_acc]))} ")
                #save_go(save_go_file, {protein_acc: set(result[protein_acc])}, mode="a")
                protein_go_dict[protein_acc] = result[protein_acc]
        #    else:
        #        with open(lack_goes, "a") as f:
        #            f.write(protein_acc + "\n")
            if e % 1000 == 0:
                logging.info(f"GO info taken from file for {protein_acc} left {e}/{number_seq}")
            if not result.get(protein_acc):
                result[protein_acc] = []
        print(len([1 for i, j in result.items() if len(j) > 0]), len(protein_list), begining,
              begining - len(protein_list))
    return result, all_go, protein_go_dict


def check_aspect(go, all_go, aspect):
    aspect_dict = dict(F="molecular_function",
                       P="biological_process",
                       C="cellular_component")
    # print(go, all_go)
    if go is None:
        return False
    if go in all_go:
        return True
    else:
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go}/"
        request = requests.get(url, timeout=10)
        request_json = request.json()
        if request_json.get("results", {}):
            aspect_go = [i for i in request_json["results"] if i["id"] == go][0]["aspect"]
            if aspect_go == aspect or aspect_go == aspect_dict.get(aspect):
                return True
    return False


def get_ancestors(
        go_list: set,
        #save_go_file: str,
        ancestors_old: dict,
        all_go: set,
        aspect: str
) -> (dict, set):
    ancestors = {}
    number_seq = len(go_list)
    for e, go in enumerate(go_list):
        if go not in ancestors_old:
            url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go.replace(':', '%3A')}/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
            print(url)
            print(f"GO ancestor info downloaded for {go} from {url} left {e}/{number_seq}")
            tries = 0
            success = False
            while tries < 10 and not success:
                try:
                    request = requests.get(url, timeout=10)
                    if request.status_code != 204:
                        request_json = request.json()
                        if request_json.get("results", {}):
                            ancestors[go] = [i.get("ancestors") for i in request_json.get("results", {}) if
                                             i["id"] == go and i.get("ancestors") is not None]
                            ancestors[go] = [i for sublist in ancestors[go] for i in sublist if
                                             i != go and check_aspect(i, all_go,
                                                                      aspect)]

                            all_go = all_go.union(set(ancestors[go]))
                            #save_go(save_go_file, {go: ancestors[go]}, "a")
                            success = True
                    else:
                        tries += 1
                except Exception as err:
                    print(err)
                    tries += 1

            if not success:
                logging.info(f"Lack of GO ancestor info for {go} left {e}/{number_seq}")
                ancestors[go] = []
    return ancestors, all_go


def save_go(file: str,
            go_set: dict,
            mode="a"
            ) -> None:
    stored_exception = None
    with open(file, mode) as f:
        for protein, goes in go_set.items():
            for go in goes:
                try:
                    f.write(f"{protein}\t{go}\n")
                except KeyboardInterrupt:
                    stored_exception = sys.exc_info()
    if stored_exception:
        raise (stored_exception[0], stored_exception[1], stored_exception[2])


def read_mapped_file(file, sign="\t"):
    result = {}
    if not os.path.exists(file):
        return result
    with open(file) as f:
        for line in f:
            if line.strip():
                line = line.split(sign)
                if line[0] not in result.keys():
                    result[line[0]] = [line[1].strip()]
                else:
                    if line[1].strip() not in result[line[0]]:
                        result[line[0]].append(line[1].strip())
    return result



def get_proteins(input_file_path):
    
    #for reading uniprot fasta file
    # with open(input_file_path) as f:
    #    for l in f.readlines():
    #        if l.startswith(">"):
    #            yield l.split("|")[1]

    with open(input_file_path, "r", encoding="utf-8") as f:
        proteins = [line.strip() for line in f]
    return proteins

def get_max_path(child: str, main_GO: str):
    url_path = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{child}/paths/{main_GO}/"
    req_path = requests.get(url_path)
    if req_path:
        max_path_len = 0
        for result_path in req_path.json()["results"]:
            if len(result_path) > max_path_len:
                max_path_len = len(result_path)
                return max_path_len
        return max_path_len


def get_paths(go: iter, path_path: str, aspect: str):
    aspect_dict = dict(F="GO:0003674",
                       molecular_function="GO:0003674",
                       P="GO:0008150",
                       biological_process="GO:0008150",
                       C="GO:0005575",
                       cellular_component="GO:0005575")
    with open(path_path, "w") as f:
        for go_id in go:
            path_len = get_max_path(go_id, aspect_dict[aspect])
            f.write(f"{go_id}\t{path_len}\n")

def add_ancestors(
        ancestors: dict,
        go_protein: dict
) -> dict:
    for e, protein in enumerate(go_protein.keys()):
        if e % 1000 == 0:
            logging.info(f"Add ancestors to {protein} {e}/{len(go_protein.keys())}")
        for old_go, anc in ancestors.items():
            if old_go in go_protein[protein]:
                go_protein[protein] += anc
                go_protein[protein] = list(set(go_protein[protein]))
    return go_protein


#old prepare_data function from go_analyse.py
def crate_annotation_file(all_go, ancestors, ouput_annotation_file):
    
    logging.info(f"Add ancestors info to protein GO")    
    all_go = add_ancestors(ancestors, all_go)


    with open(ouput_annotation_file, 'w', encoding='utf-8') as f:
        json.dump(all_go, f, indent=4, ensure_ascii=False)


def prepare_folders(input_file, exclude_IEA, ouptput_dir):

    #check if input file with proteins exists
    if not os.path.isfile(input_file):
        sys.exit("Input file with protein IDs does not exists. Exiting...")

    #create output directory    
    os.makedirs(ouptput_dir, exist_ok=True)
    
    #old files removal and creation of new empty files

    #go_max_path_file_path = os.path.join(ouptput_dir, GO_MAX_PATH_FILE)
    #if os.path.exists(go_max_path_file_path):
    #    os.remove(go_max_path_file_path)
    #path = Path(go_max_path_file_path)
    #path.touch()  

    go_names_file_path = os.path.join(ouptput_dir, GO_NAMES_FILE)
    if os.path.exists(go_names_file_path):
        os.remove(go_names_file_path)
    path = Path(go_names_file_path)
    path.touch()
    

    #go_proteins_file_path = os.path.join(ouptput_dir, GO_PROTEINS_FILE) 
    #if os.path.exists(go_proteins_file_path):
    #    os.remove(go_proteins_file_path)
    #path = Path(go_proteins_file_path)
    #path.touch()
    
    go_annotations_file_path = os.path.join(ouptput_dir, GO_ANNOTATIONS_FILE)

    if(exclude_IEA == "yes"):
        exclude_IEA = ["IEA"]
    else:
        exclude_IEA = []    


    return go_names_file_path, go_annotations_file_path, exclude_IEA

def main(options):
    
    [go_names_file_path, go_annotations_file_path, exclude_IEA] = \
        prepare_folders(options.input, options.exclude_IEA, options.output_dir)
    
    proteins = get_proteins(options.input)    

    proteins_go, all_go, protein_go_dict = get_GO(protein_list=proteins,
                                 exclude=exclude_IEA,
                                 aspect=options.aspect)
                                 #lack_goes=lack_go_file_path)

    ancestors, all_go = get_ancestors(set([item for sublist in list(proteins_go.values()) for item in sublist]),
                                      ancestors_old={},
                                      all_go=all_go,
                                      aspect=options.aspect)

    #create file with names of GO terms
    fill_names(all_go, save_file=go_names_file_path)

    #create file with max paths of GO terms
    #get_paths(all_go, path_path=go_max_path_file_path, aspect=options.aspect)

    #create json file with GO protein GO annotations including ancestors
    crate_annotation_file(protein_go_dict, ancestors, go_annotations_file_path)


def get_options():
    parser = OptionParser(description="desc")
    parser.add_option("-i", "--input", dest="input", default=None,
                      help="List of proteins for annotations", metavar="FASTA")
    parser.add_option("-e", "--exclude_IEA", dest="exclude_IEA", default="no",
                      help="Exclude GO terms with IEA? yes/no", metavar="STRING")
    parser.add_option("-s", "--aspect", dest="aspect", default="F",
                      help="Aspect of GO", metavar="STRING")
    parser.add_option('-o', '--output_dir', default='./gbsc_functional_results/', 
                      help='Project directory')
    options, args = parser.parse_args()

    return options, args


if __name__ == "__main__":
    # try:
    options, args = get_options()
    main(options)
