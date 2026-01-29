"""
GBSC Clusters GO Ontology Functional Analysis Pipeline
=======================================================

Major refactoring of original analysis scripts by Joanna Ziemska-Legiecka (2025).
GO download logic, clusters GO enrichment algorithms and s-measure caluclations preserved with fixes.

Author: Aleksandra Gruca (2026)
Original: Joanna Ziemska-Legiecka (2025)
"""

GO_ANNOTATIONS_FILE = "go_annotations.json"
ENRICHMENT_RESULTS_FILE="enrichment_results.csv"
GO_NAMES_FILE ="go_names.csv"


import os
import sys
import logging
from src.analyse_clusters import AnaliseCluster
from optparse import OptionParser
from src.go_analise import run_go_analyse
from src.utils import get_all_gbsc_proteins
from pathlib import Path


logger = logging.getLogger()
logger.setLevel(logging.INFO)


def check_folders(options):
    
    ouput_dir = options.ouput_dir
    if not os.path.isdir(ouput_dir):
        sys.exit(f"Exiting....\n Project output directory {ouput_dir} does not exist. \
            \n Run 01_download_go.py script first or check the path to project directory")
    
    gbsc_clusters_path = options.gbsc_clusters
    #check if GSBSC clusters directory path exists
    if not os.path.isdir(gbsc_clusters_path):
        sys.exit("Path to the directory with GBSC clusters does not exists. Exiting...")
    
    #check if GSBSC clusters directory path is not empty
    if not os.listdir(gbsc_clusters_path):
        sys.exit("Directory with GBSC clusters is empty. Exiting...")

    #dictuionary file with mapping GO IDs to GO names
    go_names_file_path = os.path.join(ouput_dir, GO_NAMES_FILE)
    if not os.path.exists(go_names_file_path):
        sys.exit(f"Exiting....\n File {go_names_file_path} does not exist. \
                 \n Run 01_download_go.py script first or check the path to project directory")

    #file with GO annotations - result of 01_download_go.py
    go_annotations_file_path = os.path.join(ouput_dir, GO_ANNOTATIONS_FILE)
    if not os.path.exists(go_annotations_file_path):
         sys.exit(f"Exiting....\n File {go_annotations_file_path} does not exist. \
                 \n Run 01_download_go.py script first or check the path to project directory")  


    enrichment_results_file_path = os.path.join(ouput_dir, ENRICHMENT_RESULTS_FILE)
    with open(enrichment_results_file_path, "w") as f:
        f.write(f"GBSC cluster\tGO ID\tp-value\tAll proteins\tAll proteins annotated with GO\
                \tCluster size\tProteins annotated with GO in cluster\tBonferroni corrected p-value\tBonferroni significance results alpha={options.alpha}\tBenjamini-Hochberg corrected p-value\tBenjamini-Hochberg significance results alpha={options.alpha}\n")
    
    return go_annotations_file_path, gbsc_clusters_path, go_names_file_path, enrichment_results_file_path

def main(options, args):
    
    file_handler = logging.FileHandler(options.log_file, mode='w', encoding='utf-8')
    logger.addHandler(file_handler)

    #check if requierd files with GO annotations exists
    [go_annotations_file_path, gbsc_clusters_path, go_names_file_path, enrichment_results_file_path] = check_folders(options)

    #helper function to get all GBSC protein IDs for test set for 01_download_go.py  
    # get_all_gbsc_proteins(gbsc_clusters_path)

    run_go_analyse(enrichment_results_file_path, go_annotations_file_path, options.alpha, gbsc_clusters_path)

    print(f"Estimation results saved to {enrichment_results_file_path}")     
    
    analyse_cluster = AnaliseCluster(enrichment_results_file_path, "0")
    analyse_cluster.read_enrichment_results()
    analyse_cluster.count_c()

    s_values_file = os.path.join(options.ouput_dir, "clusters_s_values.txt") 
    analyse_cluster.save(s_values_file, go_names_file_path)

    print(f"s-measure results for GBSC clusters saved to {s_values_file}")     

    
def get_options():
    parser = OptionParser(description="desc")
    parser.add_option("-c", "--gbsc_clusters", dest="gbsc_clusters", default=None,
                      help="Path to the directory with GBSC clusters", metavar="DIR")
    parser.add_option("-a", "--alpha", dest="alpha", default=0.05,
                      help="Threshold of test significance", metavar="FLOAT")
    parser.add_option('-o', '--ouput_dir', default='./gbsc_functional_results/', 
                      help='Output directory. Should be the same as used in 01_download_go.py')
    parser.add_option('-l', '--log_file', default='gbsc_functional_analysis.log', 
                      help='Log file name')
    options, args = parser.parse_args()
    return  options, args


if __name__ == "__main__":
    # try:
    options, args = get_options()
    main(options, args)