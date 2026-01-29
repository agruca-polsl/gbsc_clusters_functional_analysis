"""
GBSC Clusters GO Ontology Functional Analysis Pipeline
=======================================================

Major refactoring of original analysis scripts by Joanna Ziemska-Legiecka (2025).
GO download logic, clusters GO enrichment algorithms and s-measure caluclations preserved with fixes.

Author: Aleksandra Gruca (2026)
Original: Joanna Ziemska-Legiecka (2025)
"""

import os
import sys
import json
import logging
from multiprocessing import Pool


from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

from src.utils import get_proteins


logger = logging.getLogger(__name__)


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




def calc_hypergeometric_test(
        cluster_dict: dict,
        full_set_dict: dict,
        file: str
) -> dict:
    go_res = {}
    all_go = []
    for prot, goes in cluster_dict.items():
        all_go += list(goes)
    all_go = list(set(all_go))
    for go in all_go:
        go_res[go] = calc_hypergeometric_single_go_test(
            cluster_dict,
            full_set_dict,
            go,
            file
        )
    return go_res


def calc_hypergeometric_single_go_test(
        cluster: dict,
        full_set: dict,
        go: str,
        file: str
) -> tuple:
    # https: // alexlenail.medium.com / understanding - and -implementing - the - hypergeometric - test - in -python - a7db688a7458
    M = len(set(full_set.keys()))  # liczba wszystkich białek
    m = get_number_of_all_proteins_for_go(full_set, go)  # liczba wszystkich białek z badanym GO
    N = len(set(cluster.keys()))  # rozmiar klastra
    x = get_number_of_all_proteins_for_go(cluster, go)  # liczba białek w klastrze z badanym GO
    stat = hypergeom.sf(x - 1, M, m, N)
    logging.info(f"file: {file}, parameters for GO {go}: M={M} m={m} k={N} x={x} stat={stat}")
    return stat, M, m, N, x


def get_number_of_all_proteins_for_go(proteins_go: dict, go: str) -> int:
    return len(set({i: j for i, j in proteins_go.items() if go in j}))


def save_results(
        output_file: str,
        file: str,
        result: dict,
        bonf_correction: float,
        bh: dict,

) -> None:
    # print(correction)
    with open(output_file, "a") as f:
        if result:
            for go, value in result.items():
                bonf_correction_test_result = value[0] < bonf_correction
                f.write(
                    # f"{file}\t{go}\t{value[0]}\t{value[1]}\t{value[2]}\t{value[3]}\t{value[4]}\t{correction}\t{value[0] < correction}\t{bh[go]}\t{goes_nr}\t{value[0] < bh[go]}\t{value[0] < correction and value[0] < bh[go]}\n")
                    #Asia's version
                    #f"{file}\t{go}\t{value[0]}\t{value[1]}\t{value[2]}\t{value[3]}\t{value[4]}\t{correction}\t{correction <= alpha}\t{bh[go]}\t{goes_nr}\t{value[0] <= alpha}\n")
                    #removed goes_nr - number of GO in cluster
                    f"{file}\t{go}\t{value[0]}\t{value[1]}\t{value[2]}\t{value[3]}\t{value[4]}\t{bonf_correction}\t{bonf_correction_test_result}\t{bh[go][0]}\t{bh[go][1]}\n")


# nazwa pliku;Go name;

def select_go_for_cluster(
        cluster_proteins: list,
        all_go: dict
) -> dict:
    return {i: j for i, j in all_go.items() if i in cluster_proteins}


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


def calc_bonferroni_correction(alpha: float, go_clusters: dict) -> tuple:
    goes = []
    for go in go_clusters.values():
        goes += list(go)
    if len(set(goes)) != 0:
        return alpha / len(set(goes)), len(set(goes))
    return 100, len(set(goes))


def Benjamini_Hochberg(pvals, significance):
    pvals_test = list(pvals.values())
    if pvals_test:
        # print(pvals_test)
        rest = multipletests(pvals=[i[0] for i in pvals_test], alpha=significance, method="fdr_bh")
        # print("rest", rest)
        new_pvals = {i[0]: (rest[1][e], rest[0][e]) for e, i in enumerate(pvals.items())}
        return new_pvals
    else:
        return {}


def run_go_analyse(output_file, go_annotations_file, alpha, folder_clusters):

    files = os.listdir(folder_clusters)
    #print(f"Proteins {len([i for i, j in all_go.items() if not j])} do not have GO")
    
    #read GO annotations including ancestors for all proteins    
    with open(go_annotations_file, 'r', encoding='utf-8') as f:
        all_go = json.load(f)

    result_dict = {}
    len_files = len(files)
    runs = []
    for e, file in enumerate(files):
        runs.append((file, e, len_files, folder_clusters, all_go, alpha, output_file))
    
    #with Pool(100) as p:
    #    results = p.map(calc, runs)

    results = map(calc, runs)

    for result in results:
        data_go_results, file, test, bonf_correction, bh = result

        if output_file:
        #    logging.info(str(result))
        #    logging.info(str(test))
            save_results(output_file, file, test, bonf_correction, bh)
        result_dict[file] = data_go_results
        # print(f"{all_cl} clusters do not have any protein with GO.")
    return result_dict


def calc(data):
    file, e, len_files, folder_clusters, all_go, alfa, output_file = data
    logging.info(f"Starting calculations for {file} {e}/{len_files}")
    cluster_file = os.path.join(folder_clusters, file)
    cluster = get_proteins(open(cluster_file))
    logging.info('Selecting GO info for protein cluster')
    cluster = [i.get_acc() for i in cluster]
    cluster_go = select_go_for_cluster(cluster, all_go)
    data_go_results = []
    test, bonf_correction, bh, goes_nr = None, None, None, None
    if cluster_go:
        # print(cluster_go)
        logging.info("Running hypergeometric test")
        test = calc_hypergeometric_test(cluster_go, all_go, file)
        logging.info("Running bonferroni correction")
        bonf_correction, goes_nr = calc_bonferroni_correction(alfa, cluster_go)
        logging.info("Running Benjamini-Hochberg corection")
        bh = Benjamini_Hochberg(test, alfa)
        logging.info("Writing results to file")
        for go, value in test.items():
            data_go = dict(
                go=go,
                pvalue_hypergeom=value[0],
                all_proteins=value[1],
                all_proteins_with_go=value[2],
                cluster_size=value[3],
                proteins_with_go_in_cluster=value[4],
                bonf_correction=bonf_correction,
                result_test_bonf_adj=value[0] < bonf_correction,
                #to save info the bh dictionary us used - below just for clarity what is in bh[0] and bh[1]
                #pval_bh_adj=bh[go][0], - adjusted p-value
                #test_result_bh_adj= bh[go][1], - information if
                #goes_nr_in_cluster=goes_nr, 
                #sinificance_test=value[0] < bonf_correction and value[0] < bh[go],
            )
            data_go_results.append(data_go)
            # if file not in result.keys():
            #     result[file] = [data_go]
            # else:
            #     result[file].append(data_go)
        if not data_go_results:
            data_go_results.append([dict(
                cluster=file,
                sinificance_test=None,
                go=None,
            )])
    else:
        logging.info(f"No GO for cluster {file}")
        if not data_go_results:
            data_go_results.append(dict(
                cluster=file,
                sinificance_test=None,
                go=None,
            ))

    return data_go_results, file, test, bonf_correction, bh


