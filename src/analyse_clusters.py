"""
GBSC Clusters GO Ontology Functional Analysis Pipeline
=======================================================

Major refactoring of original analysis scripts by Joanna Ziemska-Legiecka (2025).
GO download logic, clusters GO enrichment algorithms and s-measure caluclations preserved with fixes.

Author: Aleksanra Gruca (2026)
Original: Joanna Ziemska-Legiecka (2025)
"""

def count_s_measure(cluster_sign_GO, cluster_size):
    return cluster_sign_GO / cluster_size

class AnaliseCluster:
    def __init__(self, enrichment_file, parameter_no):
        self.file_name = enrichment_file
        self.params_no = parameter_no
        self.clusters_info = {}
        self.c_value_cl = {}
        self.C_value = None
        self.cl_no = 0

    def count_c(self):
        for e, cl_name in enumerate(self.clusters_info.keys()):
            # print(cl_name, cl_name not in self.orphans, cl_name in self.clusters_main_GO)
            self.c_value_cl[cl_name] = count_s_measure(self.clusters_info[cl_name]["GO_sequences"],
                                                       self.clusters_info[cl_name]["cluster_size"])

    def read_enrichment_results(self):
        cl_names = set()
        with open(self.file_name) as f:
            next(f)
            for l in f:
                if l.strip():
                    line = l.strip().split()
                #"GBSC cluster\tGO ID\tp-value\tAll proteins\tAll proteins annotated with GO\
                #\tCluster size\tProteins annotated with GO in cluster\tBonferroni corrected p-value\t
                # \tBonferroni significance results, alpha={options.alpha}\tBenjamini-Hochberg corrected p-value\tBenjamini-Hochberg significance results, alpha={options.alpha}\n")                   
                    cluster, go, hypergeom_p_val, sequence_no, go_seq_no, cl_size, cluster_go_seq_no, bonferoni_correct, bonf_bool_test, BH_corrected, BH_bool_test = line
                    cl_size = int(cl_size)
                    if cl_size > 1 and BH_bool_test == "True":
                        cl_names.add(cluster)
                        cluster_go_seq_no=int(cluster_go_seq_no)
                        if cluster not in self.clusters_info:
                            self.clusters_info[cluster] = {"GO": go, "GO_sequences": cluster_go_seq_no,
                                                           "cluster_size": cl_size}
                        else:
                            if self.clusters_info[cluster]["GO_sequences"] < cluster_go_seq_no:
                                    self.clusters_info[cluster] = {"GO": go, "GO_sequences": cluster_go_seq_no,
                                                               "cluster_size": cl_size}
        self.cl_no = len(cl_names)

    def save(self, file, go_names_file_path):

        go_names_dict = {}    
        with open(go_names_file_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and not parts[0].startswith("#"):
                    go_names_dict[parts[0]] = parts[1]    

        with open(file, "w") as f:
            f.write(f"cluster_name;seq_no;s-measure;main GO ID;main GO name;\n")
            #for cluster, cluster_data in self.clusters_info.items():
            #    f.write(f"{cluster};{cluster_data['GO']};{cluster_data['cluster_size']};{self.c_value_cl[cluster]}\n")


            #reverse sort keys (cluster) based on value of self.c_value_cl
            sorted_clusters = sorted(
                self.clusters_info.keys(), 
                key=lambda cluster: self.c_value_cl[cluster], 
                reverse=True
            )

            for cluster in sorted_clusters:
                cluster_data = self.clusters_info[cluster]
                go_name = go_names_dict[cluster_data['GO']]
                f.write(f"{cluster};{cluster_data['cluster_size']};{self.c_value_cl[cluster]};{cluster_data['GO']};{go_name}\n")        

