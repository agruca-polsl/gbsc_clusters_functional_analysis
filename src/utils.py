"""
GBSC Clusters GO Ontology Functional Analysis Pipeline
=======================================================

Major refactoring of original analysis scripts by Joanna Ziemska-Legiecka (2025).
GO download logic, clusters GO enrichment algorithms and s-measure caluclations preserved with fixes.

Author: Aleksandra Gruca (2026)
Original: Joanna Ziemska-Legiecka (2025)
"""

from typing import IO
from typing import Iterator
import os

class Protein:
    def __init__(self, header="", sequence=""):
        self.header = header
        self.sequence = sequence

    def get_acc(self):
        return self.header.split("|")[1]

    def is_sequence_pure(self):
        global canonical_residues
        return len(self.sequence) >= 20 and set(self.sequence) <= canonical_residues


def get_proteins(db_file: IO) -> Iterator[Protein]:
    protein = None
    while True:
        line = db_file.readline()
        if not line:
            if protein is not None:
                yield protein
            break
        if line.startswith(">"):
            if protein is not None:
                yield protein
            protein = Protein(line.strip())
        elif protein is not None:
            protein.sequence += line.strip()


def get_all_gbsc_proteins(gbsc_clusters_path):
    files = os.listdir(gbsc_clusters_path)
    all_gbsc_proteins = []
    for cluster_file in files:
        cluster_file_path = os.path.join(gbsc_clusters_path, cluster_file)
        proteins = get_proteins(open(cluster_file_path))
        all_gbsc_proteins.extend([i.get_acc() for i in proteins])

    all_gbsc_proteins_uniqe = list(set(all_gbsc_proteins))
    with open('gbsc_test_protein_IDs.txt', 'w') as f:
        f.write('\n'.join(all_gbsc_proteins_uniqe))
    