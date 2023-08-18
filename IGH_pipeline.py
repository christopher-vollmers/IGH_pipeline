import sys
import os

fofn=sys.argv[1]

for line in open(fofn):
    path=line.strip()
    print path
    os.system('python index_reads.py '+ path)
    os.system('python fuzzy_index_comparison.py '+path)
    os.system('python determine_consensus_simplified_retain.py '+path)
    os.system('python combine_identicals_cluster.py '+ path)
