import collections
import time
import csv
import sys
import pandas
from Bio import SeqIO

usage = """Usage:
organize_mapping.py -d target_mapping.txt -f target.fasta -o outfile.txt

Required Arguments:
  -d            depth output file from samtools depth with read depth per 
                basepair information
  -f            fasta files containing the target sequences reads were 
                mapped onto
  -o            output file

Output is written to STDOUT

This scripts creates a tab separated table with the position, read depth,
and nucleic acid for each target mapped onto to be used in the
sliding_window.py script.
"""

try:
    input = sys.argv[sys.argv.index('-d')+1]
except ValueError:
    print("ERROR: Requires '-d' depth input file path.\n\n")
    print(usage)
    quit()
try:
    fasta_input = sys.argv[sys.argv.index('-f')+1]
except ValueError:
    print("ERROR: Requires '-f' fasta input file path.\n\n")
    print(usage)
    quit()

#parse input
#create a dictionary that holds the sequence and per base read depth per position along the sequence
def parse_mapping( input):
    mapping = collections.defaultdict(dict)
    with open(input) as f:
        for line in f:
            if "#CHROM" in line:
                pass
            else:
                target,position,read_depth = line.replace('"','').strip().split("\t")
                if target in mapping.keys():
                    mapping[target]["target"].append(str(target))
                    mapping[target]["position"].append(int(position))
                    mapping[target]["read_depth"].append(int(read_depth))
                else:
                    mapping[target]["target"] = [str(target)]
                    mapping[target]["position"] = [int(position)]
                    mapping[target]["read_depth"] = [int(read_depth)]
    return(mapping)

def compute_stats(targ,fasta_seq):
    stats = []
    targ_array = targ["target"]
    pos_array = targ["position"]
    depth_array = targ["read_depth"]
    nuc_acid_array = []
    num_pos = len(pos_array)
    for i in range(len(pos_array)):
        nuc_acid_array.append(fasta_seq[targ_array[i]].seq[pos_array[i]-1]) 
    
    stats.append({"ID":targ_array, "position":pos_array, "read_depth":depth_array, "nucleic_acid":nuc_acid_array})
    
    return(stats)


#output table
#main()
if __name__ == '__main__':
    try:
        input = sys.argv[sys.argv.index('-d')+1]
    except ValueError:
        print("ERROR: Requires '-d' depth input file path.\n\n")
        print(usage)
        quit()
    try:
        fasta_input = sys.argv[sys.argv.index('-f')+1]
    except ValueError:
        print("ERROR: Requires '-f' fasta input file path.\n\n")
        print(usage)
        quit()
    try:
        output = sys.argv[sys.argv.index('-o')+1]
    except ValueError:
        print("ERROR: Requires '-o' output file path.\n\n")
    
    mapping = parse_mapping(input)
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(fasta_input),'fasta'))
    results = []
    
    for targ in mapping.keys():
        results += compute_stats(mapping[targ],fasta_sequences)

    #write tsv 
    tsv_columns = ["ID","position","read_depth", "nucleic_acid"]
    with open(output, 'w') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=tsv_columns, delimiter = '\t')
        writer.writeheader()
        for i in range(len(results)):
            row = collections.defaultdict(dict)
            for key in results[i].keys():
                for individual in range(len(results[i]["ID"])):
                    row[individual][str(key)] = results[i][str(key)][individual]
            for data in range(len(row)):
                writer.writerow(row[data])

