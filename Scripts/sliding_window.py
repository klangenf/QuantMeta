import collections
import time
import csv
import sys

usage = """Usage:
sliding_windows_step2.py -m target_mapping.txt -n <int> -o outfile.txt

Required Arguments:
  -m            file with positional GC content and read depth per
                per basepair information
                output from 'organize_mapping.R'
  -o            output file path

Optional Arguments:
  -n            integer value for the sliding window sizes (default is 49)


Output is written to STDOUT

Nucleotide format: ATGC (non-ATGCs will be ignored) or 0 for A/T and 1 for G/C
    Non-ATGCs will be ignored from percent GC calculations
    Nucleotide cases are ignored
This scripts creates sliding windows with the average GC and read depth for
each sliding window based on the 'organize_mapping.R' output.
"""


try:
    input = sys.argv[sys.argv.index('-m')+1]
except ValueError:
    print("ERROR: Requires '-m' input file path.\n\n")
    print(usage)
    quit()
try:
    window_size = int(sys.argv[sys.argv.index('-n')+1])
except ValueError:
    window_size = 49


#parse input
#create a dictionary that holds the sequence and per base read depth per contig
def parse_input( input):
    seqs = collections.defaultdict(dict)
    with open(input) as f:
        for line in f:
            if "position" in line:
                pass
            else:
                contig_id,position,read_depth,nuc_acid = line.replace('"','').strip().split("\t")
                if contig_id in seqs.keys():
                    seqs[contig_id]["sequence"].append(nuc_acid)
                    seqs[contig_id]["read_depth"].append(int(read_depth))            
                else:
                    seqs[contig_id]["sequence"] = [nuc_acid]
                    seqs[contig_id]["read_depth"] = [int(read_depth)]
    return(seqs)

def compute_stats(seq,window_size,contig):
    stats = []
    seq_array = seq["sequence"]
    depth_array = seq["read_depth"]
    seq_len = len(seq_array)

    # check if seq is greater than sliding window
    if seq_len < window_size:
        print("does not meet length req")
        #return -1

    # get remaining windows and stats
    for i in range(seq_len - window_size):
        end = i + window_size
        seq_window = seq_array[i:end]
        depth_window = depth_array[i:end]
        avg_depth = compute_avg_depth(depth_window)
        gc = compute_gc(seq_window)
        stats.append({"ID":contig,"avg_gc":gc,"avg_depth":avg_depth})
    
    return(stats)
    

def compute_gc(seq_array):
    return((seq_array.count("C") + seq_array.count("G")) / len(seq_array))

def compute_avg_depth(depth_array):
    return(sum(depth_array)/len(depth_array))

#output table
#main()
if __name__ == '__main__':
    try:
        input = sys.argv[sys.argv.index('-m')+1]
    except ValueError:
        print("ERROR: Requires '-m' input file path.\n\n")
        print(usage)
        quit()
    try:
        window_size = int(sys.argv[sys.argv.index('-n')+1])
    except ValueError:
        window_size = 49
    try:
        output = sys.argv[sys.argv.index('-o')+1]
    except ValueError:
        print("ERROR: Requires '-o' output file path.\n\n")
    seqs = parse_input(input)
    stats = []
    
    for seq in seqs.keys():
       stats += compute_stats(seqs[seq],window_size,seq)

    #write tsv 
    tsv_columns = ["ID","avg_gc","avg_depth"]
    with open(output, 'w') as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=tsv_columns, delimiter = '\t')
        writer.writeheader()
        for row in stats:
            writer.writerow(row)
 
