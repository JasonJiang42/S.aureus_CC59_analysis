import os
import subprocess
import argparse

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Modify FASTA alignment files")
parser.add_argument("-i", "--input", required=True, help="Input file")
parser.add_argument("-o", "--output", required=True, help="Output file")
args = parser.parse_args()

# Validate input file
input_file = args.input
if not os.path.isfile(input_file):
    print(f"Error: {input_file} is not a file")
    exit(1)

# Check sequence lengths using seqkit fx2tab
seq_lengths = {}
output = subprocess.check_output(["seqkit", "fx2tab", "--name", "--length", input_file])
for line in output.decode().splitlines():
    name, length = line.split()
    seq_lengths[name] = int(length)

# Check if all sequence lengths are equal
if len(set(seq_lengths.values())) != 1:
    print("Error: sequence lengths are not equal")
    exit(1)

# Process the input file
output_file = args.output
with open(input_file) as in_file, open(output_file, "w") as out_file:
    seq_count = len(seq_lengths)
    seq_length = list(seq_lengths.values())[0]
    out_file.write(f"{seq_count}    {seq_length}\n")
    for line in in_file:
        if line.startswith(">"):
            out_file.write(line[1:])
        else:
            out_file.write(line)
