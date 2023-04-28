import argparse
from Bio import SeqIO
from sys import stdin as std_in
from sys import stdout as std_out

def mutate_sequence(fasta_file, out_file):
    """
    Mutates every base of each DNA sequence in the input FASTA file and outputs a FASTA of all mutated sequences.
    
    Args:
        fasta_file (str): The path to a FASTA file containing one or more DNA sequences.
        out_file (str): The path to save the output FASTA file.
    Returns:
        None.
    """
    mutations = ['A', 'C', 'G', 'T']
    
    # Use SeqIO to parse all sequences from the FASTA file
    records = SeqIO.parse(fasta_file, 'fasta')
    
    # Initialize FASTA file to store the mutated sequences
    fasta_file_output = open(f"{out_file}.fa", 'w')
        
    for record in records:
        input_seq = str(record.seq)
        
        # Save original sequence
        seq_header = f"{record.name}_ORIGINAL"
        fasta_file_output.write(">" + seq_header + "\n")
        fasta_file_output.write(input_seq + "\n")
        
        for i in range(len(input_seq)):
            for mutation in mutations:
                mutated_seq = input_seq[:i] + mutation + input_seq[i+1:]
                seq_header = f"{record.name}_{i}_{mutation}"
                fasta_file_output.write(">" + seq_header + "\n")
                fasta_file_output.write(mutated_seq + "\n")
    
    # Close the output FASTA file
    fasta_file_output.close()
    

def get_options():

    parser = argparse.ArgumentParser(description='ISM_FASTA')
    
    parser.add_argument('-I', metavar='input', nargs='?', default=std_in, required=True,
                        help='path to the input FASTA file, defaults to standard in')
    parser.add_argument('-O', metavar='output', nargs='?', default=std_out, required=True,
                        help='path to the output FASTA file, defaults to standard out')
    args = parser.parse_args()

    return args

def main():

    # Get options
    args = get_options()
    
    # Get FASTAs
    mutate_sequence(fasta_file=args.I, out_file=args.O)
    
if __name__ == '__main__':
    main()    
