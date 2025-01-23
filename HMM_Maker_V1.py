import os
import argparse
import subprocess
import csv
from Bio import SeqIO
from collections import defaultdict

def extract_sequences(genbank_file):
    gene_sequences = defaultdict(list)
    gene_products = {}
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                gene = feature.qualifiers.get("gene", ["Unknown"])[0]
                # Skip genes named "Unknown"
                if gene == "Unknown":
                    continue
                product = feature.qualifiers.get("product", [""])[0]
                if feature.qualifiers.get("translation"):
                    aa_seq = feature.qualifiers["translation"][0]
                    gene_sequences[gene].append(aa_seq)
                    gene_products[gene] = product
    return gene_sequences, gene_products

def create_fasta_files(gene_sequences, output_dir):
    for gene_name, sequences in gene_sequences.items():
        with open(os.path.join(output_dir, f"{gene_name}.fasta"), "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">{gene_name}_{i+1}\n{seq}\n")

def align_sequences(fasta_files, output_dir):
    aligned_files = []
    for fasta_file in fasta_files:
        gene_name = os.path.splitext(os.path.basename(fasta_file))[0]
        
        # Check number of sequences in the FASTA file
        with open(fasta_file, 'r') as f:
            sequence_count = sum(1 for line in f if line.startswith('>'))
        
        if sequence_count > 1:
            # Multiple sequences, perform alignment
            try:
                subprocess.run(["muscle", "-align", fasta_file, "-output", fasta_file], 
                               check=True, 
                               stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                print(f"Alignment failed for {fasta_file}")
                print(f"STDOUT: {e.stdout.decode()}")
                print(f"STDERR: {e.stderr.decode()}")
        
        aligned_files.append(fasta_file)
    return aligned_files

def build_hmm_profiles(fasta_files, output_dir):
    hmm_dir = os.path.join(output_dir, "hmms")
    os.makedirs(hmm_dir, exist_ok=True)
    for fasta_file in fasta_files:
        # Check sequence length
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        if not sequences or len(sequences[0].seq) < 50:  # Check if file is empty or sequence is too short
            print(f"Skipping HMM build for {fasta_file}: sequence too short (< 50 amino acids)")
            continue
        
        gene_name = os.path.splitext(os.path.basename(fasta_file))[0]
        hmm_file = os.path.join(hmm_dir, f"{gene_name}.hmm")
        try:
            subprocess.run(["hmmbuild", "--amino", hmm_file, fasta_file], check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"Error building HMM for {fasta_file}: {e.stderr}")

def create_metadata_sheet(gene_products, output_dir, metadata_filename):
    metadata_file = os.path.join(output_dir, metadata_filename)
    with open(metadata_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['ID', 'Function', 'Gene'])
        for gene_name, product in gene_products.items():
            writer.writerow([gene_name, product, gene_name])

def concatenate_hmm_profiles(output_dir, final_hmm_name):
    hmm_dir = os.path.join(output_dir, "hmms")
    hmm_files = [f for f in os.listdir(hmm_dir) if f.endswith('.hmm')]
    final_hmm = os.path.join(output_dir, final_hmm_name)
    with open(final_hmm, 'wb') as outfile:
        for hmm_file in hmm_files:
            with open(os.path.join(hmm_dir, hmm_file), 'rb') as infile:
                outfile.write(infile.read())
    return final_hmm_name

def process_genbank_files(input_dir, output_dir):
    all_gene_sequences = defaultdict(list)
    all_gene_products = {}
    
    for filename in os.listdir(input_dir):
        if filename.endswith(".gb") or filename.endswith(".gbk"):
            genbank_file = os.path.join(input_dir, filename)
            gene_sequences, gene_products = extract_sequences(genbank_file)
            
            for gene_name, sequences in gene_sequences.items():
                all_gene_sequences[gene_name].extend(sequences)
            
            all_gene_products.update(gene_products)
    
    fasta_dir = os.path.join(output_dir, "fasta")
    os.makedirs(fasta_dir, exist_ok=True)
    create_fasta_files(all_gene_sequences, fasta_dir)
    fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith('.fasta')]
    
    aligned_files = align_sequences(fasta_files, fasta_dir)
    build_hmm_profiles(aligned_files, output_dir)
    
    final_hmm_name = "MiBiG_4.hmm"
    concatenate_hmm_profiles(output_dir, final_hmm_name)
    
    metadata_filename = os.path.splitext(final_hmm_name)[0] + ".tsv"
    create_metadata_sheet(all_gene_products, output_dir, metadata_filename)

def main():
    parser = argparse.ArgumentParser(description="Build HMM profiles and metadata sheet from GenBank files")
    parser.add_argument("input_dir", help="Directory containing GenBank files")
    parser.add_argument("output_dir", help="Directory to store output files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    process_genbank_files(args.input_dir, args.output_dir)
    print("Process completed. HMM profiles, metadata sheet, and concatenated HMM file have been created.")

if __name__ == "__main__":
    main()
