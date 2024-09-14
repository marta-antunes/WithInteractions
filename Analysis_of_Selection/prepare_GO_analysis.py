import requests
import Bio
from Bio import SeqIO

def read_gff(file_path):
    """
    Read GFF file and return a dictionary with CDS entries.
    """
    cds_info = {}
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == "CDS":  # Filter by CDS
                    gene_id = fields[8].split(';')[0].split('=')[1]
                    cds_info[gene_id] = fields
    return cds_info

def get_proteinName_for_selectedGenes(gff, subset_file):
    """
    Subset GFF information based on a list of gene IDs.
    """
    gff_values = gff.values()
    subset_info = {}
    with open(subset_file, 'r') as file:
        subset_genes = [line.strip() for line in file]
        for gene_id in subset_genes:
            for item in gff_values:  
                print(item)
                if any(gene_id in sub_item for sub_item in item):
                    for sub_item in item:
                        if "product=" in sub_item:
                            product_info = sub_item.split("protein_id=")[1]
                            subset_info[gene_id] = product_info
                            break 
    return subset_info

def get_unique_in_dict(loc_dict):
    printed_locs = set()

    for loc, value in loc_dict.items():
        if loc not in printed_locs:
            print(f'{loc}: {value}')
            #printed_locs.add(loc)


# Example usage:
# Read GFF file
gff_info = read_gff("GCF_008121235.1_UCBerk_Dsub_1.0_genomic.gff")

# Subset GFF file based on another file
subset_info = get_proteinName_for_selectedGenes(gff_info, "CandidatesLowLat.txt")

#print unique in dict
get_unique_in_dict(subset_info)
