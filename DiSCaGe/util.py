import os
import networkx as nx
import pandas as pd
import collections
import operator


'''
Input: Binary Mutation Matrix (BMM) or weight Mutation Matrix (WMM). Both are dataframes
Output: Mutation Matrix in a text file
'''
def create_txt_file_mutation_matrix(matrix, output_file_name):
    matrix.to_csv(output_file_name, sep='\t')


def get_mutation_matrix_from_file(mm_file_name):
    mm = pd.read_csv(mm_file_name, sep='\t', index_col=0)
    return mm


def get_mutation_score_from_file(ms_file_name):
    gene_scores_mutations = pd.read_csv(ms_file_name, sep='\t', index_col=0, squeeze=True).to_dict()
    return gene_scores_mutations


'''
Input: Alteration matrix
Output: A simple report about the samples
'''
def get_sample_report(am, output_file_name):
    output_file = open(output_file_name, 'w')
    
    num_samples = len(am)
    output_file.write("# Number of samples: {}\n".format(num_samples))
    output_file.write("sample\tmutated_genes\n")
    for sample in am:
        output_file.write("{}\t{}\n".format(sample, len(am[sample])))
        
    output_file.close()
    return 0

'''
Input: Binary Mutation matrix
Output: Output: A simple report about the genes
'''
def get_gene_report(bmm, output_file_name):
    output_file = open(output_file_name, 'w')
    
    genes = list(bmm)
    output_file.write("# Number of genes: {}\n".format(len(genes)))
    output_file.write("gene\tmutated_samples\tfreq\n")
    for gene in genes:
        q = bmm[gene].sum()
        f = round(q / len(list(bmm.index)), 3)
        output_file.write("{}\t{}\t{}\n".format(gene, q, f))
    
    output_file.close()
    return 0


'''
Input: Binary Mutation matrix
Output: Output: A simple report about the genes
'''
def generate_gene_score_mutation_file(gene_scores_mutations, output_file_name):
    output_file = open(output_file_name, 'w')
    
    sorted_gene_scores_mutations = sorted(gene_scores_mutations.items(), key=lambda kv: kv[1], reverse=True)
    output_file.write("gene\tmutation_score\n")
    for gene, score in sorted_gene_scores_mutations:
        output_file.write("{}\t{}\n".format(gene, score))
    output_file.close()
    return 0

'''
Input: Related Gene Network
Output: Information about the RGN
'''
def get_network_report(rgn, output_file_name):
    output_file = open(output_file_name, 'w')
    output_file.write("Number of nodes (genes) in the network:\t{}\n".format(len(list(rgn.nodes))))
    output_file.write("Number of edges (interactions) in the network:\t{}\n".format(len(list(rgn.edges))))
    output_file.close()
    return 0

   
'''
Input: output folder and the name of each mutation matrix"
Output: true if all files exists
'''
def files_exists(folder, wmm, bmm, am):
    if os.path.isfile(folder + "/" + wmm) and os.path.isfile(folder + "/" + bmm) and os.path.isfile(folder + "/" + am):
        return True
    return False
    
