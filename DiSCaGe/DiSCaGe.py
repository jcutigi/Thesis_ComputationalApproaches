import datetime
import os
import sys
import configparser
import pandas as pd
import networkx as nx
import operator as op
import util
import logging

# constants
NONSENSE = "Nonsense_Mutation"
MISSENSE = "Missense_Mutation"
SPLICE ="Splice_Site"
FRAMESHIFT_DEL = "Frame_Shift_Del"
FRAMESHIFT_INS = "Frame_Shift_Ins"
INFRAME_DEL = "In_Frame_Del"
INFRAME_INS = "In_Frame_Ins"
_3UTR = "3'UTR"
_5UTR = "5'UTR"
NONSTOP_MUTATION = "Nonstop_Mutation"
TRANSLATION_START_SITE = "Translation_Start_Site"

'''
Input: Pandas DataFrame with the MAF file information
OutPut: Binary Mutation Matrix (BMM), in a pandas DataFrame

Binary Mutation matrix format:
Rows are patients ((samples) and columns are genes.
1: gene is mutated in patient; 0: gene is not mutated in patient.

Example:
Hugo_Symbol g1  g2  g3  ... gm
p1  0   0   1   0
p2  1   0   0   1
p3  1   0   1   0
...
pn  1   0   0   0

'''
def get_binary_mutation_matrix(maf):
    bmm = pd.crosstab(maf.Tumor_Sample_Barcode, maf.Hugo_Symbol).clip(upper=1)
    return bmm


'''
Input: Pandas DataFrame with the MAF file information
Output: Weighted Mutation Matrix (WMM), in a pandas DataFrame
'''
def get_weighted_mutation_matrix(maf, mutation_weights):
    patients = list(set(maf.Tumor_Sample_Barcode))
    genes = list(set(maf.Hugo_Symbol))
    mutation_types = list(set(maf.Variant_Classification))
    
    mutation_weights = {k: v for k, v in mutation_weights.items() if k in mutation_types}
    
    # cross tab to get the number of mutations of each pair gene-patient
    mut = pd.crosstab(maf.Hugo_Symbol, [maf.Tumor_Sample_Barcode, maf.Variant_Classification], values=maf.Variant_Classification, aggfunc=len, dropna=False)

    wmm = pd.DataFrame(index=patients, columns=genes) # Creating Weighted Mutation Matrix to keep the scores.
    wmm = wmm.fillna(0.0) # fill with 0s rather than NaNs
    
    # calcuting the scores for each pait gene-patient: s(p_i, g_j)
    for p in patients:
        mut[p] = mut[p].fillna(0) # fill with 0s rather NaNs
        for g in genes:
            mut_dict = mut[p].loc[g].to_dict() # getting dictionary of number of mutations for the pair patient-gene, organized by mutation type
            
            number_of_mutations = sum(mut_dict.values())
            sum_product = sum(mut_dict[k] * mutation_weights[k] for k in mut_dict)
            score = 0
            if number_of_mutations > 0:
                score =  round(sum_product / number_of_mutations, 3)
            wmm.at[p, g] = score
    return wmm



'''
Input: Binary Mutation Matrix (BMM)
Output: An Alteration Matrix, in a dictionary

Each row contais patient (sample) in the first columns. The follow columns have the genes which are mutated in the sample.
Example:
p1  g1  g4  g13
p2  g7  g16 g53 g104
p3  g10 g12
...
pn  g1  g3
'''
def get_alteration_matrix(bmm):
    matrix = dict()
    sample_list = bmm.index.tolist()
    for sample in sample_list:
        matrix[sample] = list(bmm.columns.values[bmm.loc[sample] == 1])

    return matrix


def get_genes_from_bmm(bmm):
    return list(bmm)

'''
Input: Binary Mutation Matrix (BMM) and Weighted Mutation Matrix (WMM)
Output: A score for each gene
'''
def get_score_genes(bmm, wmm):
    samples = bmm.index.tolist()
    genes = list(bmm)
    gene_scores = dict()
    
    score_max = 0
    coeff = 0
    for g in genes:
        score = round(wmm[g].sum() / len(samples), 3) # getting the 'weighted frequence' of the scores in the gene
        gene_scores[g] = score
        if score > score_max: # keeping the biggest score
            score_max = score
    
    # putting scores between 0 and 1
    for g in gene_scores:
        gene_scores[g] = round(gene_scores[g] / score_max, 3)
      
    return gene_scores


'''
Input: the gene network file (reading using rb)
Output: a networkX graph
'''
def read_gene_network_nx(gene_network_file):
    mx = nx.read_edgelist(gene_network_file, delimiter='\t')
    return mx

'''
Input: a gene network (networkX), a threshold and a a list of genes considered in the analisys
Output: a Gene Strength Spreading Network (GSSN) and the biggest weight considering all the edges
'''
def create_gene_strength_spreading_network(gene_network_nx):
    gssn = nx.DiGraph()
    gene_network_nx.remove_edges_from(nx.selfloop_edges(gene_network_nx))
    
    
    max_weight = 0
    for g_i, g_j in gene_network_nx.edges:
        neighbors_i = dict(gene_network_nx[g_i])
        neighbors_j = dict(gene_network_nx[g_j])
        
        neighbors_j_out = {k:v for k,v in neighbors_j.items() if k not in neighbors_i}
        neighbors_j_out.pop(g_i, None)
        r_i = sum([neighbors_i[g]["weight"] for g in neighbors_i])
        r_j_out = sum([neighbors_j_out[g]["weight"] for g in neighbors_j_out])


        neighbors_i_out = {k:v for k,v in neighbors_i.items() if k not in neighbors_j}
        neighbors_i_out.pop(g_j, None)
        r_j = sum([neighbors_j[g]["weight"] for g in neighbors_j])
        r_i_out = sum([neighbors_i_out[g]["weight"] for g in neighbors_i_out])

        weight_ij = gene_network_nx[g_i][g_j]["weight"]
        s_ij = (1 + (r_i * r_j_out)) * weight_ij
        s_ji = (1 + (r_j * r_i_out)) * weight_ij

        gssn.add_edge(g_i, g_j, weight=s_ij)
        gssn.add_edge(g_j, g_i, weight=s_ji)

        if max(s_ij, s_ji) > max_weight:
            max_weight = max(s_ij, s_ji)

    return gssn, max_weight


def combine_gene_network(gene_network_nx_list):
    consensus_gene_network = nx.Graph()
    
    max_weight = 1
    for gene_network in gene_network_nx_list:
        for g_i, g_j in gene_network.edges:
            if not consensus_gene_network.has_edge(g_i, g_j):
                consensus_gene_network.add_edge(g_i, g_j, weight=1)
            else:
                consensus_gene_network[g_i][g_j]["weight"] += 1
                if consensus_gene_network[g_i][g_j]["weight"] > max_weight:
                    max_weight = consensus_gene_network[g_i][g_j]["weight"]
    
    for g_i, g_j in consensus_gene_network.edges:
        consensus_gene_network[g_i][g_j]["weight"] /= max_weight

    return consensus_gene_network


'''
Input: gene scores, a directed and weighted gene network (networkX), the biggest weight considering all the edges of network
Output: a Gene Strength Spreading Network (GSSN) and the biggest weight considering all the edges
'''
def get_score_genes_neighbors(gene_scores, gssn, max_weight, consensus_gene_network):
    gene_scores_from_neighbors = {}
   
    genes_out = list(set(gssn) - set(gene_scores)) # genes of gssn not in gene_scores
    for g in genes_out:
        gene_scores[g] = 0
    max_score_neighbors = 0
    for g_i in gene_scores:
        score_from_neighbors = 0
        if gssn.has_node(g_i):
            g_neighbors = gssn.neighbors(g_i)
            for g_j in g_neighbors:
                w_ji = round(gssn[g_j][g_i]["weight"] / max_weight, 3)
                score_from_neighbors = score_from_neighbors + (gene_scores[g_j] * w_ji)
            if score_from_neighbors > max_score_neighbors:
                max_score_neighbors = score_from_neighbors
        gene_scores_from_neighbors[g_i] = round(score_from_neighbors, 3)
    for g_i in gene_scores_from_neighbors:
        gene_scores_from_neighbors[g_i] = round(gene_scores_from_neighbors[g_i] / max_score_neighbors, 3)
    return gene_scores_from_neighbors
      

def get_gene_scores(gene_scores_mutations, gene_scores_neighbors):
    gene_scores = {}
    
    for g in gene_scores_neighbors:
        scores = []
        scores.append(gene_scores_mutations[g])
        scores.append(gene_scores_neighbors[g])
        scores.append(gene_scores_mutations[g] + gene_scores_neighbors[g])
        gene_scores[g] = scores
    return gene_scores


def main():
    input_parameters_file = sys.argv[1]
    cp = configparser.ConfigParser() 
    cp.read(input_parameters_file)
    
    cp_input_type = cp["INPUT_TYPE"]
    input_type = int(cp_input_type["TYPE"])
    
    cp_input = cp["INPUT"]
    
    if input_type == 1:
        input_maf_file_name = cp_input["MAF_FILE_NAME"]
    elif input_type == 2:
        input_bmm_file_name = cp_input["BMM_FILE_NAME"]
        input_wmm_file_name = cp_input["WMM_FILE_NAME"]
    input_gene_network_files_name = cp_input["GENE_NETWORK_FILE_NAME"].split(" ")
    
    cp_mutation_weights = cp["VARIANT_CLASSIFICATION_WEIGHTS"]
    mutation_weights = {}
    mutation_weights["Nonsense_Mutation"] = float(cp_mutation_weights["NONSENSE"])
    mutation_weights["Missense_Mutation"] = float(cp_mutation_weights["MISSENSE"])
    mutation_weights["Splice_Site"] = float(cp_mutation_weights["SPLICE"])
    mutation_weights["Frame_Shift_Del"] = float(cp_mutation_weights["FRAMESHIFT_DEL"])
    mutation_weights["Frame_Shift_Ins"] = float(cp_mutation_weights["FRAMESHIFT_INS"])
    mutation_weights["In_Frame_Del"] = float(cp_mutation_weights["INFRAME_DEL"])
    mutation_weights["In_Frame_Ins"] = float(cp_mutation_weights["INFRAME_INS"])
    mutation_weights["3'UTR"] = float(cp_mutation_weights["_3UTR"])
    mutation_weights["5'UTR"] = float(cp_mutation_weights["_5UTR"])
    mutation_weights["Nonstop_Mutation"] = float(cp_mutation_weights["NONSTOP_MUTATION"])
    mutation_weights["Translation_Start_Site"] = float(cp_mutation_weights["TRANSLATION_START_SITE"])
    
    cp_output = cp["OUTPUT"]
    output_folder = cp_output["OUTPUT_FOLDER"]
    output_file = cp_output["OUTPUT_FILE"]
    output_wmm_file_name = output_file + ".wmm"
    output_bmm_file_name = output_file + ".bmm"
    output_am_file_name = output_file + ".am"
    output_gene_mutations_score = output_file + ".mutation.score"
    output_rgn_file_name = output_file + ".sgn"
    output_score_file_name = output_file + ".score"
    output_sample_report_file_name = output_file + ".sample.report"
    output_gene_report_file_name = output_file +  ".gene.report"
    output_network_report_files_name = []
    for gene_network_file_name in input_gene_network_files_name:
        output_network_report_files_name.append(output_file + '_' + gene_network_file_name + ".network.report")
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', filename=output_folder + '/' + 'log.txt')
    logging.info('******** Starting ********')
    
    if input_type == 1:
        maf = pd.read_csv(input_maf_file_name, sep="\t", comment='#', usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification"])
        wmm = get_weighted_mutation_matrix(maf, mutation_weights)
        bmm = get_binary_mutation_matrix(maf)
        util.create_txt_file_mutation_matrix(wmm, output_folder + "/" + output_wmm_file_name)
        util.create_txt_file_mutation_matrix(bmm, output_folder + "/" + output_bmm_file_name)
        am = get_alteration_matrix(bmm)
        util.get_sample_report(am, output_folder + "/" + output_sample_report_file_name)
        util.get_gene_report(bmm, output_folder + "/" + output_gene_report_file_name)
        mutated_genes = get_genes_from_bmm(bmm)
        logging.info('Input Type 1 - END matrix generations')
        
        gene_scores_mutations = get_score_genes(bmm, wmm)
        util.generate_gene_score_mutation_file(gene_scores_mutations, output_folder + "/" + output_gene_mutations_score)
        logging.info('END - getting gene_scores_mutations')
    elif input_type == 2:
        wmm = util.get_mutation_matrix_from_file(input_wmm_file_name)
        bmm = util.get_mutation_matrix_from_file(input_bmm_file_name)
        am = get_alteration_matrix(bmm)
        util.get_sample_report(am, output_folder + "/" + output_sample_report_file_name)
        util.get_gene_report(bmm, output_folder + "/" + output_gene_report_file_name)
        mutated_genes = get_genes_from_bmm(bmm)
        logging.info('Input Type 2 - END matrix generations')
        
        gene_scores_mutations = get_score_genes(bmm, wmm)
        util.generate_gene_score_mutation_file(gene_scores_mutations, output_folder + "/" + output_gene_mutations_score)
        logging.info('END - getting gene_scores_mutations')
    elif input_type == 3:
        gene_scores_mutations = util.get_mutation_score_from_file(output_folder + "/" + output_gene_mutations_score)
        mutated_genes = list(gene_scores_mutations)
    logging.info('END - matrix generations')
    
    gene_network_nx_list = []
    all_genes = set()
    
    for input_gene_network_file_name in input_gene_network_files_name:
        gene_network_file = open(input_gene_network_file_name, 'rb')
        gene_network_nx = read_gene_network_nx(gene_network_file)
        gene_network_nx_list.append(gene_network_nx)
        gene_network_file.close()
        
    consensus_gene_network = combine_gene_network(gene_network_nx_list)
    util.get_network_report(consensus_gene_network, output_folder + "/" + output_file + ".consensus.network.report")
    
    gssn, max_edge_weight = create_gene_strength_spreading_network(consensus_gene_network)
    consensus_gene_score_neighbors = get_score_genes_neighbors(gene_scores_mutations, gssn, max_edge_weight, consensus_gene_network)

    logging.info("END - " + input_gene_network_file_name +  "processing")


    logging.info("END - creating Gene Strength Spreading Networks (GSSNs)")
    logging.info("END - getting scores from neighbors")

    gene_scores = get_gene_scores(gene_scores_mutations, consensus_gene_score_neighbors)
    
    sorted_gene_scores = sorted(gene_scores.items(), key=lambda kv: kv[1][2], reverse=True)
    
    logging.info('END - sorting results')
       
    score_file = open(output_folder + "/" + output_score_file_name + ".mutatedGenes", "w")
    score_file.write("gene\tmutation_score\tscore_from_neighbors\tfinal_score\n")
    for g, scores in sorted_gene_scores:
        if g in mutated_genes:
            score_file.write("{}\t{:.3f}\t{:.3f}\t{:.3f}\n".format(g, scores[0], scores[1], scores[2]))
    score_file.close()

main()
