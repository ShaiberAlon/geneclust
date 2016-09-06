#!/usr/bin/env python
import matplotlib.pyplot as plt

import numpy as np
import sklearn.mixture as mixture



def run_GMM(gene_coverage):
	 d1 = mixture.GMM(n_components=5, covariance_type='full', random_state=None, thresh=None, tol=0.001, verbose=1, min_covar=0.001, n_iter=10, params='wmc', init_params='wmc')
	 gene_coverage=np.reshape(gene_coverage,(-1,1))
	 C=d1.fit_predict(gene_coverage)
	 return C,d1

def create_within_samples_clusters():
    file_name="/Users/alonshaiber/github/genelearn/sandbox/Bfrag_positive_gene_coverage_normalized_to_100.txt"
    gene_coverage=np.loadtxt(file_name, delimiter="\t", skiprows=1)
    Ngenes , Nsamples = gene_coverage.shape
    Nsamples=Nsamples-1 # first column is the gene number
    GMM_results=np.empty([Ngenes,Nsamples]) # generate an empty matrix for GMM results
    gene_ids=gene_coverage[:,0]
    #file_name="/Users/alonshaiber/github/genelearn/sandbox/one_sample_example.txt" # TODO: will change to a file with all metagenomic samples and run a for loop on the columns
    #gene_ids=np.loadtxt(file_name, delimiter="\t", skiprows=1, usecols=[0]) # gene id array
    #gene_coverage=np.loadtxt(file_name, delimiter="\t", skiprows=1, usecols=[1]) # TODO: will change usecol=[1] to usecol=[i] in a for loop going over all metagenomic samples
    for idx in range(Nsamples): 
    	idx=idx+1 # the coverage data starts at the second column
	# all zero coverage genes form one cluster
	zero_coverage_genes_indexes=np.where(gene_coverage[:,idx] == 0)
    	zero_coverage_gene_ids=gene_ids[zero_coverage_genes_indexes[0]]
    	# Running the GMM only on the genes that have non-zero coverage
    	gene_coverage_non_zero_indexes=gene_coverage[:,idx].nonzero()
    	gene_coverage_non_zero=gene_coverage[gene_coverage_non_zero_indexes,idx]
    	Clusters,GMM1=run_GMM(gene_coverage_non_zero)
	GMM_results[zero_coverage_genes_indexes,idx-1]=0
	GMM_results[gene_coverage_non_zero_indexes,idx-1]=Clusters+1 #adding 1 since the cluster number 0 is for genes with zero coverage
    output_file='GMM_out.txt' #FIXME: change so the name of output is an input argument
    np.savetxt(output_file,GMM_results,delimiter="\t")
    return GMM_results
# End of create_within_samples_clusters()
def metric1(a1,a2):
    # takes two gene_cluster_assignment arrays (numpy arrays) and computes the distance between the genes
    # The input arrays should contain the cluster assignment per sample
    # The distance is calculated as the percent of the samples in which the two genes appear in the same cluster
    N=len(a1)
    d=sum(a1==a2)*1.0/N
    return d

def cluster_genes(GMM_results):
    metric = sklearn.neighbors.DistanceMetric.get_metric('pyfunc', func=metric1)	
if __name__ == '__main__':
    #parser = argparse.ArgumentParser(description="Run GMM on gene coverage data")
    #parser.add_argument("", help="")
    #parser.add_argument("-s", "--length_sort", help="sort sequences based on length (excluding trailer)", 
    #    action="store_true")
#
    args = parser.parse_args()
    GMM_results=create_within_samples_clusters
    gene_clustering_results=cluster_genes(GMM_results)
    ## 
    #plt.hist(Clusters)
    #plt.title("Clusters Histogram")
    #plt.xlabel("Value")
    #plt.ylabel("Frequency")

    #fig = plt.gcf()
    #plt.show()
