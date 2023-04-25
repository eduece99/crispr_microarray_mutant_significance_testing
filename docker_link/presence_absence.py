#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:16:21 2023

@author: Edmund Duesbury

==INFO==

The premise of this script is to create mutation-gene knock-out pairs and find those
which for a given mutation, has significantly different
cell counts vs the absence of the mutation.

The data is based on cell models, and two source data files.  The first has 
the mutants listed per cell model type, and the other has model type-gene KOs.  These
after being melted to long type, can be merged to yield the fold counts for each
mutant-KO and thus tested on.

The Brunner-Munzel test was used as opposed to the Mann-Whitney-U test, as it does
not assume equivariance in the two samples, and is otherwise standardised for sample
size (whereas the U test is not)

"""

#import pyspark
import os
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests



# read data
work_dir = "/home/edmund/Documents/programming/tests/Mosaic_TX/docker_link/"
input_data_dir =  "./input_data/"

def read_data_file( workdir : str, filename : str ) -> pd.DataFrame :
    
    """
    Convenience function for reading in file in defined work directory

    Parameters
    ----------
    workdir : str
        Directory where file is located
    filename : str
        File Name with extension

    Returns
    -------
    Pandas DataFrame 

    """
    return( pd.read_table( workdir + filename ) )


def pairwise_mutations_kos( mt : pd.DataFrame, ko : pd.DataFrame ) -> pd.DataFrame :
    """
    

    Parameters
    ----------
    mt : pd.DataFrame
        Pandas DataFrame containing presence/absence of mutations (rows) per cell model (columns)
    ko : pd.DataFrame
        Pandas DataFrame depicting cell count fold changes for each gene knockout (columns) per cell model (rows)

    Returns
    -------
    Long format Pandas DataFrame showing the cell count fold change for each mutation-knockout pair.

    """
    
    # convert tables to long format for joining and later grouping for pairwise tests
    m_long = mt.melt( id_vars=["Mutation"], var_name="Model", value_name="Mutation_Presence")
    ko_long = ko.melt( id_vars=["Model"], var_name="gene_KO", value_name="qn_cellcount")
    
    #m_ko_wide = m_long.merge( ko_pdf, on="Model" ).drop("Model", axis=1)
    m_ko = m_long.merge( ko_long, on="Model" ).drop("Model", axis=1)
    
    return( m_ko )
    


def test_presence_vs_absence( df : pd.DataFrame, plot_hist : bool = False  ) -> pd.Series :
    """
    

    Parameters
    ----------
    df : pd.DataFrame
        A dataframe which contains the data only for a specific mutation-knockout pair.
        Must have "Mutation_Presence" (0 or 1) and qn_cellcount (double) columns
    plot_hist : bool, optional
        Set to True, to plot the presence vs absence distributions. The default is False.

    Returns
    -------
    Summary statistics (including Brunner-Munzel test p-value) comparing presence of mutation vs absence.

    """
    
    #print(df)
    sample_p0 = df.loc[ (df["Mutation_Presence"] == 0), "qn_cellcount" ] 
    sample_p1 = df.loc[ (df["Mutation_Presence"] == 1), "qn_cellcount" ]
    
    if plot_hist:
        sample_p0.hist()
        sample_p1.hist()
        
    mwu_stats = stats.brunnermunzel(sample_p0, sample_p1)
    #mwu_stats.absence_n = 0
    
    pval = mwu_stats.pvalue 
    
    # replace NaN with very small values
    if np.isnan( pval ) :
        pval = 1e-09
    
    output = {
            "absence_n" : len(sample_p0), 
            "presence_n" : len(sample_p1), 
            "absence_mean" : sample_p0.mean(), 
            "presence_mean" : sample_p1.mean(), 
            "pvalue" : pval
    }
    output_series = pd.Series( output )
    
    return( output_series )



def fdr_pvalues( df : pd.DataFrame ) -> pd.DataFrame :
    """
    

    Parameters
    ----------
    df : pd.DataFrame
        input dataframe with the "pvalue" column representing p values of a statistical test.

    Returns
    -------
    Pandas DataFrame as input but with the "fdr_corrected_p" column, representing 
    p values corrected using the Benjamini-Hochberg false discovery rate calculation.

    """
    
    fdr_pvals = multipletests( df["pvalue"], alpha=0.05, method="fdr_bh" )[1]
    df["fdr_corrected_p"] = fdr_pvals
    
    test_results_asc = df.sort_values(by="fdr_corrected_p")
    
    return( test_results_asc )
    


def run_example( plot_hist_row : int = -1 ) -> str :
    
    """

    Parameters
    ----------
    plot_hist_row : int
        plot histograms of an example (row index)

    Returns
    -------
    str 
        Text and table depicting mutation-gene KO combinations whose cell fold counts
        are statistically significant when comparing presence vs absence of mutation.

    """
    
    #os.chdir( "." )  # switch to relative paths from hereon
    
    mutations_pdf = read_data_file( input_data_dir, "Mutations.tsv" )
    ko_pdf = read_data_file( input_data_dir, "Gene_KOs.tsv" )  
    
    
    m_ko = pairwise_mutations_kos( mutations_pdf, ko_pdf )
    
    # how is the data distributed? 
    #ko_pdf.set_index("Model").transpose().hist()
    
    # the histograms do not look normally distributed, and it is not really clear 
    # what distribution to use, so non-parametric tests are likely better than parametric
    
    
    # non paired nonparametric test - brunnermunzel test
    # would have used Mann-Whitney-U but it assumes equivalent variance between samples
    # Brunner-Munzel does not have this assumption
    # This is despite the data being quantile-normalised, which is supposed to make statistical properties the same
    # because the mutant/non-mutant distributions have not been normalised, as they are subsets here
    #sample_p0 = m_ko.loc[ (m_ko["Mutation"] == "Gene1_mut") & (m_ko["Mutation_Presence"] == 0), "GeneA_KO"]
    #sample_p1 = m_ko.loc[ (m_ko["Mutation"] == "Gene1_mut") & (m_ko["Mutation_Presence"] == 1), "GeneA_KO"]
    
    #stats.mannwhitneyu(sample_p0, sample_p1)
    
    
    test_results = m_ko.groupby( ["Mutation", "gene_KO"] ).apply( test_presence_vs_absence )
    
    # p value correction 
    test_results_asc = fdr_pvalues( test_results )
    
    # choose significant matches and then print presence vs absence means & sample sizes
    
    
    results = test_results[ test_results["fdr_corrected_p"] < 0.05 ]
    
    # plot some significant examples to see why they differ
    
    if plot_hist_row >= 0:
            
        test_presence_vs_absence(  m_ko.loc[ 
            (m_ko["Mutation"] == test_results_asc.index[plot_hist_row][0]) & 
            (m_ko["gene_KO"] == test_results_asc.index[plot_hist_row][1]), 
        ], plot_hist = True )

    # one sided tests

    
    return( results )


if __name__ == "__main__":
    
    # for testing only
    #s.chdir( work_dir )
    
    print("Mutation-Gene KO combinations with significant p values \n")
    sig_pairs = run_example( 0 )
    print( sig_pairs )
    
    
