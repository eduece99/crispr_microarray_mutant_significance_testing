#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  1 18:35:36 2023

@author: Edmund Duesbury
"""

"""
UNIT test ideas

"""


import presence_absence
import unittest
from pandas.api.types import is_numeric_dtype

data_dir = "./input_data/"

# are the files there, and readable?
class TestMutationsFileContent(unittest.TestCase):
    
    def setUp(self):
        self.df = presence_absence.read_data_file( data_dir, "Mutations.tsv")
        
    def test_mutations_file_not_empty(self):
        
        # size is not 0
        actual = len(self.df)
        expected = 0
        self.assertNotEqual(actual, expected)
    
    # Are there missing values?
    def test_mutations_file_no_nan(self):
        # check missing values
        self.assertTrue( all( self.df.isna().sum() == 0 ) )
    
    # not strictly boolean but no idea what else to call this test
    def test_mutations_file_booltype(self):
        # values should be 0 or 1
        f_df = self.df.loc[ : , ~self.df.columns.isin(['Mutation']) ]
        
        self.assertTrue( all( f_df.isin([0,1]) ) )
        
    
class TestKOFileContent(unittest.TestCase):
    
    def setUp(self):
        self.df = presence_absence.read_data_file( data_dir, "Gene_KOs.tsv")
        
    def test_ko_file_not_empty(self):
        
        # size is not 0
        actual = len(self.df)
        expected = 0
        self.assertNotEqual(actual, expected)
    
    def test_ko_file_no_nan(self):
        # check missing values
        self.assertTrue( all( self.df.isna().sum() == 0 ) )
    
    # all data (apart from first column) should be numeric
    def test_ko_file_numerictype(self):
        # values should be numeric
        f_df = self.df.loc[ : , ~self.df.columns.isin(['Model']) ]
        
        self.assertTrue( all(f_df.apply(is_numeric_dtype)) )
        
    

# On merge are there any data?  Are the columns and data types correct?
class TestMergedFileContent(unittest.TestCase):
    
    def setUp(self):
        self.mutations_pdf = presence_absence.read_data_file( data_dir, "Mutations.tsv")
        self.ko_pdf = presence_absence.read_data_file( data_dir, "Gene_KOs.tsv")
        self.m_ko = presence_absence.pairwise_mutations_kos( self.mutations_pdf, self.ko_pdf )
        
    def test_merge_file_not_empty(self):
        
        # size is not 0
        actual = len(self.m_ko)
        expected = 0
        self.assertNotEqual(actual, expected)
        

class TestExampleContent(unittest.TestCase):
    
    def setUp(self):
        self.results = presence_absence.run_example()
        
    def test_results_pairreference(self):
        # Should be 5 specific examples appearing with a pvalue < 0.05 in specific order
        
        gene_ref_list = ['Gene1_mut', 'Gene1_mut', 'Gene1_mut', 'Gene2_mut', 'Gene4_mut', 'Gene5_mut', 'Gene8_mut']
        mutation_ref_list = ['GeneB_KO', 'GeneH_KO', 'GeneI_KO', 'GeneB_KO', 'GeneB_KO', 'GeneD_KO', 'GeneB_KO']
        ref_pairs = list( zip(gene_ref_list, mutation_ref_list) )
        
        result_pairs = list( self.results.index.values )
        
        print(ref_pairs)
        print(result_pairs)
        
        # expected mutation-knockout pairs appear
        self.assertEqual(ref_pairs, result_pairs)