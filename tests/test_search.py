import unittest

from engage import search
from engage.search import TF, Cluster

import pandas as pd

##################################################

class TestTFInitialization(unittest.TestCase):

    def test_TF__init__(self): 
        tf = TF('ZICL','CAGCTGTG')
        self.assertEqual(tf.name,'ZICL')
        self.assertEqual(tf.motifs,['CAGCTGTG'])

        tf = TF('ZICL',['CAGCTGTG'])
        self.assertEqual(tf.name,'ZICL')
        self.assertEqual(tf.motifs,['CAGCTGTG'])
        
    def test_TF__init__errors(self):
        # Verify TypeErrors with improper initialization 
        with self.assertRaises(ValueError):
            tf = TF('','CAGCTGTG')
        with self.assertRaises(TypeError):
            tf = TF(123,'CAGCTGTG')
        with self.assertRaises(TypeError):
            tf = TF('ZICL',123)
        with self.assertRaises(TypeError):
            tf = TF('ZICL',['CAGCTGTG',123])
        with self.assertRaises(ValueError):
            tf = TF('ZICL',[])


class TestClusterInitialization(unittest.TestCase):

    def test_Cluster__init__(self):
        cluster = Cluster('Test',window=500)
        self.assertEqual(cluster.name,'Test')
        self.assertEqual(cluster.window,500)

        tf = TF('ZICL','CAGCTGTG')
        cluster.add_member(tf,2)
        self.assertEqual(len(cluster.members),1)
        self.assertEqual(type(cluster.members),dict)

    def test_Cluster__init__errors(self):
        cluster = Cluster('Test',window=500)
        with self.assertRaises(TypeError):
            cluster.add_member('ZICL',2)

        tf = TF('ZICL','CAGCTGTG')
        cluster.add_member(tf,2)
        # Verify that we can't add a transcription factor twice
        with self.assertRaises(ValueError):
            cluster.add_member(tf,2)


class TestHelperMethods(unittest.TestCase):

    def test__read_genome(self):
        filename = 'datasets/ciona_2008_KhC1-short.fasta'
        self.assertEqual(type(search._read_genome(filename)), dict)

    def test__motif_to_regex(self):
        # Validate IUPAC replacement and concatenantion between multiple motifs
        motifs = ['CANS','CACAGCTG']
        self.assertEqual(search._motif_to_regex(motifs), 
            'CA(A|C|T|G)(G|C)|CACAGCTG')

        motifs = 'ZED'
        with self.assertRaises(ValueError):
            search._motif_to_regex(motifs)

    def test__find_motif_locations(self):
        # Define the test genome and cluster information to use 
        genome = search._read_genome('datasets/ciona_2008_KhC1-short.fasta')
        cluster = Cluster('Test_Parameters')
        tf = TF('Test_Member','TNCNT')
        cluster.add_member(tf,2)
        
        # Validate proper initialization of the storage objects
        a,b = search._find_motif_locations(genome,cluster)
        self.assertEqual(len(a), 1) 
        self.assertEqual(len(a['KhC1']), 1500)
        self.assertEqual(b, {'Test_Member':2})

        # Check if our method found the match we found through
        # https://regexr.com for our pattern of interest
        self.assertEqual(list(a['KhC1'][132])[0],'Test_Member') 

    def test__consolidate_overlapping_regions(self): 
        test_data = [ 1, 4, 5, 6, 10, 15, 16, 17, 18, 22, 25, 26, 27, 28 ]
        self.assertEqual(len(search._consolidate_overlapping_regions(test_data)), 6) 


class TestMainMethods(unittest.TestCase):

    def test_find_motif_cluster(self):
        filename = 'datasets/ciona_2008_KhC1-short.fasta'
        cluster = Cluster('Test_Parameters')
        tf = TF('Test_Member','TNCNT')
        cluster.add_member(tf,2)

        result_df = search.find_motif_cluster(filename,cluster)
        self.assertEqual(result_df.shape[0], 4)
        