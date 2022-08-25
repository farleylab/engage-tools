from Bio import SeqIO

import pandas as pd 
import re 

from itertools import groupby
from operator import itemgetter

from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import natsort

##################################################
# Class Definitions 
##################################################

class TF: 
    name = "" 
    motifs = []

    def __init__(self,name,motifs):
        if type(name) != str: 
            raise TypeError("The transcription factor name should be type(str)")
        elif len(name) == 0:
            raise ValueError("Transcription factor names cannot be empty")

        if type(motifs) == list and all(type(e) == str for e in motifs) == False: 
            raise TypeError(
                "All transcription factor motifs in a list should be type(str)")
        elif type(motifs) == list and len(motifs) == 0:
            raise ValueError("Transcription factor motifs cannot be empty")
        elif type(motifs) != list and type(motifs) != str: 
            raise TypeError(
                "A single transcription factor motif should be type(str)")

        self.name = name
        self.motifs = motifs
        if type(motifs) == str:
            self.motifs = [motifs]

class Cluster: 
    name = ""
    window = 0
    members = {}

    def __init__(self,name,window=100,members=None):
        self.name = name
        self.window = window

        if members == None:
            self.members = {}

    def add_member(self,tf,number):
        """Adds a transcription factor to the cluster object"""
        if type(tf) != TF: 
            raise TypeError("Cluster members should be type(TF)")
        if type(number) != int: 
            raise TypeError("The number of cluster members should be type(int)")

        if tf not in self.members: 
            self.members[tf] = number
        else: 
            raise ValueError("The cluster member already exists in the object")

    def __str__(self): 
        output =  f"[ {self.name} ] Cluster Parameters:\n"
        output += "====================================\n"
        output += f"Window Size = {self.window}\n"
        output += "Transcription Factors =\n"
        if len(self.members) == 0:
            output += "* None Defined\n"
        else:
            for tf in self.members: 
                output += f"* {tf.name}: {self.members[tf]}x Sites\n"
                output += f"- Motifs: {tf.motifs}\n"

        return(output)

##################################################
# Helper Methods
##################################################

def _read_genome(filename):
    """Opens a FASTA genome file then covert contents to a SeqIO dict object"""
    return(SeqIO.to_dict(SeqIO.parse(open(filename,'r'), "fasta")))

def _motif_to_regex(motifs):
    """Converts a list of motifs to a regular expression"""
    motifs_string = "|".join(motifs).upper()
    
    # Conversions between non-A,C,T,G bases to their A,C,T,G counterparts 
    # Referenced from: https://www.bioinformatics.org/sms/iupac.html
    iupac_conversions = { "A": "A", "C": "C",
                          "T": "T", "G": "G",
                          "R": "(A|G)", "Y": "(C|T)",
                          "S": "(G|C)", "W": "(A|T)",
                          "K": "(G|T)", "M": "(A|C)",
                          "B": "(C|G|T)", "D": "(A|G|T)",
                          "H": "(A|C|T)", "V": "(A|C|G)",
                          "N": "(A|C|T|G)"}
    
    # Translate between IUPAC bases and the traditional nucleotide bases 
    updated_motifs_string = ""
    for nucleotide in motifs_string:
        nt = nucleotide
        if nt in iupac_conversions: 
            nt = iupac_conversions[nt]
        elif nt not in iupac_conversions and nt.isalpha() is True:
            raise ValueError("Invalid nucleotide detected")

        updated_motifs_string += nt 

    return(updated_motifs_string)

def _search_chromosome(search_pattern_parameters,chromosome_sequence): 
    """Search through a chromosome for clusters"""
    # Define the chromosome sequence and regions that can contain a site
    search_member_locations = {e:{} for e in range(len(chromosome_sequence))}

    # Search through each chromosome and define the location of sites
    for member in search_pattern_parameters: 
        regions_found = search_pattern_parameters[member].finditer(chromosome_sequence)
        for site in regions_found: 
            if member not in search_member_locations[site.start()]: 
                search_member_locations[site.start()][member] = tuple()
            search_member_locations[site.start()][member] = \
                (site.start(),site.end(),chromosome_sequence[site.start():site.end()])
    
    return(search_member_locations)

def _find_motif_locations(genome_dictionary,cluster_parameters,num_processes=4): 
    """Find the locations of all desired motifs in the genome"""
    # For each transcription factor motif present in the cluster, generate
    # a separate regular expression to search for it within the string later
    cluster_members = cluster_parameters.members
    search_pattern_parameters = {} # regex patterns for each member
    search_number_parameters = {} # number of appearances for each member
    for member in tqdm(cluster_members,
                       desc="Processing Search Parameters"): 
        search_pattern_parameters[member.name] = re.compile(
            _motif_to_regex(member.motifs),re.IGNORECASE)
        search_number_parameters[member.name] = cluster_members[member]

    # Go through each of the chromosomes and search for the clusters 
    # with multiprocessing on default 4 processes
    search_member_locations = {}
    chromosome_list = natsort.natsorted(genome_dictionary.keys())
    chromosome_sequences = [str(genome_dictionary[c].seq) for c in chromosome_list]

    # !TODO : Fix multiprocessing.pool + tqdm tracking 
    # The current method of checking if each chromosome is done processing 
    # is incorrect -- tqdm currently tallies /when/ chromosome is initiated

    with Pool(processes=num_processes) as pool: 
        search_results = pool.imap(
            partial(_search_chromosome,search_pattern_parameters),
            tqdm(chromosome_sequences,desc="Searching Chromosomes for Clusters"))
        for chromosome,result in zip(chromosome_list,search_results): 
            search_member_locations[chromosome] = result

    return(search_member_locations,search_number_parameters)

def _consolidate_overlapping_regions(cluster_indices): 
    """
    Determines groupings of sequential values for a list of integers

    **Referenced from: 
    https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
    Answer from Dominic Rodger (Top Answer) with Python 3 mods from Kevin & Euler_Salter
    """
    output = []
    for k,g in groupby(enumerate(sorted(cluster_indices)), lambda ix : ix[0] - ix[1]): 
        output.append(list(map(itemgetter(1), g)))
    return(output) 

##################################################
# Main Method
##################################################

def find_motif_cluster( filename,cluster_parameters,
                        minimum_flanking=5,exact_number=False,
                        num_processes=4 ): 
    """
    Given a genome file, find all clusters of desired motifs in the genome.

    find_motif_cluster() loads a genome .fasta/.fa file and searches chromosome
    by chromosome for clusters of user-defined transcription factor motifs

    Input: 
        filename - Location of the genome .fasta/.fa file 
        cluster_parameters - User-defined Cluster[] object containing the 
            motifs that are present within the cluster and the number of 
            times each motif appears within a given window size
    
    Returns: 
        Pandas DataFrame object containing the locations of the clusters
    """
    genome_dictionary = _read_genome(filename)
    if type(cluster_parameters) != Cluster: 
        raise TypeError("Cluster parameter provided should be type(Cluster)")

    # Find the location of all desired motifs in the genome 
    # and the number of times they should appear within a given window
    search_member_locations,search_number_parameters = \
        _find_motif_locations(genome_dictionary,cluster_parameters,num_processes)

    window_size = cluster_parameters.window

    # Go through each chromosome and evaluate a sliding window to see if we 
    # meet the minimum requirements for the number of motif appearances

    # !TODO : Incorporate multiprocessing when collapsing overlapping regions

    data = []
    for chromosome in tqdm(search_member_locations,
                           desc="Consolidating Overlapping Regions",
                           total=len(search_member_locations)):
        clusters_found = set()
        chromosome_search_results = search_member_locations[chromosome]
        for chromosome_index_start in range(len(chromosome_search_results)-window_size+1): 
            motif_tallies = { member:0 for member in search_number_parameters }
            chromosome_index_end = chromosome_index_start + window_size 

            # Evaluate the window for the motifs that are present and tally them
            for chromosome_index in range(chromosome_index_start,chromosome_index_end): 
                for member in chromosome_search_results[chromosome_index]: 
                    # Ensure that the beginning of the motif has some 
                    # breathing room (set by minimum_flanking) at the start of 
                    # the window, and that the end of the motif isn't outside
                    if chromosome_search_results[chromosome_index][member][0] < \
                        chromosome_index_start + minimum_flanking:
                        continue
                    if chromosome_search_results[chromosome_index][member][1] > \
                        chromosome_index_end - minimum_flanking:
                        continue 
                    motif_tallies[member] += 1

            # Save the region if we meet the motif appearance requirements for
            # either exact matches or "minimum match" requirements
            if exact_number is True:
                if motif_tallies == search_number_parameters: 
                    clusters_found = clusters_found | \
                        set(range(chromosome_index_start,chromosome_index_end))
            else: 
                total_number_parameters_passed = 0
                for member in search_number_parameters: 
                    if motif_tallies[member] >= search_number_parameters[member]: 
                        total_number_parameters_passed += 1
                if total_number_parameters_passed == len(search_number_parameters): 
                    clusters_found = clusters_found | \
                        set(range(chromosome_index_start,chromosome_index_end))

        for cluster_region in _consolidate_overlapping_regions(clusters_found): 
            if len(cluster_region) < window_size:
                continue
            else:
                data.append({ "Chromosome": chromosome, 
                              "Start": cluster_region[0], 
                              "End": cluster_region[-1] })

    return(pd.DataFrame(data))
