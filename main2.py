#!/usr/bin/python3
# -*- coding: utf-8 -*-
#
# Version: 1.1
#
# The Script retruns all possible cuts and how likely they are for a specific protease and sequence.
#
# This script was developed and tested on Ubuntu 18.04. It requires the installation of pyhton3-biopython and python3-mysqldb
# The meropsrefs.sql file needs to be imported into MYSQL into the database "merops121" and the mysql username and password
# need to be adjusted accordingly. The test data works without the mysql database as the necessary data is cached in the file
# data.py

from Bio import SeqIO
import operator
import MySQLdb
import data

WEIGHTS = [4,8,16,35,20,10,5,2]

CODES = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
         "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
         "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U",
         "Pyl": "O", "Xaa": "X", "Asx": "B", "Glx": "Z", "Xle": "J", "-": "-"}


# WRITES PROTEINASE DATA FOR INTO data.py FILE
# The relative usage of amino acids for the 8 amino acids surrounding
# the cleavage position (-4 to +4) is given in %. The format of data.py
# is a dictionary of dictionaries, where the outer dictionary is a dictionary
# of proteinases and the inner dictionary a dictionary of amino acid residues.
# The amino acid residue is the key and the value is a list of 8 numbers,
# which respresent the  4 amino acids preceeding and following the cleavage
# site.
def _write_data(data): #
    f = open("data.py", "w")
    f.write("data = " + data)
    f.close()

def translate(_3letter): #translates 3-letter-code into 1-letter-code
    if _3letter in CODES.keys():
        return CODES[_3letter]
    else:
        return "X"

def print_results(r): #sorts results according to score, prints them
    sorted_r = sorted(r.items(), key=operator.itemgetter(1))
    print(sorted_r[::-1])
    print('\n')

# This function is not used anymore.
# Assigns each cut in data_list a score according to how likely that cut is.
#
def score(data_list):
    scoredict = {}
    for k in data_list.keys():
        #print(data_list[k])
        scoredict[k] = data_list[k][3] + data_list[k][4]
    return scoredict

# This is the heart of the program. It generates a dictionary of all possible cuts within
# the amino acid sequence.
# It assigns each cut in data_list a score according to how likely that cut is (new system)
def new_score(data_list):
    scoredict = {}
    for k in data_list.keys():
        score = 0
        for i in range(len(data_list[k])):
            score += data_list[k][i]*WEIGHTS[i]
        scoredict[k] = score
    return scoredict

def get_data(protease, data): #checks if data is in data.py, adds data to data.py if it's not, returns data
    if protease in data.data:
        print('FOUND IN DATA {0}'.format(protease))
    else:
         table = calc_table(protease)
         data.data[protease] = table
         _write_data(str(data.data))
         print('ADDED TO DATA {0}'.format(protease))
    return data.data[protease]


# Gets data from the sql database, returns a dictionary of how many times a specific amino acid residue
# occurs at the positons -4 to +4 relative to the cut site.
#
#
def calc_table(protease):
    conn = MySQLdb.connect(host='localhost', user='merops_user', passwd='12345678', db='merops121')
    cursor = conn.cursor()
    # Returns the cleavage context (-4 to +4) for all substrates of a specific proteinase
    command = "SELECT Site_P4, Site_P3, Site_P2, Site_P1, Site_P1prime, Site_P2prime, Site_P3prime, Site_P4prime FROM Substrate_search WHERE Protease = '{0}'".format(protease)
    cursor.execute(command)
    substrate_list = cursor.fetchall()
    print('Received {0} rows from the following database request:\n{1}'.format(len(substrate_list), command))

    if substrate_list: #checks if sql returned something
        p_data = {} #writes data from sql to a dictionary
        for aminoacid in CODES.values():
            p_data[aminoacid] = [0, 0, 0, 0, 0, 0 , 0, 0]
        # Loop through list of substrates
        for substrate in substrate_list:
            # Loop through aminoacid positions -4 to +4 (0 does not exist)
            for i in range(8):
                p_data[translate(substrate[i])][i] += 1
        # Because different proteinases have different numbers of known
        # substrates, we need to divide by the number of substrates
        # in order to make values comparable between proteinases.
        for aminoacid in p_data.keys():
            for i in range(8):
                p_data[aminoacid][i] /= len(substrate_list)
    else:
        raise ValueError('Protease named {0} not found'.format(protease))

    print(p_data)

    return p_data

def get_sequence_from_str(string): #isolates sequence (e.g. "---MHLLG") from sequence string (e.g. "1/---MHLLG/8")
    sequence = ""
    for char in string:
        if char in CODES.values():
            sequence += char
    if "-" in sequence[:4]:
        sequence = sequence[1:]
    return sequence

def get_number_from_str(string): #isolates number (e.g. 1) from sequence string (e.g. "1/---MHLLG/8")
    half = string[:int(len(string)/2)]
    i_list = list(filter(str.isdigit, half))
    number = ""
    for i in i_list:
        number += i
    return int(number)

def lookup(protease, sequence_file, file_format="fasta", **kwargs): #main function
    sequence = str(next(SeqIO.parse(open(sequence_file), file_format)).seq)
    sequence = "---" + sequence + "---"
    protease_data = get_data(protease, data) #gets data
    #print(len(sequence))

    cut_dict = {} #creates a dict of all possible cuts and a list of how many times the protease cuts the sequence after that amino-acid
    for i in range(3, len(sequence)-4):
        seq_string = str(i-6) + "/" + sequence[i-3:i+5] + "/" + str(i+2)
        #print(seq_string)
        p_string = get_sequence_from_str(seq_string)
        #print(p_string)
        l = []
        for j in range(len(p_string)):
            l.append(protease_data[p_string[j]][j])
        cut_dict[seq_string] = l
    #print(cut_dict)

    scoredict = new_score(cut_dict) #scores each cut according to likelihood
    #print(scoredict)

    if "nocut" in kwargs:
        filtered = {}
        for (key, value) in scoredict.items():
            if int(get_number_from_str(key)) >= kwargs["nocut"]-2:
                filtered[key] = value
        scoredict = filtered

    print_results(scoredict) #prints the results
    return scoredict

if __name__ == "__main__":
    # Demo with cached proteinase
    data_list = lookup('kallikrein-related peptidase 3', "P49767.fasta", nocut=21)
    # Demo with mysql lookup
    data_list = lookup('matrix metallopeptidase-9', "P49767.fasta", nocut=21)
    # Demo with mysql lookup
    data_list = lookup('thrombin', "P49767.fasta", nocut=21)
