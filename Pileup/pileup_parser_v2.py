import os 
import re
from selectors import BaseSelector
import pandas as pd
import numpy as np

def open_and_cut(directory):
#this function opens pileup files, trims their headers, and deposits a trimmed TSV file
#this file can then be opened and then trimmed with pandas
#it is important to create the deposition folders beforehand and manually sort the files
    list_of_files = os.listdir(directory)
    print(list_of_files)
    for pileup in list_of_files:
        print(pileup)
        f = open((directory+pileup), "r")
        data = f.readlines()
        length = len(data)
        start = 0
        for line_no, line in enumerate(data):
            if "#CHROM" in line:
                start = line_no
                break
        trimmed_data = data[start:length]
        with open("trimmed_pileups/trm_"+pileup+".tsv", "w") as trimmed:
            for line in trimmed_data:
                trimmed.write(line)
            
def isolate_bases(directory):
#this function isolates reads with specific base pairs and exports a TSV file
#tsv files were used due to some output data including commas
    list_of_files = os.listdir(directory)
    for pileup in list_of_files:
        infile = directory+pileup
        print(infile)
        input_data = pd.read_csv(infile, sep='\t', header=0)
        clipped_input = input_data
        clipped_input['INFO'] = clipped_input['INFO'].str.split(pat=';')
        print(clipped_input["INFO"])
        output = pd.DataFrame()
        output['CHROM'] = clipped_input["#CHROM"]
        output['POS'] = clipped_input["POS"]
        output['REF'] = clipped_input["REF"]
        output['ALT'] = clipped_input["ALT"]
        output["BASES"] = clipped_input['REF'] + "," + clipped_input['ALT']
        print(output["BASES"])
        #these specific positions may change based on the pileup output that you generate
        #make sure to change these based on your specific mcftools mpileup output
        output["DEPTHS"] = clipped_input["INFO"].str[1]
        print(output["DEPTHS"])
        output["DEPTHS"] = output["DEPTHS"].str[3:]
        print(output["DEPTHS"])
        output["MOST_ABUNDANT_EDIT"] = output["BASES"].str[2]
        output = output.reset_index()
        math = pd.DataFrame()
        math["DEPTH"] = pd.DataFrame(output["DEPTHS"].str.split(','))
        math2 = pd.DataFrame()
        math2 = pd.DataFrame(math["DEPTH"].tolist())
        print(math2)
        math2 = math2.fillna(0)
        math2 = math2.astype(int)
        sums = pd.DataFrame()
        sums["TOTAL_DEPTH"] = math2.sum(axis=1)
        math2_freqs = pd.DataFrame()
        math2_freqs = math2.iloc[:,[0,1,2,3]].div(sums["TOTAL_DEPTH"].values,axis=0)
        math2_freqs.columns = ["FREQ1","FREQ2","FREQ3","FREQ4"]
        o1 = pd.DataFrame()
        o1 = output.join(sums)
        o2 = o1.join(math2_freqs)
        o2.to_csv("subsetted_pileups/clipped_"+pileup,sep='\t')

open_and_cut()
isolate_bases()

#outputs should be manually collated based on condition before analysis with the .Rmd script
