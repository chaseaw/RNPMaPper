#!/usr/bin/env python

'''
--------------------------------------------------------------------- 
This file is a part of RNP-MaPper, and is licensed under the terms   
of the MIT license. Copyright 2022 Chase Weidmann.                   
---------------------------------------------------------------------

Created April 2021
@author: chasew
Parses a shapemapper2 output profile for RNP-MaP analysis

Usage: python SM_parser.py [options] --SMtxt
'''
import os,sys,argparse
import numpy as np
import pandas as pd

## shortcut output for dataframe to useable .csv file ##
def write_RNP_csv(df,outname):
	
	if not outname.endswith(".csv"):
		outfile = outname + ".csv"
	else:
		outfile = outname
	
	df = df.replace([np.inf, -np.inf], np.nan)

	df.to_csv(outname+".csv",header=True,index=True,na_rep='nan')

## shortcut output to create nucleotide list .csv files (for warnings) ##
def write_index_list(df,colname,outname):
	if not df.empty:
		indices = list(df.index.values)
		values = list(df[colname])
		assert len(indices) == len(values),"Column and index lengths do not match."
	
		with open(outname + ".csv",'w') as outfile:
			for i in range(len(indices)):
				outfile.write("{},{}\n".format(indices[i],values[i]))

## Converts shapemapper information into dataframe compatible with downstream RNP-MaP analysis steps ##
def importshapemapper(profile,dir_prefix=""):

	# columns needed for RNP-MaP
	info_cols = ['Sequence','Modified_mutations','Modified_effective_depth','Untreated_mutations','Untreated_effective_depth'] 

	# import data
	df = pd.read_csv(profile,index_col=0,sep='\t')[info_cols]

	# rename columns 
	df = df.rename(columns={"Modified_mutations": "XL_mutations", "Modified_effective_depth": "XL_depth","Untreated_mutations": "NoXL_mutations", "Untreated_effective_depth": "NoXL_depth"}) # rename columns

	# Changes NoXL values of 0 to 1 where there are mutations in the XL sample
	# (allows taking a conservative reactivity measurement that would otherwise fail at divide by 0)
	replace_mask = ((df['NoXL_mutations'] == 0) & (df['XL_mutations'] > 0))
	write_index_list(df[replace_mask],"Sequence",os.path.join(dir_prefix,"Rounded_NoXL_Rate_NTs"))
	df.loc[replace_mask, 'NoXL_mutations'] = 1 

	# add rate columns
	df['XL_rate'] = df['XL_mutations']/df['XL_depth']
	df['NoXL_rate'] = df['NoXL_mutations']/df['NoXL_depth']

	# get the log of rates and add columns
	with np.errstate(divide='ignore', invalid='ignore'): # ignores divide by zero warnings
		df['XL_log10rate'] = np.log10(df['XL_rate'])
		df['NoXL_log10rate'] = np.log10(df['NoXL_rate'])


	return df

## filters out nucleotides that don't pass background and depth checks ##
def bgfilter(df,LoDepth = 5000,HiBG = 0.05,dir_prefix=""):
	write_index_list(df[(df['NoXL_rate'] >= (1-HiBG*2))],"Sequence",os.path.join(dir_prefix,"Check_FASTA")) # outputs a csv list of nucelotides that are likely incorrect in the input FASTA that will affect analysis

	# removes nucleotides not passing the depth filter
	depth_mask = ((df['NoXL_depth'] >= LoDepth) & (df['XL_depth'] >= LoDepth))
	depth_pass = df[depth_mask] 

	# removes nucleotides with specified high background rates
	bg_mask = (depth_pass['NoXL_rate'] < HiBG)
	bg_pass = depth_pass[bg_mask]

	# removes infinite and NaN-containing rows
	nan_pass = bg_pass[~bg_pass.isin([np.nan, np.inf, -np.inf]).any(1)] 

	return nan_pass

## splits a dataframe into dictionaries subset by specified nucleotides ##
def split_nts(df,NTs,dontsplit=False):
	NT_dfs = {}
	
	if dontsplit:
		NT_dfs['N'] = df.copy()

	else:
		for letter in NTs:
			cases = [letter.upper(),letter.lower()]
			NT_dfs[letter] = (df[df['Sequence'].isin(cases)])

	return NT_dfs

if __name__ == '__main__': # allows another python script to import the functions

	parser = argparse.ArgumentParser(description="Parses a shapemapper2 output profile for RNP-MaP analysis", usage='%(prog)s [options] --SMtxt ',add_help=False) # a description of the function

	required = parser.add_argument_group('Required Input', 'These specifications are necessary to run.')

	required.add_argument("--SMtxt",help="Input shapemapper2 profile.txt", required=True)

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the program runs.')
	data_opt.add_argument('-o',"--outname",default=None,help="Specify the name of the output files.")
	data_opt.add_argument('-d',"--LoDepth",default=5000,help="Specify depth threshold for inclusion of a nucleotide. Default = 5000.")
	data_opt.add_argument('-b',"--HiBG",default=0.05,help="Specify the background rate to exclude a nucleotide. Default = 0.05.")

	help_opt = parser.add_argument_group('Help')
	help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")

	args = parser.parse_args()

	parsedInfo = importshapemapper(args.SMtxt)
	BG_pass = bgfilter(parsedInfo,LoDepth=args.LoDepth,HiBG=args.HiBG)

	if args.outname is not None:
		outfile = args.outname
	else:
		outfile = os.path.splitext(args.SMtxt)[0]

	write_RNP_csv(parsedInfo,outfile)
	write_RNP_csv(BG_pass,"{}_bgpass".format(outfile))
	
	nts = ["A","C","G","U"]
	NTs = split_nts(BG_pass,nts)
	for nt in nts:
		write_RNP_csv(NTs[nt],"{}_bgpass_{}s".format(outfile,nt))	

