#!/usr/bin/env python

'''
--------------------------------------------------------------------- 
This file is a part of RNP-MaPper, and is licensed under the terms   
of the MIT license. Copyright 2022 Chase Weidmann.                   
---------------------------------------------------------------------

Created April 2021
@author: chasew
Calls protein binding sites based on shapemapper2 output of an RNP-MaP experiment

Usage: python RNP_MaP.py [options] --SMtxt <shapemapper2 profile txt>
'''
import os
import sys
import argparse
import pandas as pd
import numpy as np
import SM_parser as SMparse

## Sets global prefix ##
def cure_prefix(infile, prefix=None):
	return prefix if prefix else os.path.splitext(infile)[0]

## Establishes specified main and temporary directories ##
def Make_Dirs(prefix):
	main_folder = "{}_RESULTS".format(prefix)
	os.makedirs(main_folder, exist_ok=True)
	return main_folder

## Call to shapemapper_parser.py ##
def parse_SM(SMprofile, LoDepth=5000, HiBG=0.05, nts=["A","C","G","U"], dir_prefix=""):
	parsedInfo = SMparse.importshapemapper(SMprofile, dir_prefix=dir_prefix)
	BG_pass = SMparse.bgfilter(parsedInfo, LoDepth=LoDepth, HiBG=HiBG, dir_prefix=dir_prefix)
	NTs = SMparse.split_nts(BG_pass, nts)
	return [parsedInfo, BG_pass, NTs]

## Original RNP-MaP analysis ##
def RNP_MaP(NT_dfs, Zthresh=2.575, threshes=[0.29, 0.93, 0.78, 0.59], min_muts_above_bg=50):
	nts = ["A", "C", "G", "U"]
	nt_factors = {nt: thresh for nt, thresh in zip(nts, threshes)}

	new_dfs = {}
	NTthreshes = {}
	
	for NT in NT_dfs:
		new_dfs[NT] = NT_dfs[NT].copy()
		new_dfs[NT]['XL_Error'] = np.sqrt((new_dfs[NT]['XL_rate'])) / np.sqrt((new_dfs[NT]['XL_depth']))
		new_dfs[NT]['NoXL_Error'] = np.sqrt((new_dfs[NT]['NoXL_rate'])) / np.sqrt((new_dfs[NT]['NoXL_depth']))
		new_dfs[NT]['Zfactor'] = 1 - ((Zthresh * (new_dfs[NT]['XL_Error'] + new_dfs[NT]['NoXL_Error'])) / (new_dfs[NT]['XL_rate'] - new_dfs[NT]['NoXL_rate']).abs())
		new_dfs[NT]['Zfactor'] = new_dfs[NT]['Zfactor'].replace([np.inf, -np.inf], -100)
		new_dfs[NT]['RelDiff'] = new_dfs[NT]['XL_rate'] / new_dfs[NT]['NoXL_rate']
		new_dfs[NT]['mts_abv_bg'] = new_dfs[NT]['XL_mutations'] - (new_dfs[NT]['NoXL_rate'] * new_dfs[NT]['XL_depth'])

		SDs = new_dfs[NT]['RelDiff'].std(ddof=0)
		NTthreshes[NT] = new_dfs[NT]['RelDiff'].median() + nt_factors[NT] * SDs

		new_dfs[NT]['RNPsite'] = np.where(
			(new_dfs[NT]['mts_abv_bg'] >= min_muts_above_bg) &
			(new_dfs[NT]['Zfactor'] > 0) &
			(new_dfs[NT]['RelDiff'] > NTthreshes[NT]), 1, 0
		)

	return [new_dfs, NTthreshes]

## recombines results to produce final output ##
def recombine_results(RNP_dfs, AllInfo):
	results = pd.DataFrame()
	for NT in RNP_dfs:
		results = results.append(RNP_dfs[NT])
	results.sort_index(inplace=True)
	newcols = results[['XL_Error', 'NoXL_Error', 'Zfactor', 'RelDiff', 'RNPsite']]
	catted = pd.concat([AllInfo, newcols], axis=1)
	catted['RNPsite'] = catted['RNPsite'].fillna(0)
	return catted

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Runs legacy RNP-MaP analysis.",usage='%(prog)s [options] --SMtxt ',add_help=False)
	parser.add_argument("--SMtxt", help="Input shapemapper2 profile.txt with RNP-MaP data.", required=True)
	parser.add_argument("-o", "--OutPrefix", default=None, help="Specify Prefix of output files.")
	parser.add_argument('-d', "--LoDepth", type=int, default=5000, help="Specify depth threshold for inclusion.")
	parser.add_argument('-b', "--HiBG", type=float, default=0.05, help="Specify the background rate to exclude a nucleotide.")
	parser.add_argument('-z', "--zthresh", type=float, default=2.575, help="Specify Z threshold. Default is 2.575 (99% CL)")
	parser.add_argument('-h', '--help', action='help', help="show this help message and exit")

	args = parser.parse_args()

	prefix = cure_prefix(args.SMtxt, prefix=args.OutPrefix)
	mainf = Make_Dirs(prefix)

	nts = ["A", "C", "G", "U"]
	AllInfo, BG_PassInfo, NT_dfs = parse_SM(args.SMtxt, LoDepth=args.LoDepth, HiBG=args.HiBG, nts=nts)

	RNP_dfs, NTthreshes = RNP_MaP(NT_dfs,Zthresh=args.zthresh)
	result = recombine_results(RNP_dfs, AllInfo)

	results_file = "{}/{}_RESULTS".format(mainf,prefix)
	SMparse.write_RNP_csv(result,results_file)

	thresh_txt = "{}/{}_NT_thresholds.txt".format(mainf, prefix)
	with open(thresh_txt, "w") as threshfile:
		for NT in NTthreshes:
			threshfile.write("{}-nt threshold difference is {}\n".format(NT, NTthreshes[NT]))