import numpy as np
import getopt
import pandas as pd
import sys
from subprocess import call
from percentile_helper_functions import getP,LN_mean

def getNetTiters(sampleInfo, FiltCol="TREATMENT",filtVal=["Cas9neg","CAS9-NEGATIVE","K","KP"]):
	'''Having read in the sample info files, return the net titer for each experiment, splititng on the specified treatment
	'''
	sampleInfodict = {key: group for key, group in sampleInfo.groupby('DATASET')}
	data = {}
	cas9neg = {}
	for key, value in sampleInfodict.items():

		main = value[~value[FiltCol].isin(filtVal)] #value[FiltCol] is the column of the df
		filt = value[value[FiltCol].isin(filtVal)]
		data[key] = main["TITERS"].sum()
		cas9neg[key] = filt["TITERS"].sum()
	return(data,cas9neg)


def input_data_adaptive_genotype_sort(infname, sampleInfo, sgids_to_exclude, sortOn=["Cas9neg","V1_Cas9neg","V3_Cas9neg","KT","K","KP","RIT1","Rit1","CAS9-NEGATIVE"], convert_sgRNA=False):
	'''Function that reads in tumor file and a list of sgids and samples to exclude and 
	returns a nested dictionary where keys are 
	sample --> sgRNA --> list of tumor sizes (cell #)
	infname is tumor file

	Returns two dictionaries
	'''
	#print("in input data, sgids to exclude is {}, samples to exclude are {} and RCtrim is {}, filterOn is {}\n".format(sgids_to_exclude,samples_to_exclude,RCtrim, filterOn))
	sgids_to_exclude.append("Spi") #always exclude reads from spike-in cells
	mouse_dict = {}
	mouse_dict2 = {}
	f = open(infname, 'rt')
	header = f.readline().strip().split(",")
	samples = list(sampleInfo['SAMPLES'].unique())


	for l in f:
		fields = l.strip().split(",")
		mouse = fields[5]
		sgid = fields[0]
		RC = int(fields[2])
		CN = float(fields[3])
		if sgid not in sgids_to_exclude and sgid_to_sgRNA(sgid) not in sgids_to_exclude:
			if mouse in samples:
				if convert_sgRNA:
					sgRNA = sgid_to_sgRNA(sgid)
				else:
					sgRNA = sgid

				GT = fields[6]
				if GT in sortOn: #this is data from a mouse that goes into second dict
					if mouse in mouse_dict2.keys():
						if sgRNA in mouse_dict2[mouse].keys():
							mouse_dict2[mouse][sgRNA].append(CN)
						else:
							mouse_dict2[mouse][sgRNA] = [CN]
					else:
						mouse_dict2[mouse] = {}
						mouse_dict2[mouse][sgRNA] = [CN]
				else: #this is data from a mouse that goes into first dict
					if mouse in mouse_dict.keys():
						if sgRNA in mouse_dict[mouse].keys():
							mouse_dict[mouse][sgRNA].append(CN)
						else:
							mouse_dict[mouse][sgRNA] = [CN]
					else:
						mouse_dict[mouse] = {}
						mouse_dict[mouse][sgRNA] = [CN]

	f.close()
		
	return(mouse_dict, mouse_dict2)


def get_sgID_ratio(dataset, inerts):
	''' Function that returns the breakdown of the viral pool
	Inputs:
		dataset: mouse -> sgID -> cell numbers. Should be from cas9-negative mice.
		inerts: list of inert sgIDs (because I calculated % across all inerts)

	Outputs:
		sample (or median) -> sgID (or inerts) -> proportion
	'''
	ratio = {}
	sgids=[]
	mice = []
	ratio["Median"]={}
	for m in dataset:
		inert_total = 0
		mice.append(m)
		ratio[m] = {}
		totalTN = sum([len(s) for s in dataset[m].values()])
		for s in dataset[m]:
			sgids.append(s)
			ratio[m][s] = len(dataset[m][s])/totalTN
		for i in inerts:
			try:
				inert_total += len(dataset[m][i])
			except:
				print("Possible ienrt sgID {} not found in dataset\n".format(i))
				pass
		ratio[m]["Inert"] = inert_total/totalTN


	sgids_final = list(set(sgids))
	if "Spi" in sgids_final:
		sgids_final.remove("Spi")
	sgids_final.append("Inert")

	for s in sgids_final:
		vals = []
		for m in mice:
			try:
				vals.append(ratio[m][s])
			except:
				print("I am appending 0 for sgid {} for mouse {}\n".format(s,m))
				vals.append(0)
		ratio["Median"][s] = np.median(vals)

	return(ratio)

def get_adaptive_dataset_sizes(datasets, focal_dataset, focal_sgID, focal_number, titer_dict, sgID_ratios):
	'''Calculates the # of tumors that should be sampled for each dataset based on
	1. # of tumors sampled for focal sgID in focal dataset
	2. titers used in all datasets
	3. sgID ratios in all datasets
	Inputs:
		datasets: vector of dataset IDs
		focal_dataset: dataset ID of focal dataset
		focal_sgID: focal sgID (all sampling relative to this in the focal dataset)
		focal_number: # of tumors to sample for focal sgID in focal dataset
		titer_dict: datasetID --> overall titer
		sgID_ratios: datasetID --> mouse/median --> sgID --> proportion in pool
	Outputs:
		datasetID -> sgID --> # of tumors to sample
	'''
	adaptive_dataset_sizes = {}
	adaptive_dataset_sizes[focal_dataset] = {}
	adaptive_dataset_sizes[focal_dataset][focal_sgID] = focal_number

	non_focal_datasets = [i for i in datasets if i !=focal_dataset]

	for d in non_focal_datasets:
		adaptive_dataset_sizes[d] = {}
		adaptive_dataset_sizes[d][focal_sgID] = focal_number * titer_dict[d]/titer_dict[focal_dataset] * sgID_ratios[d][focal_sgID]/sgID_ratios[focal_dataset][focal_sgID]

	for d in sgID_ratios:
		for s in sgID_ratios[d]:
			if s != focal_sgID:
				adaptive_dataset_sizes[d][s] = adaptive_dataset_sizes[d][focal_sgID] * sgID_ratios[d][s]/sgID_ratios[d][focal_sgID]
	return(adaptive_dataset_sizes)

def sgid_to_sgRNA(string):
	splitter = string.strip().split("_")
	return("_".join([splitter[0],splitter[-1]]))

def take_top_N_weighted(tumor_dict, adaptive_dataset_size):
	'''Takes the top N tumors per cohort by dividing N among all mice in the cohort
	based on Tumor Number (so weighting the contributions of each mouse based on tumor number)
	'''
	###Figure out how you're going to distribute the N tumors across each mouse
	count_dict= {}
	ratio_dict = {}
	TTN = 0
	for m in tumor_dict:
		nT = len(sum(list(tumor_dict[m].values()), []))
		count_dict[m] = nT
		TTN += nT
	for m in count_dict:
		ratio_dict[m] = count_dict[m]/TTN

	#########
	top_dict = {}
	T_sampled = {}
	for sgID in adaptive_dataset_size:
		T_sampled[sgID] = 0
	
	for m in tumor_dict: # looping over mice
		top_dict[m] = {}
		min_t_size=min(sum(list(tumor_dict[m].values()), [])) # minimum tumor size for this sample
		for s in tumor_dict[m]: #looping over sgIDs
			tumors_all = sorted(tumor_dict[m][s], reverse=True) # per mouse
			overallN = int(adaptive_dataset_size[s])
			mouseN = round(overallN * ratio_dict[m])
			T_sampled[s] += mouseN
			if mouseN <=len(tumors_all): #if there are enough tumors for you to take the top N...
				top_dict[m][s] = tumors_all[:mouseN] # take first N items from list (note that list has been sorted by size)
			else:
				n_extra = mouseN-len(tumors_all)
				padding = [min_t_size]*n_extra
				for_dict = tumors_all[:mouseN] + padding
				top_dict[m][s] = for_dict

	return(top_dict)

def percentile_bootstrapping_preparation(size_point_estimates):
	'''Prepapre data structures to store information from bootstrapping.

	Note that this code is written to be flexible based upon the formating of
	the input dictionary, size_point_estimates. Can be sgID -> statistic -> value
	or just sgID -> value
	'''
	bootstrapped_percentiles = {}
	CI_upper = {}
	CI_lower = {}
	P = {}
	try: # additional layer of keying beyond sgID
		for sgid in size_point_estimates.keys():
			bootstrapped_percentiles[sgid] = {}
			CI_upper[sgid] = {}
			CI_lower[sgid] = {}
			P[sgid] = {}
			for p in size_point_estimates[sgid].keys():
				bootstrapped_percentiles[sgid][p] = []
				CI_upper[sgid][p] = 0
				CI_lower[sgid][p] = 0
				P[sgid][p] = 0
			bootstrapped_percentiles[sgid]["LNmean"] = []
			CI_upper[sgid]["LNmean"] = 0
			CI_lower[sgid]["LNmean"] = 0
			P[sgid]["LNmean"] = 0
	except: #Only one value per sgID
		for sgid in size_point_estimates.keys():
			bootstrapped_percentiles[sgid] = []
			CI_upper[sgid] = 0
			CI_lower[sgid] = 0
			P[sgid] = 0
	return(bootstrapped_percentiles,CI_upper,CI_lower,P)

def percentile_bootstrap_resulting_processing(bootstrapped_percentiles, CI_upper, CI_lower, P, upper, lower):
	'''Also hopefully written flexibly
	'''
	try:
		for sgid in bootstrapped_percentiles.keys():
			for p in bootstrapped_percentiles[sgid].keys():
				CI_upper[sgid][p] = np.percentile(bootstrapped_percentiles[sgid][p], upper)
				CI_lower[sgid][p] =np.percentile(bootstrapped_percentiles[sgid][p], lower)
				P[sgid][p] = getP(bootstrapped_percentiles[sgid][p])
	except:
		for sgid in bootstrapped_percentiles.keys():
			CI_upper[sgid] = np.percentile(bootstrapped_percentiles[sgid], upper)
			CI_lower[sgid] =np.percentile(bootstrapped_percentiles[sgid], lower)
			P[sgid] = getP(bootstrapped_percentiles[sgid])	
	return(CI_upper,CI_lower,P)

