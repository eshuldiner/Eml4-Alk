import numpy as np
import random
import sys
sys.path.append('/home/eshuldin/private_python/')

def LN_mean(t_sizes):
	LN_x = np.log(t_sizes)
	m = np.mean(LN_x)
	v = np.var(LN_x, ddof=1)
	output = np.exp(m + 0.5*v) 
	return(output)

def trimmed_mean(S, prop):
    from scipy.stats import trim_mean
    return trim_mean(S, prop) if prop > 0 else np.mean(S) #prop is proportion of values to trim for trimmed mean

def fdrCorrection(pval_dict, geneList=False):
	'''This function performs an FDR-adjustment of bootstrapped pvalues for an analysis where the input is structured as follows:
	sgID -> statisic -> pvalue. E.g., if you have calculated different percentiles.

	FDR correction is performed per statistic, across all sgids.

	If data is also being analyzed at the gene-level, a separate FDR-adjustment is performed across genes.

	'''
	import statsmodels.stats.multitest as smm
	if geneList:
		print("Doing separate FDR adjustment for sgIDs and genes")
		sgIDList = [i for i in pval_dict.keys() if i not in geneList]

		corrected_gene_level = {}
		corrected_sgid_level = {}

		stats =list(pval_dict[random.choice(list(pval_dict.keys()))].keys())
		for stat in stats:
			pvals_to_correct = [pval_dict[i][stat] for i in geneList]
			corrected_gene_level[stat] = list(smm.fdrcorrection(np.array([pval_dict[i][stat] for i in geneList])))[1]
			corrected_sgid_level[stat] = list(smm.fdrcorrection(np.array([pval_dict[i][stat] for i in sgIDList])))[1]

		p_sgid_out = {}
		p_gene_out = {}

		for i,gene in enumerate(geneList):
			p_gene_out[gene] = {}
			for stat in corrected_gene_level:
				p_gene_out[gene][stat] = corrected_gene_level[stat][i]

		for i,sgid in enumerate(sgIDList):
			p_sgid_out[sgid] = {}
			for stat in corrected_sgid_level:
				p_sgid_out[sgid][stat] = corrected_sgid_level[stat][i]				

		both_out = {**p_sgid_out, **p_gene_out}
		return(both_out)
	else:
		p_sgid_out = {}
		sgIDlist = list(pval_dict.keys())
		print("Only analyzing data the sgID-level; therefore not performing an FDR-adjustment for genes separately")
		corrected_sgid_level = {}
		stats =list(pval_dict[random.choice(sgIDlist)].keys())
		for stat in stats:
			corrected_sgid_level[stat] = list(smm.fdrcorrection(np.array([pval_dict[i][stat] for i in sgIDlist])))[1]
			print("corrected_sgid_level for stat {} is {}, len {}\n".format(stat, corrected_sgid_level[stat], len(corrected_sgid_level[stat])))

		for i,sgid in enumerate(sgIDlist):
			p_sgid_out[sgid] = {}
			for stat in corrected_sgid_level:
				p_sgid_out[sgid][stat] = corrected_sgid_level[stat][i]
		return(p_sgid_out)



def size_trim(mouse_dict, min_t_size):
	'''takes mouse_dict (result of input_data)
	and returns dictionary only containing tumors
	above a size threshold'''
	mouse_dict_trim = {}
	for sample in mouse_dict.keys():
		mouse_dict_trim[sample] = {}
		for sgid in mouse_dict[sample].keys():
			mouse_dict_trim[sample][sgid] = []
			for t in mouse_dict[sample][sgid]:
				if t > min_t_size:
					mouse_dict_trim[sample][sgid].append(t)
	return(mouse_dict_trim)

def input_data(infname, list_of_samples, samples_to_exclude, sgids_to_exclude):
	'''Function that reads in tumor file and a list of samples and 
	returns a nested dictionary where keys are 
	sample --> sgid --> list of tumor sizes (cell #)
	infname is tumor file
	'''

	list_of_samples_to_include = [i for i in list_of_samples if i not in samples_to_exclude]
	sgids_to_exclude.append("Spi")
	print("In input_data function, list of samples to include is {}, and sgids NOT being included are {}\n".format(list_of_samples_to_include, sgids_to_exclude))
	mouse_dict = {}
	f = open(infname, 'rt')
	f.readline()
	for l in f:
		fields = l.strip().split(",")
		mouse = fields[5]
		sgid = fields[0]
		if sgid not in sgids_to_exclude and mouse in list_of_samples_to_include:

			GT = fields[6]
			CN = float(fields[3])
			#print("tumor size is {}\n".format(CN))
			if mouse in mouse_dict.keys():
				if sgid in mouse_dict[mouse].keys():
					mouse_dict[mouse][sgid].append(CN)
				else:
					mouse_dict[mouse][sgid] = [CN]
			else:
				mouse_dict[mouse] = {}
				mouse_dict[mouse][sgid] = [CN]

	f.close()
		
	return(mouse_dict)


def merge_mice(mouse_dict):
	merged_mouse_dict = {}

	for mouse in mouse_dict.keys():
		for sgid in mouse_dict[mouse].keys():
			if sgid in merged_mouse_dict.keys():
				merged_mouse_dict[sgid].extend(mouse_dict[mouse][sgid])
			else:
				merged_mouse_dict[sgid] = mouse_dict[mouse][sgid]

	return(merged_mouse_dict)


def getP(stats, alternative=2, nullValue=1):
	if alternative == 2:
		n_boot = len(stats)
		boolean_stat_ge = []
		boolean_stat_le = []
		for i in stats:
			boolean_stat_ge.append(i>=nullValue)
			#boolean_stat_le.append(i<=nullValue)
		pval = 2 * min((sum(boolean_stat_ge) / n_boot), 1 - (sum(boolean_stat_ge) / n_boot))
		if pval > 1:
			sys.exit("ERROR! pval is {}. sum(boolean_stat_ge) is {}, nboot is {}\n".format(pval, sum(boolean_stat_ge), n_boot))

		return(pval)
	elif alternative == "greater":
		# arbitrarily testing if X >= 1
		n_boot = len(stats)
		boolean_stat_ge = []
		for i in stats:
			boolean_stat_ge.append(i >= nullValue)
		return(1 - sum(boolean_stat_ge) / n_boot)
	elif alternative == "less":
		n_boot = len(stats)
		boolean_stat_le = []
		for i in stats:
			boolean_stat_le.append(i <= nullValue)
		return(1 - sum(boolean_stat_le) / n_boot)	
	else:
		sys.exit("Invalid alternative specified, {}\n".format(alternative))	



def adj_percentiles(data, percentiles, inerts, genelevel = False):
	'''April 27 2023, finally adding error catching: should
	throw an error telling you what sgid it failed on

	Data here is already merged across mice. So it is just dict of sgid --> lists of CellNums

	First, collect information for your inert sgRNAs. Then, normalize everything to either weighted mean or inerts
	or median value of sgRNAs.
	'''
	from scipy.stats.mstats import gmean
	#print("in adj percentiles, inerts are {}\n".format(inerts))

	perc_out = {}
	inert_tumors = []

	inert_percs={}
	inert_medians={}

	for p in percentiles:
		perc=p*100 
		inert_percs[p] = []

		for i in inerts:
			if i in data.keys():
				i_tumors = data[i]
				inert_percs[p].append(np.percentile(i_tumors,perc))

	inert_percs["LNmean"] = []

	for i in inerts:
		if i in data.keys():
			i_tumors = data[i]
			inert_percs["LNmean"].append(LN_mean(i_tumors))
	for p in inert_percs.keys():
		inert_medians[p] = np.median(inert_percs[p])

	for sgid in data.keys():
		try:
			perc_out[sgid] = {}
			tumors = data[sgid]
			for p in percentiles:
				perc = p * 100
				perc_out[sgid][p] = np.percentile(tumors, perc) / inert_medians[p]

			perc_out[sgid]["LNmean"] = LN_mean(tumors) / inert_medians["LNmean"]
		except:
			print("Error on sgid {}. This sgid has {} tumors and there are {} inert tumors.\n".format(sgid, len(tumors), len(inert_tumors)))


	if genelevel:
		#print("In Adj percentiles, I recognize gene levels")
		gene_out = {}
		both_out = {}

		TN_num_sgID = {k : len(v) for (k,v) in data.items()}
		genes = list(set([sgRNA_to_gene(i) for i in TN_num_sgID]))
		TN_num_gene = {g:0 for g in genes}
		TN_gene_to_sgid_weight = {g:{} for g in genes}

		for sgid,count in TN_num_sgID.items():
			TN_num_gene[sgRNA_to_gene(sgid)] += count

		for sgid, count in TN_num_sgID.items():
			gene = sgRNA_to_gene(sgid)
			TN_gene_to_sgid_weight[gene][sgid] = count / TN_num_gene[gene]

		for gene in TN_gene_to_sgid_weight.keys():
			gene_out[gene] = {}
			sgids = list(TN_gene_to_sgid_weight[gene].keys())
			#print("sgids are {}\n".format(sgids))
			weights = [TN_gene_to_sgid_weight[gene][s] for s in sgids] # vector of weights
			try:
				for p in perc_out[sgids[0]].keys():
					vals = [perc_out[i][p] for i in sgids]
					gene_out[gene][p]=sum(x*y for x, y in list(zip(weights, vals)))
			except:
				print("Error on gene {}, which is associated with sgids {}\n".format(gene, sgids))
		both_out = {**perc_out, **gene_out}
		return(both_out)
	else:
		return(perc_out)

def sgRNA_to_gene(sgRNA):
	pre = sgRNA.split("_")[0]
	if pre in ["NT","Neo","Safe","V3.NT","V3.Neo","V3.Safe","V3.sgNeo","V3.sgNT","V3.sgSafe"]:
		return("Inert")
	elif pre in ["V1.NT","V1.Neo","V1.Safe","V1.sgNeo","V1.sgNT","V1.sgSafe"]:
		return("V1Inert")
	else:
		return(pre)




def pull_focal_stat(multi_stat_dataset, focal_stat):
	'''Input: dictionary where keys are sgID --> stat --> values.
	   Output: dictionary where keys are sgID --> value

	   Helper function for rank_calculation
	'''
	single_stat_dataset = {}
	for s in multi_stat_dataset:
		print("{}, {}\n".format(s, multi_stat_dataset[s]))
	for sgID in multi_stat_dataset:
		single_stat_dataset[sgID] = multi_stat_dataset[sgID][focal_stat]
	return(single_stat_dataset)


def rank_calculation(dataset, exclude=[], include=False):
	'''
	Input: dictionary where keys are sgIDs and values are statistics.
	Output: dictionary where keys are sgIDs and values are ranks (where 1 = strongest effect)

	Can optionally only rank certain sgids, as specified using exclude and include parameters
	'''

	dataset_filt = {}

	if not include:
		for sgid in dataset:
			if sgid not in exclude:
				dataset_filt[sgid] = dataset[sgid]
	else:
		for sgid in dataset:
			if sgid not in exclude:
				if sgid in include:
					dataset_filt[sgid] = dataset[sgid]
		
	return({key: rank for rank, key in enumerate(sorted(dataset_filt, key=dataset_filt.get, reverse=True), 1)})


def bootstrap_mice(mouse_dict, n_mice):
	'''This function resamples <n_mice> mice from dictionary where the keys are mice
	and the values are all the tumors in the mice (keyed sgID --> list of tumor sizes)'''
	bs_mouse_dict = {}
	mouse_list = list(mouse_dict.keys())
	bs_mouse_list = np.random.choice(mouse_list,n_mice, replace = True)

	for i,m in enumerate(bs_mouse_list):
		bs_mouse_dict[i] = mouse_dict[m]

	return(bs_mouse_dict)


def bootstrap_tumors(data):
	'''Bootstrap tumors within each mouse by sgID
	'''
	new_data = {}
	for m in data.keys(): #for each mouse
		new_data[m] = {}
		for s in data[m].keys():
			t_indices = range(len(data[m][s]))
			new_t_indices = np.random.choice(t_indices, len(t_indices), replace = True)
			new_t = [data[m][s][i] for i in new_t_indices]
			new_data[m][s] = new_t
	return(new_data)
