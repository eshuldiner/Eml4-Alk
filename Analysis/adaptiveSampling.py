import getopt
import sys
import os
import numbers
import numpy as np
import pandas as pd
#from scipy.stats.mstats import gmean
from subprocess import call
from helper_functions import read_project_info, read_parameter_info
from percentile_helper_functions import merge_mice, adj_percentiles, bootstrap_mice, bootstrap_tumors, getP, LN_mean, pull_focal_stat, rank_calculation, fdrCorrection
from adaptive_helper_functions import input_data_adaptive_genotype_sort, get_sgID_ratio, get_adaptive_dataset_sizes, take_top_N_weighted, percentile_bootstrapping_preparation, percentile_bootstrap_resulting_processing, getNetTiters
from ranking_helper_functions import randomRankingSimulation, calculateRankStat

# Step 0.0: collect arguments from the command line

# argument defaults
datasets=None
percentiles = [0.5,0.6,0.7,0.75,0.8,0.9,0.95]
upper = 97.5
lower = 2.5
focal_number = None
focal_sgID = None
focal_dataset = None
sgID_ratios = None
analysis_type="None"
root = "~/"
incList=False
includeStr="All"
exclude=[]
mainNetT = {}
cas9negNetT = {}
geneLevel=False
rankOn=0.95
sgids_to_exclude_pre=""
sgids_to_exclude=[]

try:
    opts, args = getopt.getopt(sys.argv[1:], "", ["root=",
        "analysis=", 
        "project_name=",
        "info_file=",
        "sgids_to_exclude=",
        "focal_dataset=",
        "focal_sgID=",
        "focal_number=",
        "datasets=",
        "datafiles=",
        "inerts=",
        "nboot=",
        "incList=",
        "exclude=",
        "geneLevel=",
        "rankOn="])
except getopt.GetoptError:
    print("no arguments recognized\n")
    sys.exit(2)
for opt, arg in opts:
    if opt in ("--root"):
        root = arg
        print("found root, {}\n".format(root))
    elif opt in ("--project_name"):
        project_id = arg
        print("found project_name, {}\n".format(project_id))
    elif opt in ("--analysis"):
        analysis_type = arg
        print("found analysis, {}\n".format(analysis_type))
    elif opt in ("--info_file"):
        infofile=arg.strip()
    elif opt in ("--sgids_to_exclude"):
    	sgids_to_exclude_pre=arg
    elif opt in ("--focal_dataset"):
    	focal_dataset = arg
    elif opt in ("--focal_sgID"):
    	focal_sgID = arg
    elif opt in ("--focal_number"):
    	focal_number = int(arg)
    elif opt in ("--datasets"):
    	datasets = arg.split(",")
    elif opt in ("--datafiles"):
    	datafiles=arg.strip().split(",")
    elif opt in ("--inerts"):
        inerts=[str(x) for x in arg.strip().split(",")]
    elif opt in ("--nboot"):
    	n_boot = int(arg)
    elif opt in ("--incList"):
        splitter = arg.strip().split(":")
        print("arg is {}, splitter is {}\n".format(arg,splitter))
        includeStr = splitter[0]
        incList = [str(x) for x in splitter[1].strip().split(",")]
    elif opt in ("--exclude"):
        exclude=[str(x) for x in arg.strip().split(",")]
    elif opt in ("--geneLevel"):
        geneLevel=True
    elif opt in ("--rankOn"):
        try:
            rankOn = float(arg)
        except:
            rankOn=str(arg)
    else:
        assert False, "unhandled option"

################################################################################################################################################
# Step 0.1 Process arguments
if datasets == None:
	sys.exit("ERROR, datasets not provided\n")


if ":" in sgids_to_exclude_pre:
	sgids_to_exclude_lists = sgids_to_exclude_pre.strip().split(":")
	if len(sgids_to_exclude_lists) != len(datasets):
		sys.exit("Error, # of sgid exclusions lists {} does not match # of datasets {}\n".format(len(sgids_to_exclude_lists), len(datasets)))
	else:
		sgids_to_exclude = dict(map(lambda i,j : (i,j.strip().split(",")) , datasets, sgids_to_exclude_lists))
else:
	sgids_to_exclude = dict(map(lambda i,j : (i,j) , datasets, [sgids_to_exclude_pre.strip().split(",")]*len(datasets)))


print("datasets are {}\nsgIDs being excluded from analysis are {}\n".format(datasets, sgids_to_exclude))

print("project_id is {}, focal_dataset {}, focal_sgID {}, focal_number is {}\n".format(project_id, focal_dataset, focal_sgID, focal_number))

#Ensure output path exists
if not os.path.exists(root):
    os.makedirs(root)
    print(f"Directory created: {root}")

################################################################################################################################################
# Step 0.2 Make sure that sampling method is appropriately specified and that proper information to implement it was provided.

if [x for x in (focal_number, focal_sgID, focal_dataset) if x is None]:
    sys.exit("You must explicitly specific a focal number, focal sgID and focal dataset. One of these is currently None\n")
if not isinstance(focal_number, numbers.Number):
    sys.exit("Focal number is not correctly specified {}\n".format(focal_number))

parameter_summary_str = project_id + "_" + focal_dataset + "_" + focal_sgID + "_" + str(focal_number)

sampleInfo = pd.read_csv(infofile, sep='\t')

mainNetT,cas9negNetT = getNetTiters(sampleInfo)


################################################################################################################################################
#Step 1.1 and 1.2 Input raw data, calculate sgID ratios for each dataset
#INPUT data and do necessary pre-calculations
all_data = {} 
all_data_cas9neg = {}
sgID_ratios = {}


for i,d in enumerate(datasets):
    print("Loading data for dataset {}\n".format(d))
    x,xn = input_data_adaptive_genotype_sort(infname = datafiles[i], sampleInfo = sampleInfo, sgids_to_exclude = sgids_to_exclude[d], convert_sgRNA=True)
    all_data[d] = x
    all_data_cas9neg[d] = xn

    full=get_sgID_ratio(xn, inerts)
    sgID_ratios[d] = full["Median"]


################################################################################################################################################

def singleDatasetAnalysis(data, adaptive_dataset_size=None, genelevel=False):
    '''For a given dataset: adaptively samples tumors, merges data across mice, and calculates
    the percentile tumor sizes for that dataset.
    '''
    adaptive_d = take_top_N_weighted(data, adaptive_dataset_size)

    if any(isinstance(i,dict) for i in adaptive_d.values()):
        merged_mouse_dict = merge_mice(adaptive_d)

    else:
        merged_mouse_dict= adaptive_d

    size_point_estimates = adj_percentiles(merged_mouse_dict, percentiles, inerts, genelevel)

    return(size_point_estimates, merged_mouse_dict)

def analysis_rank_stat(root, focal_number, datasets, sgID_ratios, net_titer_dict, focal_dataset, focal_sgID, focalStat, nboot, summary_str, includeStr, exclude, GL, incList=False):
    '''Currently can only rank based on percentile values (or LNmean)
    focalStat is the percentile you want to rank on.

    Iterates through pairs of datasets to do the comparisons
    '''
    rankstat_outname = root + parameter_summary_str + "_" + includeStr + "_" + str(focalStat) + "_adaptive_sampling_testing_ranking.txt"
    rank_out = open(rankstat_outname, 'wt')
    rank_out.write("Dataset1,Dataset2,IncludeList,Focal_N,focal_sgID,RankingOn,nboot,point_est,CI_median,CI_lower,CI_upper,P_pe\n")

    adaptive_dataset_sizes=get_adaptive_dataset_sizes(datasets, focal_dataset, focal_sgID, focal_number, net_titer_dict, sgID_ratios)

    dataset_list = list(all_data.keys())


    for i,dataset_id1 in enumerate(dataset_list):
        for j,dataset_id2 in enumerate(dataset_list):
            if j>=i:
                print("In ranking comparison analysis, working on datasets {} and {}\n".format(dataset_id1, dataset_id2))
                data1 = datasets[dataset_id1]
                data2 = datasets[dataset_id2]
                adaptive_sizes1 = adaptive_dataset_sizes[dataset_id1]
                adaptive_sizes2 = adaptive_dataset_sizes[dataset_id2]
                pe_results1, sampleddata1 = singleDatasetAnalysis(data1, adaptive_sizes1, GL)
                pe_results2, sampleddata2= singleDatasetAnalysis(data2, adaptive_sizes2, GL)

                pe_ranks1 = rank_calculation(pull_focal_stat(pe_results1, focalStat), exclude=exclude, include=incList)
                pe_ranks2 = rank_calculation(pull_focal_stat(pe_results2, focalStat), exclude=exclude, include=incList)

                pe_RankStat, p_RankStatPE = calculateRankStat(pe_ranks1, pe_ranks2)

                bootstrapped_stats = []
                CI_upper=0
                CI_lower=0
                for j in range(nboot):
                    if j%25==0:
                        print("working on bootstrap {}\n".format(j))
                    bs_mice1 = bootstrap_mice(data1,len(data1.keys()))
                    bs_mice2 = bootstrap_mice(data2,len(data2.keys()))
                    bs_mice_bs_tumors1 = bootstrap_tumors(bs_mice1)
                    bs_mice_bs_tumors2 = bootstrap_tumors(bs_mice2)
                    bs_results1, bs_samp_data1 = singleDatasetAnalysis(bs_mice_bs_tumors1, adaptive_sizes1, GL)
                    bs_results2, bs_samp_data2 = singleDatasetAnalysis(bs_mice_bs_tumors2, adaptive_sizes2, GL)
                    bs_ranks1 = rank_calculation(pull_focal_stat(bs_results1, focalStat), exclude=exclude, include=incList)
                    bs_ranks2 = rank_calculation(pull_focal_stat(bs_results2, focalStat), exclude=exclude, include=incList)
                    bs_RankStat, bs_RankStatPE = calculateRankStat(bs_ranks1, bs_ranks2)
                    bootstrapped_stats.append(bs_RankStat)
                CI_upper = np.percentile(bootstrapped_stats, upper)
                CI_lower = np.percentile(bootstrapped_stats, lower)
                med = np.percentile(bootstrapped_stats, 50)


                rank_out.write("{},{},{},{},{},{},{},{},{},{},{},{}\n".format(dataset_id1, dataset_id2, includeStr, focal_number, focal_sgID, focalStat, nboot, pe_RankStat, med, CI_upper, CI_lower, p_RankStatPE))

    Random_results=[randomRankingSimulation(nGenes=len(incList)) for _ in range(nboot)]
    random_median = np.percentile(Random_results, 50)
    random_upper = np.percentile(Random_results, 97.5)
    random_lower = np.percentile(Random_results, 2.5)
    rank_out.write("RandomSim,RandomSim,{},{},{},{},{},{},{},{},{},{}\n".format(includeStr, focal_number, focal_sgID, focalStat, nboot, random_median, random_median, random_upper, random_lower, "NA"))

    rank_out.close() 


def analysis(root, focal_number, datasets, sgID_ratios, net_titer_dict, focal_dataset, focal_sgID, nboot, summary_str, GL):
    print("in analysis")
    print("Datasets are {}\n".format(datasets.keys()))
    outname = root + parameter_summary_str + "_adaptive_sampling_percentiles.txt"
    outfname = open(outname, 'wt')
    outfname.write("Dataset,Focal_N,focal_sgID,N_sampled,nboot,sgID,percentile,datatype,point_est,CI_lower,CI_upper,P,P_FDR_corr\n")

    adaptive_dataset_sizes=get_adaptive_dataset_sizes(datasets, focal_dataset, focal_sgID, focal_number, net_titer_dict, sgID_ratios)

    
    for dataset_id in all_data:
        print("in analysis function, working on dataset {}\n".format(dataset_id))
        data = datasets[dataset_id]
        adaptive_sizes = adaptive_dataset_sizes[dataset_id]

        pe_results, sampled_data = singleDatasetAnalysis(data, adaptive_sizes, GL)
        
        bootstrapped_percentiles,CI_upper,CI_lower,P = percentile_bootstrapping_preparation(pe_results)

        for j in range(nboot):
            if j%25==0:
                print("working on bootstrap {}\n".format(j))
            bs_mice = bootstrap_mice(data,len(data.keys()))
            bs_mice_bs_tumors = bootstrap_tumors(bs_mice)
            bs_results,bs_sampled_data=singleDatasetAnalysis(bs_mice_bs_tumors, adaptive_sizes, GL)


            for sgid in bs_results.keys():
                for p in bs_results[sgid].keys():
                    bootstrapped_percentiles[sgid][p].append(bs_results[sgid][p])
                    bootstrapped_percentiles[sgid]["LNmean"].append(bs_results[sgid]["LNmean"])

        CI_upper,CI_lower,P = percentile_bootstrap_resulting_processing(bootstrapped_percentiles, CI_upper, CI_lower, P, upper, lower)
        if GL:
            test_gene_list = [i for i in list(P.keys()) if "_" not in i]
            P_FDR_corr=fdrCorrection(P, geneList=test_gene_list)
        else:
            P_FDR_corr=fdrCorrection(P)

        for sgid in pe_results.keys():
            for p in pe_results[sgid].keys():
                try:
                    outfname.write("{},{},{},{},{},{},{},Size,{},{},{},{},{}\n".format(dataset_id, focal_number,focal_sgID, adaptive_sizes[sgid], nboot, sgid, p, pe_results[sgid][p], CI_lower[sgid][p], CI_upper[sgid][p], P[sgid][p], P_FDR_corr[sgid][p]))
                except:
                    outfname.write("{},{},{},{},{},{},{},Size,{},{},{},{},{}\n".format(dataset_id, focal_number,focal_sgID, "NA-gene", nboot, sgid, p, pe_results[sgid][p], CI_lower[sgid][p], CI_upper[sgid][p], P[sgid][p], P_FDR_corr[sgid][p]))
    outfname.close()


if analysis_type == "basic":
    analysis(root=root, focal_number=focal_number, datasets=all_data, sgID_ratios=sgID_ratios, net_titer_dict=mainNetT, focal_dataset=focal_dataset, focal_sgID=focal_sgID, nboot=n_boot, summary_str=parameter_summary_str, GL=geneLevel)
elif analysis_type == "rank_stat":
    analysis_rank_stat(root=root, focal_number=focal_number, datasets=all_data, sgID_ratios=sgID_ratios, net_titer_dict=mainNetT, focal_dataset=focal_dataset, focal_sgID=focal_sgID, focalStat=rankOn, nboot=n_boot, summary_str=parameter_summary_str, exclude=exclude, includeStr=includeStr, GL=geneLevel, incList=incList)
else:
    sys.exit("Analysis type specified ({}) is invalid\n".format(analysis_type))



