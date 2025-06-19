import numpy as np
import random
import sys
from scipy.stats import spearmanr, ttest_ind

sys.path.append('/home/eshuldin/private_python/')


def randomRankingSimulation(nGenes):
	from scipy.stats import spearmanr

	x = list(range(1, nGenes+1))
	y = random.sample(x, k=len(x))
	return(spearmanr(np.array(x), np.array(y)).correlation)


def calculateRankStat(ranks1, ranks2):

	ranks_numerical1 = []
	ranks_numerical2 = []
	for k,v in ranks1.items():
		ranks_numerical1.append(v)
		ranks_numerical2.append(ranks2[k])
	spearmanobj = spearmanr(np.array(ranks_numerical1), np.array(ranks_numerical2))
	return((spearmanobj.correlation,spearmanobj.pvalue))









