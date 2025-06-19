import numpy as np
import sys

def unique(list1):
    x = np.array(list1)
    return(np.unique(x))

def revcom(text):
    complement={'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    return ''.join(complement[base] for base in text[::-1])


def hamming_distance(str1, str2):
    if len(str1)!=len(str2):
        print("ALERT! strings must be of same length to calculate hamming distance, but str1 has length {} and str2 has length {}\n".format(len(str1), len(str2)))
    return sum(str1[i]!=str2[i] for i in range(len(str1)))

def get_quality_string(sample_str):
    list_of_quals = [int(ord(c))-33 for c in sample_str]
    return list_of_quals

def get_mean_quality_read(list_of_quals):
    mean_qual = np.average(list_of_quals)
    return mean_qual

def get_gc(seq):
    A=float(seq.count('A')/len(seq))
    G=float(seq.count('G')/len(seq))
    T=float(seq.count('T')/len(seq))
    C=float(seq.count('C')/len(seq))
    GC=float((seq.count('G')+seq.count('C'))/len(seq))
    #return([A,T,G,C,GC])
    return(GC)

def read_project_info(project_info_file):
    print("in read_project_info function. Project_info_file is {}\n".format(project_info_file))
    class project_info():
        def __init__(self, _project_name = None, _project_type = None, _sample_ids = None, _sample_paths = None, _sample_to_sample_paths = None, _genotypes = None, _sample_to_gt = None,  _sample_to_titer = None, _control_gt = None, 
            _sgids = None, _genes = None, _sgRNAs = None, _inerts = None, _sgids_to_gene = None,
            _sample_to_spike_bcs = None, _sample_to_spike_sizes = {},
            _sample_to_treatment = None):
            
            self.project_name = _project_name
            self.project_type = _project_type

            self.sample_ids = _sample_ids
            self.sample_paths = _sample_paths
            self.sample_to_sample_paths = _sample_to_sample_paths
            self.genotypes = _genotypes
            self.sample_to_gt = _sample_to_gt
            self.sample_to_titer = _sample_to_titer
            self.control_gt = _control_gt
            self.sample_to_treatment = _sample_to_treatment
            self.sgids = _sgids
            self.genes = _genes
            self.sgRNAs = _sgRNAs
            self.inerts = _inerts
            self.sgids_to_gene = _sgids_to_gene
            self.sample_to_spike_bcs = _sample_to_spike_bcs
            self.sample_to_spike_sizes = _sample_to_spike_sizes

            self.populate_project_info(project_info_file)

        def populate_project_info(self, project_info_file):
            project_info={}
            f = open(project_info_file, 'rt')
            for l in f:

                if l[0]!="#":

                    fields = l.strip().split("\t")
                    print("fields is {}\n".format(fields))
                    project_info[fields[0]]=fields[1]

            self.project_name = project_info["PROJECT"].strip()
            self.sample_ids = project_info["SAMPLES"].strip().split(",")
            self.sample_paths = project_info["SAMPLES_PATHS"].strip().split(",")
            
            genotypes_raw = project_info["GENOTYPES"].strip().split(",")
            treatments_raw = project_info["TREATMENT"].strip().split(",")

            if "," in project_info["SIZE_SPIKE"]:
                spike_sizes_raw = project_info["SIZE_SPIKE"].strip().split(",")
            else:
                spike_sizes_raw = [project_info["SIZE_SPIKE"].strip()] * len(genotypes_raw)

            if "," in project_info["SPIKE_BC"]:
                spike_bc_raw = project_info["SPIKE_BC"].strip().split(",")
            else:
                spike_bc_raw = [project_info["SPIKE_BC"].strip()] * len(genotypes_raw)
            titers_raw = project_info["TITERS"].strip().split(",")


            self.control_gt = project_info["CAS9NEG_GT"].strip()
            self.sample_to_sample_paths = dict(zip(self.sample_ids, self.sample_paths))
            self.sample_to_gt = dict(zip(self.sample_ids, genotypes_raw))
            self.sample_to_titer = dict(zip(self.sample_ids, titers_raw))
            self.sample_to_treatment = dict(zip(self.sample_ids, treatments_raw))
            self.sample_to_spike_bcs = dict(zip(self.sample_ids, [i.strip().split(":") for i in spike_bc_raw]))
            self.sample_to_spike_sizes = dict(zip(self.sample_ids, [list(map(int, s.split(':'))) for s in spike_sizes_raw]))

            self.genotypes = unique(genotypes_raw)
            self.sgids = project_info["SGIDS"].strip().split(",")
            self.sgRNAs = project_info["SGRNAS"].strip().split(",")
            genes_raw = project_info["GENES"].strip().split(",")

            self.sgids_to_gene = {}
            for i in range(len(self.sgids)):
                self.sgids_to_gene[self.sgRNAs[i]] = genes_raw[i]
            self.genes = unique(genes_raw)

            self.inerts = project_info["INERT"].strip().split(",")


    return(project_info())


def read_parameter_info(parameter_info_file):
    class parameter_info():
        def __init__(self, _parameter_set_name = None, _dist_between_bc_for_read = None, _R1_regex_looseness = None, _R2_regex_looseness = None, _read_error_size_threshold = None, _bc_similarity_threshold = None,
            _contam_removal_threshold = None):
            self.dist_between_bc_for_read = _dist_between_bc_for_read
            self.parameter_set_name = _parameter_set_name
            self.R1_regex_looseness = _R1_regex_looseness
            self.R2_regex_looseness = _R2_regex_looseness
            self.read_error_size_threshold = _read_error_size_threshold
            self.bc_similarity_threshold = _bc_similarity_threshold 
            self.contam_removal_threshold = _contam_removal_threshold
            self.populate_parameter_info(parameter_info_file)

        def populate_parameter_info(self, parameter_info_file):
            parameter_info={}
            f = open(parameter_info_file, 'rt')
            for l in f:
                fields = l.strip().split("\t")
                parameter_info[fields[0]]=fields[1]

            self.parameter_set_name = parameter_info["PARAMETER_SET"].strip()
            self.dist_between_bc_for_read = parameter_info["DIST_BETWEEN_BC_FOR_READ"].strip()
            self.R1_regex_looseness = parameter_info["R1_REGEX_LOOSENESS"].strip()
            self.R2_regex_looseness = parameter_info["R2_REGEX_LOOSENESS"].strip()
            self.read_error_size_threshold = parameter_info["READ_ERROR_SIZE_THRESHOLD"].strip()
            self.bc_similarity_threshold = parameter_info["BC_SIMILARITY_THRESHOLD"].strip()
            self.contam_removal_threshold = parameter_info["CONTAMINATION_REMOVAL_THRESHOLD"].strip()
    return(parameter_info())

def write_input_file(root, project_info, parameter_info):
    fname = root + "/tubaseq_inp_files/" + project_info.project_name + "_" + parameter_info.parameter_set_name + ".inp"
    f = open(fname, 'wt')
    for sample in project_info.sample_ids:
        f.write("--sample={} --sample_path={} --root={} --project_name={} --parameter_name={}\n".format(sample, project_info.sample_to_sample_paths[sample], root, project_info.project_name, parameter_info.parameter_set_name))

def make_regexes(R1_regex_looseness, R2_regex_looseness):
    forward_structure = "GA" + "(........)" + "(GC.....TA.....GC.....TA.....GC)" + "(ATGCCCA){e<"+str(R1_regex_looseness) +"}"
    reverse_structure ='(........)' +'(GC.....TA.....GC.....TA.....GC)' + '(ATGCCCA){e<'+ str(R2_regex_looseness) +"}"
    return([forward_structure,reverse_structure])

def make_sgIDDict(sgids, sgRNAs):
    sgIDdict = {}
    for i, sgid in enumerate(sgids):
        sgIDdict[sgid] = sgRNAs[i]
    return(sgIDdict)

def getsgID(sgIDdict, sgid, sgiddist):
    '''flexbile getsgID returns match if sgid is within sgiddist
    from valid sgid'''
    for valid_sgid in sgIDdict.keys():
        if hamming_distance(valid_sgid, sgid)<sgiddist:
            return(sgIDdict[valid_sgid])
    return("None")

def pull_sgid(sgid_bc):
    return(sgid_bc.rsplit("_",1)[0])

def pull_bc(sgid_bc):
    return(sgid_bc.split("_")[-1])

def ListToDict(_list):
    d = {}
    for i in _list:
        fields = i.strip().split(",")
        d[fields[1]]=fields[2]
    return d

def onlyGoodSpikeIns(SpikeInData, expected_spi_bc):
    '''
    Returns true if the barcodes with the highest read counts
    among those with the spike-in sgID are the expected sequences.
    '''
    print("SpikeInData is {}\n".format(SpikeInData))
    if isinstance(SpikeInData, dict):

        for bc in SpikeInData:
            if bc not in expected_spi_bc:
                return False
        return True
    elif isinstance(SpikeInData, list):
        for i in SpikeInData:
            fields = i.strip().split(",")
            bc = fields[1]
            if bc not in expected_spi_bc:
                return False
        return True
    else:
        sys.exit("SpikeInData not correct object type\n")

def get_conversion_factor_unweighted(bc_to_sizes_theoretical, actual_data):
    '''
    Returns RC --> CellNumber conversion factor. aka cells per read
    '''
    if not actual_data: #you hit this if there were spike-ins but none of them were the correct barcode.
        factor = 1

    ratios = [int(bc_to_sizes_theoretical[bc])/int(actual_data[bc]) for bc in actual_data]

    return np.mean(ratios)
