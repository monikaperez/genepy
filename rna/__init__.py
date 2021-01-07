# Jeremie Kalfon
# for BroadInsitute
# in 2019

from __future__ import print_function
import warnings
from matplotlib import pyplot as plt
from bokeh.palettes import *
from bokeh.plotting import *
from JKBio.rna import pyDESeq2
from JKBio.utils import helper as h
import pdb

import pandas as pd
import numpy as np

from taigapy import TaigaClient
tc = TaigaClient()


def filterProteinCoding(listofgenes, from_idtype='ensembl_gene_id'):
    """
    Given a list of genes, provide the args where the genes are protein coding genes:

    This functtion will use a file in taiga, you need taigapy installed

    Args:
    -----
      listofgenes: list of genes
      from_idtype: one of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name", the gene name format

    Returns:
    -------
      the args where the genes are protein coding
    """
    tokeep = []
    b = 0
    print("you need access to taiga for this (https://pypi.org/project/taigapy/)")
    gene_mapping = tc.get(name='hgnc-87ab', file='hgnc_complete_set')
    for i, val in enumerate(listofgenes):
        if from_idtype == "ensembl_gene_id":
            val = val.split(".")[0]
        elif from_idtype == "hgnc_id":
            val = "HGNC:" + str(val)
        a = gene_mapping["locus_group"][gene_mapping[from_idtype]
                                        == val].values
        if len(a) > 0:
            if a[0] == "protein-coding gene":
                tokeep.append(i)
        else:
            b += 1
    print(str(b))
    return(tokeep)


def convertGenes(listofgenes, from_idtype="ensembl_gene_id", to_idtype="symbol"):
    """
    Given a list of genes, provide the args where the genes are protein coding genes:

    This functtion will use a file in taiga, you need taigapy installed

    Args:
    -----
      listofgenes: list of genes
      from_idtype: one of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name", the gene name format
      to_idtype: one of "symbol","uniprot_ids","pubmed_id","ensembl_gene_id","entrez_id","name", the gene name format

    Returns:
    -------
      1: the new names for each genes that were matched else the same name
      2: the names of genes that could not be matched
    """
    print("you need access to taiga for this (https://pypi.org/project/taigapy/)")
    gene_mapping = tc.get(name='hgnc-87ab', file='hgnc_complete_set')
    not_parsed = []
    renamed = []
    b = 0
    to = {}
    for i, val in gene_mapping.iterrows():
        to[val[from_idtype]] = val[to_idtype]
    for i, val in enumerate(listofgenes):
        if from_idtype == "ensembl_gene_id":
            val = val.split(".")[0]
        elif from_idtype == "hgnc_id":
            val = "HGNC:" + str(val)
        try:
            a = to[val]
            renamed.append(int(a) if to_idtype == 'entrez_id' else a)
        except KeyError:
            b += 1
            not_parsed.append(val)
            renamed.append(val)
    print(str(b) + " could not be parsed... we don't have all genes already")
    return(renamed, not_parsed)


def getSpikeInControlScales(refgenome, fastq=None, fastQfolder='', mapper='bwa', pairedEnd=False, cores=1,
                            pathtosam='samtools', pathtotrim_galore='trim_galore', pathtobwa='bwa',
                            totrim=True, tomap=True, tofilter=True, results='res/', toremove=False):
    """
    Will extract the spikeInControls from a fastq file (usefull for, let say ChIPseq data with spike ins)

    Count based sequencing data is not absolute and will be normalized as each sample will be sequenced at a specific depth. 
    To figure out what was the actual sample concentration, we use Spike In control
    You should have FastQfolder/[NAME].fastq & BigWigFolder/[NAME].bw with NAME being the same for the same samples


    Args:
    -----
      refgenome: str the file path to the indexed reference genome
      FastQfolder: str the folder path where the fastq files are stored (should be named the same as files in BigWigFolder)
      BigWigFolder: str the folder path where the bigwig files are stored (should be named the same as files in FastQfolder)
      mapper: str flag to 'bwa', ...
      pairedEnd: Bool flat to true for paired end sequences. if true, You should have FastQfolder/[NAME]_1|2.fastq

    Returns:
    --------
      dict(file,float) the scaling factor dict

    """
    if len(fastQfolder) > 0:
        print('using all files from folder')
        fastqs = os.listdir(fastQfolder)
        fastqs = [i for i in fastqs if '.fq.gz' ==
                  i[-6:] or '.fastq.gz' == i[-9:]]
        fastqs.sort()
        if pairedEnd and (tomap or totrim):
            print("need to be name_*1, name_*2")
            fastqs = [i for i in grouped(fastqs, 2)]
    elif fastq is None:
        raise Error('you need input files')
    else:
        if type(fastq) is list:
            print('your files need to be all in the same folder')
            fastQfolder = '/'.join(fastq[0].split('/')[:-1]) + '/'
            if not totrim and not tomap:
                fastqs = [f.split('/')[-1] for f in fastq]
            else:
                print("need to be name_*1, name_*2")
                fastqs = [[f[0].split('/')[-1], f[1].split('/')[-1]]
                          for f in grouped(fastq, 2)]
        else:
            fastQfolder = '/'.join(fastq.split('/')[:-1]) + '/'
            fastqs = [fastq.split('/')[-1]]
    print(fastqs)
    if not totrim:
        print("you need to have your files in the " + results + " folder")
    if totrim and tomap:
        print("\n\ntrimming\n\n")
        if pairedEnd:
            cmds = []
            rm = []
            for file in fastqs:
                cmd = pathtotrim_galore + ' --paired --fastqc --gzip ' + fastQfolder + \
                    file[0] + ' ' + fastQfolder + file[1] + " -o " + results
                if toremove:
                    rm.append('rm ' + fastQfolder +
                              file[0] + ' ' + fastQfolder + file[1])
                cmds.append(cmd)
            print(cmds)
            parrun(cmds, cores, add=rm)
            fastqs = [[file[0].split('.')[
                0] + '_val_1.fq.gz', file[1].split('.')[0] + '_val_2.fq.gz'] for file in fastqs]
    if tomap:
        print("\n\nmapping\n\n")
        if pairedEnd:
            cmds = []
            rm = []
            for file in fastqs:
                cmd = pathtobwa + ' mem ' + refgenome + ' ' + results + file[0] + ' ' + results +\
                    file[1] + ' | ' + pathtosam + ' sort - -o ' + \
                    results + file[0].split('.')[0] + '.sorted.bam'
                if toremove:
                    rm.append('rm ' + results +
                              file[0] + ' ' + results + file[1])
                cmds.append(cmd)
            parrun(cmds, cores, add=rm)
            fastqs = [file[0].split('.')[0] + '.sorted.bam' for file in fastqs]

    if tofilter:
        print("\n\nfiltering\n\n")
        cmds = []
        rm = []
        parrun([pathtosam + ' index ' + results + file.split('.')
                [0] + '.sorted.bam' for file in fastqs], cores)
        parrun([pathtosam + ' flagstat ' + results + file.split('.')[0] + '.sorted.bam > ' +
                results + file.split('.')[0] + '.sorted.bam.flagstat' for file in fastqs], cores)
        parrun([pathtosam + ' idxstats ' + results + file.split('.')[0] + '.sorted.bam > ' +
                results + file.split('.')[0] + '.sorted.bam.idxstat' for file in fastqs], cores)
        fastqs = [file.split('.')[0] + '.sorted.bam' for file in fastqs]
    else:
        print("files need to be named: NAME.sorted.bam")
        fastqs = [file for file in fastqs if '.sorted.bam' == file[-11:]]
    mapped = {}
    norm = {}
    unique_mapped = {}
    print("\n\ncounting\n\n")
    for file in fastqs:
        mapped[file.split('.')[0]] = int(os.popen(pathtosam + ' view -c -F 0x004 -F 0x0008 -f 0x001 -F 0x0400 -q 1 ' + results +
                                                  file + ' -@ ' + str(cores)).read().split('\n')[0])
       # unique_mapped[file.split('.')[0]] = int(re.findall("Mapped reads: (\d+)", os.popen('bamtools stats -in '+results +
        #                                             file + '.sorted.bam').read())[0])
    nbmapped = np.array([i for i in mapped.values()])
    nbmapped = np.sort(nbmapped)[0] / nbmapped.astype(float)
    for i, val in enumerate(mapped.keys()):
        norm[val] = nbmapped[i]
    return norm, mapped,  # unique_mapped


def GSEAonExperiments(data, experiments, res={}, savename='', scaling=[], geneset='GO_Biological_Process_2015',
                      cores=8, cleanfunc=lambda i: i.split('(GO')[0]):
    """

    Will run GSEA on a set of experiment

    Args:
    -----
      data: a pandas.df rows: gene counts; columns: [experimentA_X,..., experimentD_X..., control_X] where X is the replicate number
      experiments: a list of experiment names (here experimentA,.. experimentD)
      scaling: a dict(experiment:(mean,std)) of scaling factors and their associated standard error for each experiments
      res: you can provide a dict containing results from
      savename: if you want to save the plots as pdfs, provides a location/name
      geneset: the geneset to run it on. (can be a filepath to your own geneset)
      cores: to run GSEA on
      cleanfunc: a func applied to the names of the gene sets to change it in some way (often to make it more readable)
    Returns
    -------
      plots the results 
      1: returns a matrix with the enrichment for each term for each experiment
      2: returns a dict(experiment:pd.df) with dataframe being the output of GSEA (with pvalues etc..) for each experiments
    """
    for i, val in enumerate(experiments):
        print(val)
        totest = data[[v for v in data.columns[:-1]
                       if val + '-' in v or 'AAVS1' in v]]
        cls = ['Condition' if val + '-' in v else 'DMSO' for v in totest.columns]
        if scaling:
            if abs(scaling[val.split('_')[1]][0]) > scaling[val.split('_')[1]][1]:
                print("rescaling this one")
                cols = [i for i in totest.columns if val + '-' in i]
                totest[cols] = totest[cols] * \
                    (2**scaling[val.split('_')[1]][0])
        if val in res:
            print(val + " is already in set")
            continue
        res[val] = gseapy.gsea(data=totest, gene_sets=geneset,
                               cls=cls, no_plot=False, processes=cores)
        res[val].res2d['Term'] = [i for i in res[val].res2d.index]
        for i, v in res.items():
            res[i].res2d['Term'] = [cleanfunc(i) for i in v.res2d['Term']]
        plt.figure(i)
        sns.barplot(data=res[val].res2d.iloc[:25], x="es", y="Term",
                    hue_order="geneset_size").set_title(val)
    a = set()
    for k, val in res.items():
        a.update(set(val.res2d.Term))
    a = {i: [0] * len(res) for i in a}
    for n, (k, val) in enumerate(res.items()):
        for i, v in val.res2d.iterrows():
            a[v.Term][n] = v.es
    pres = pd.DataFrame(a, index=res.keys())
    a = sns.clustermap(figsize=(25, 20), data=res, vmin=-1,
                       vmax=1, yticklabels=res.index, cmap=plt.cm.RdYlBu)
    b = sns.clustermap(-res.T.corr(), cmap=plt.cm.RdYlBu, vmin=-1, vmax=1)
    if savename:
        res.to_csv(savename + ".csv")
        a.savefig(savename + "_genesets.pdf")
        b.savefig(savename + "_correlation.pdf")
    return pres, res


def runERCC(ERCC, experiments, featurename="Feature", issingle=False, dilution=1 / 100,
            name="RNPv2", spikevol=1, control="AAVS1", fdr=0.1, totalrnamass=0.5):
    """
    Runs the ERCC dashboard Rpackage on your notebook

    you will need to run this function from ipython and to have the R package erccdashboard installed

    Args:
    ----
      ERCC: a pandas.df rows: ERCC counts columns: [experimentA_X,..., experimentD_X..., control_X] where X is the replicate number
      experiments: a list of experiment names (here experimentA,.. experimentD)
      featurename: columns where the ERCC pseudo gene names are stored
      issingle: ERCC parameters to choose between Single and RatioPair
      dilution: ERCC dilution parameter
      name: the name of the experiment set
      spikevol: ERCC spikevol parameter
      control: the control name (here control)
      fdr: ERCC fdr parameter
      totalrnamass: ERCC totalrnamass parameter

    Returns:
    -------
      a dict(experimentName:(val, ste)) a dict containing the scaling factor and its standard error for each experiment

    Raises:
    ------
      RuntimeError: if you are not on ipython
    """
    try:
        ipython = get_ipython()
    except:
        raise RuntimeError('you need to be on ipython')
    ipython.magic("load_ext rpy2.ipython")
    ipython.magic("R library('erccdashboard')")
    # "count" for RNA-Seq data, "array" for microarray data
    ipython.magic("R datType = 'count'")
    # flag to indicate if input expression measures are already normalized, default is FALSE
    ipython.magic("R isNorm = F")
    # user defined filename prefix for results files
    ipython.magic("R -i name filenameRoot = name")
    # name for sample 2 in the experiment
    ipython.magic("R -i control sample2Name = control")
    # name of ERCC mixture design, "RatioPair" is default
    ipython.magic(
        "R -i issingle erccmix <- if(issingle) 'Single' else 'RatioPair'")
    # dilution factor used for Ambion spike-in mixtures
    ipython.magic("R -i dilution erccdilution = dilution")
    # volume (in microliters) of diluted spike-in mixture added to total RNA mass
    ipython.magic("R -i spikevol spikeVol = spikevol")
    # user defined false discovery rate (FDR), default is 0.05
    ipython.magic("R -i fdr choseFDR = fdr")
    ipython.magic("R exDat = ''")
    ipython.magic("R -i totalrnamass totalRNAmass <- totalrnamass")

    cols = list(ERCC.columns)
    cols.sort()
    res = {}
    for val in experiments:
        d = {}
        e = 0
        c = 0
        d.update({
            featurename: 'Feature'
        })
        for i in cols:
            if val + '-' in i:
                e += 1
                d.update({i: val.split('_')[-1] + '_' + str(e)})
            if control + "-" in i:
                c += 1
                d.update({i: control + "_" + str(c)})
        a = ERCC[list(d.keys())].rename(columns=d)
        a.to_csv('/tmp/ERCC_estimation.csv', index=None)
        val = val.split('_')[-1]

        ipython.magic("R -i val print(val)")

        ipython.magic("R print(sample2Name)")
        ipython.magic("R a <- read.csv('/tmp/ERCC_estimation.csv')")
        ipython.magic("R print(head(a))")

        try:

            ipython.magic("R -i val exDat = initDat(datType=datType, isNorm=isNorm, exTable=a,\
                                 filenameRoot=filenameRoot, sample1Name=val,\
                                 sample2Name=sample2Name, erccmix=erccmix,\
                                 erccdilution=erccdilution, spikeVol=spikeVol,\
                                 totalRNAmass=totalRNAmass, choseFDR=choseFDR)")
            ipython.magic("R exDat = est_r_m(exDat)")
            ipython.magic("R exDat = dynRangePlot(exDat)")

        except Warning:
            print("failed for " + val)
            continue
        except:
            print('worked for ' + val)

        ipython.magic("R print(summary(exDat))")
        ipython.magic("R grid.arrange(exDat$Figures$dynRangePlot)")
        ipython.magic("R grid.arrange(exDat$Figures$r_mPlot)")
        ipython.magic("R grid.arrange(exDat$Figures$rangeResidPlot)")

        ipython.magic("R -o rm rm <- exDat$Results$r_m.res$r_m.mn")
        ipython.magic("R -o se se <- exDat$Results$r_m.res$r_m.mnse")
        ipython.magic("R write.csv(c(rm,se), file = '/tmp/JKBIO_ercc.csv')")
        ipython.magic("R print(se,rm)")
        l = h.fileToList("/tmp/JKBIO_ercc.csv")
        res[val] = (float(l[1][4:]), float(l[2][4:]))
    for i, v in res.items():
        if abs(v[0]) > v[1]:
            print(i, v[0])
    return res


def mergeSplicingVariants(df, defined='.'):
    df = df.sort_index()
    foundpoint = False
    # pdb.set_trace()
    torename = {}
    todrop = []
    for i, v in enumerate(df.index.tolist()):
        showcount(i, len(df))
        if foundpoint:
            if foundpoint in v:
                tomerge.append(v)
            else:
                if foundpoint not in df.index:
                    if len(tomerge) > 1:
                        #print("merging "+str(len(tomerge)))
                        df.loc[foundpoint] = df.loc[tomerge].sum()
                        todrop.extend(tomerge)
                    else:
                        torename.update({tomerge[0]: foundpoint})
                else:
                    todrop.extend(tomerge)
                    tomerge.append(foundpoint)
                    df.loc[foundpoint] = df.loc[tomerge].sum()
                foundpoint = False
        elif defined in v:
            foundpoint = v.split(defined)[0]
            tomerge = [v]
    if len(torename) > 0:
        df = df.rename(index=torename)
    df = df.drop(index=todrop)
    return df

def readFromSlamdunk(loc='res/count/', flag_var=100, convertTo='symbol',
                     minvar_toremove=0, mincount_toremove=5):
    """
    
    Args:
    -----

    Returns:
    --------

    """
    files = os.listdir(loc)
    files = [file for file in files if file.endswith('.tsv')]
    data = {}
    for file in files:
        data[file.split('/')[-1].split('.')[0]] = pd.read_csv(loc +
                                                              file, sep='\t', comment='#', header=0)
    prev = -2
    print("found " + str(len(data)) + ' files:' + str(data.keys()))
    for k, val in data.items():
        if len(set(val.Name)) != prev and prev != -2:
            raise ValueError(
                'we do not have the same number of genes in each file')
        prev = len(set(val.Name))
    readcounts = {i: [0] * len(data) for i in val.Name.unique()}
    tccounts = {i: [0] * len(data) for i in val.Name.unique()}
    for n, (_, val) in enumerate(data.items()):
        val = val.sort_values(by="Name")
        j = 0
        #print('              ',end='\r')
        readcount = [val.iloc[0].ReadCount]
        tccount = [val.iloc[0].TcReadCount]
        prevname = val.iloc[0].Name
        for _, v in val.iloc[1:].iterrows():
            if v.Name == 4609:
                print(v.ReadCount, v.TcReadCount)
                print(readcount, tccount)
            if v.Name == prevname:
                readcount.append(v.ReadCount)
                tccount.append(v.TcReadCount)
            else:
                readcounts[prevname][n] = np.sum(readcount)
                tccounts[prevname][n] = np.sum(tccount)
                # if np.var(readcount) > flag_var:
                #    print("pb with "+str(v.Name))
                prevname = v.Name
                j += 1
                # print(j,end='\r')
                readcount = [v.ReadCount]
                tccount = [v.TcReadCount]
    files = [*data]
    readcounts = pd.DataFrame(
        data=readcounts, columns=val.Name.unique(), index=data.keys()).T
    tccounts = pd.DataFrame(
        data=tccounts, columns=val.Name.unique(), index=data.keys()).T
    if convertTo:
        names, _ = convertGenes(readcounts.index.tolist(
        ), from_idtype="entrez_id", to_idtype=convertTo)
        readcounts.index = names
        names, _ = convertGenes(tccounts.index.tolist(
        ), from_idtype="entrez_id", to_idtype=convertTo)
        tccounts.index = names
    nottodrop = np.argwhere(tccounts.values.var(1) >= minvar_toremove).ravel()
    tccounts = tccounts.iloc[nottodrop]
    readcounts = readcounts.iloc[nottodrop]
    nottodrop = np.argwhere(readcounts.values.max(1) >=
                            mincount_toremove).ravel()
    tccounts = tccounts.iloc[nottodrop]
    readcounts = readcounts.iloc[nottodrop]
    return readcounts, tccounts


def DESeqSamples(data, experiments, scaling=None, keep=True, rescaling=None, results={}, controlcontain='RNP_AAVS1',
                 spikecontrolscontain="ERCC-", threshforscaling=2):
    """
    Args:
    ----

    Returns:
    ------
    
    """
    if "gene_id" not in list(data.columns):
        print("using index as gene_id")
        data['gene_id'] = data.index
    warnings.simplefilter("ignore")
    if type(controlcontain) is str:
        controlcontain = [controlcontain]*len(experiments)
    for j, val in enumerate(experiments):
        print(val)
        cond = [1 if val+'-' in i else 0 for i in data.columns[:-1]]
        contr = [1 if controlcontain[j] in i else 0 for i in data.columns[:-1]]
        design = pd.DataFrame(index=data.columns[:-1], columns=['DMSO', 'Target'],
                              data=np.array([contr, cond]).T)
        design.index = design.index.astype(str).str.replace('-', '.')
        deseq = pyDESeq2.pyDESeq2(count_matrix=data, design_matrix=design,
                                  design_formula='~DMSO + Target', gene_column="gene_id")
        if type(scaling) is bool:
            print("scaling using ERCC")
            if scaling:
                deseq.run_estimate_size_factors(
                    controlGenes=data.gene_id.str.contains(spikecontrolscontain))
        elif type(scaling) is list or type(scaling) is set:
            print("scaling using a gene set")
            deseq.run_estimate_size_factors(controlGenes=data.gene_id.isin(scaling))
        elif type(scaling) is dict:
            if val in scaling:
                print("auto scaling from ERCCdashboard mean/std values")
                if abs(scaling[val][0]) > threshforscaling*scaling[val][1]:
                    print("  estimating sizeFactors for this one")
                    deseq.run_estimate_size_factors(
                        controlGenes=data.gene_id.str.contains(spikecontrolscontain))
                    if rescaling is not None:
                        if val in rescaling:
                            sizeFact = deseq.getSizeFactors()
                            sizeFact[np.where(cond)[0]
                                     ] *= rescaling[val.split('_')[1]]
                            deseq.setSizeFactors(sizeFact)
            else:
                print("  not in scaling dict")
        elif val in results and keep:
            continue
        deseq.run_deseq()
        deseq.get_deseq_result()
        r = deseq.deseq_result
        r.pvalue = np.nan_to_num(np.array(r.pvalue), 1)
        r.log2FoldChange = np.nan_to_num(np.array(r.log2FoldChange), 0)
        results[val] = r
        print('____________________________________________')
    return results


def gsva(data, geneset_file, pathtoJKBio, method='ssgsea'):
  data.to_csv('/tmp/data_JKBIOhelper_gsva.csv', index_col=0)
  cmd = "Rscript "+pathtoJKBio + \
      "/helper/ssGSEA.R /tmp/data_JKBIOhelper_gsva.csv " + geneset_file + " " + ssgsea
  res = subprocess.run(cmd, shell=True, capture_output=True)
  if res.returncode != 0:
    raise ValueError('issue with the command: ' + str(res.stderr))
  print(res)
  res = pd.read_csv("/tmp/res_JKBIO_ssGSEA.csv")
  return res
