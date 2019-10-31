import os
import pandas as pd
import pysam
import numpy as np
import Helper
import re
from pybedtools import BedTool
import seaborn as sns
import pyBigWig
import matplotlib.pyplot as plt
import pdb
from scipy.optimize import curve_fit
from sklearn.preprocessing import normalize
# from scipy.misc import factorial
from scipy.stats import poisson


size = {"GRCh37": 2864785220,
        "GRCh38": 2913022398}

cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
         'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
         'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']


def extractPairedSingleEndFrom(folder, pattern='R1/R2', sep='-', namepos=2):
    """
    given a folder, find using a specific patterns, single end and paired end files
    return a list of single end files and a list of list of paired end files [[R1,R2]]
    """
    single = []
    paired = {}
    for val in os.listdir(folder):
        if 'R1' in val:
            name = val.split(sep)[namepos]
            paired[name] = {'R1': val}
        elif 'R2' in val:
            name = val.split(sep)[namepos]
            paired[name].update({'R2': val})
        else:
            single.append(val)

    return single, pd.DataFrame(paired)


def findReplicates(folder, sep='-', namings='-r([0-9])', namepos=2):
    """
    creates a dict of name and replicate files
    """
    rep = {}
    for val in os.listdir(folder):
        if val[-4:] == '.bam':
            match = re.search(namings, val)
            if match:
                number = match.groups()[0]
                name = val.split(sep)[namepos]
                if name in rep:
                    rep[name].append(val)
                else:
                    rep[name] = [val]

    return rep


def computeSingleEnd(singlend, folder="data/seqs/", numthreads=8, peaksFolder="peaks/",
                     ismapped=False, mappedFolder='mapped/', refFolder='data/reference/index'):
    """
    run the singleEnd pipeline
    for alignment etc, one can use pysam ready made implementation of samtools
    """
    for val in singlend:
        out1 = folder + mappedFolder + val.split('.')[0] + ".mapped.sam"
        if not ismapped:
            in1 = folder + val
            os.system("bowtie2 -x " + refFolder + " --threads " + str(numthreads) +
                      " -t -k 1 --very-sensitive -U " + in1 + " -S " + out1)
        out2 = folder + peaksFolder + val.split('.')[0]
        print(out1)
        os.system("macs2 callpeak -f SAM -t " + out1 + " --outdir " + out2)
        # it can take many TB so better delete


def computePairedEnd(pairedend, folder="data/seqs/", numthreads=8, peaksFolder="peaks/",
                     ismapped=False, mappedFolder='mapped/', refFolder='data/reference/index'):
    """
    # run the paired end pipeline
    """
    for key, val in pairedend.items():
        out1 = folder + mappedFolder + val[0].split('.')[0] + ".mapped.sam"
        in1 = folder + val[0]
        in2 = folder + val[1]
        os.system("bowtie2 -x " + refFolder + " --threads " + str(numthreads) + " -t -k 1 \
      --very-sensitive -1 " + in1 + " -2 " + in2 + " - S " + out1)
        out2 = folder + peaksFolder + val[0].split('.')[0]
        print(out1)
        changefrom = out1
        changeto = out1[:-4] + '.bam'
        os.system("samtools view -b " + changefrom + " -o " + changeto)
        os.system("macs2 callpeak --format 'BAMPE' --treatment " + changeto + " --outdir " + out2)
        # it can take many TB so better delete


def bigWigFrom(bams, folder="data/seqs/", numthreads=8, genome='GRCh37'):
    """
    run the bigwig command line for a set of bam files in a folder
    """
    for i in bams:
        in1 = folder + i
        out1 = folder + "bigwig/" + i + split('.')[0] + '.bw'
        os.system("bamCoverage --effectiveGenomeSize " + size[genome] + " -p " + str(numthreads) +
                  " -b " + in1 + "-of bigwig -o " + out1)


def mergeBams(rep):
    """
    uses samtools to merge a set of replicates considered into one file
    """
    for i, val in rep.items():
        out1 = i + '.merged.bam'
        for bam in val:
            in1 += ' ' + bam
        os.system("samtools merge " + out1 + in1)


def loadNarrowPeaks(peakfolder="data/peaks/", isMacs=True, CTFlist=[], skiprows=0):
    # for each data peak type file listed in a MACS2 way in a given MACS2 output folder, will merge them
    # all into one dataframe and output the dataframe
    bindings = pd.DataFrame()
    for folder in os.listdir(peakfolder):
        if isMacs:
            if any(tf in folder for tf in CTFlist) or not CTFlist:
                binding = pd.read_csv(peakfolder + folder + "/NA_peaks.narrowPeak", sep='\t', header=None)
                binding['name'] = folder.replace('.narrowPeak', '')
                bindings = bindings.append(binding)
        else:
            if folder[-11:] == ".narrowPeak" and (any(tf in folder for tf in CTFlist) or not CTFlist):
                binding = pd.read_csv(peakfolder + folder, sep='\t', header=None, skiprows=skiprows)
                binding['name'] = folder.replace('.narrowPeak', '')
                bindings = bindings.append(binding)
    bindings = bindings.drop(5, 1).drop(4, 1)
    bindings = bindings.rename(columns={
        0: "chrom",
        1: 'start',
        2: 'end',
        3: 'peak_number',
        6: "foldchange",
        7: "-log10pvalue",
        8: "-log10qvalue",
        9: 'relative_summit_pos'})
    bindings = bindings.sort_values(by=["chrom", "start", "end"], axis=0)
    bindings.start = bindings.start.astype('int')
    bindings.end = bindings.end.astype('int')
    bindings.size = bindings.size.astype('int')
    bindings.relative_summit_pos = bindings.relative_summit_pos.astype('int')
    bindings.foldchange = bindings.foldchange.astype('float')
    bindings["-log10pvalue"] = bindings["-log10pvalue"].astype('float')
    bindings['-log10qvalue'] = bindings['-log10qvalue'].astype('float')
    return bindings


def Pysam_computePeaksAt(peaks, bams, folder='data/seqs/', window=1000, numpeaks=1000, numthreads=8):

    # get pysam data
    # ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
    # for each counts, do a rolling average (or a convolving of the data) with numpy
    # append to an array
    # return array, normalized
    loaded = {}
    res = {i: np.zeros((len(peaks), window * 2)) for i in bams}
    peaks = peaks.sort_values(by="foldchange", ascending=False).iloc[:numpeaks]
    peaks.chrom = peaks.chrom.astype(str)
    for val in bams:
        loaded.update({val: pysam.AlignmentFile(folder + val, 'rb', threads=numthreads)})
    for k, bam in loaded.items():
        for num, (i, val) in enumerate(peaks.iterrows()):
            print(num / len(peaks), end='\r')
            center = int((val['start'] + val['end']) / 2)
            for pileupcolumn in bam.pileup(val['chrom'], start=center - window,
                                           stop=center + window, truncate=True):
                res[k][num][pileupcolumn.pos - (center - window)] = pileupcolumn.n
    fig, ax = plt.subplots(1, len(res))
    for i, (k, val) in enumerate(res.items()):
        sns.heatmap(val, ax=ax[i])
        ax[i].set_title(k.split('.')[0])
    fig.show()
    return res, fig


def Bedtools_computePeaksAt(peaks, bams, folder='data/seqs/', window=1000, numpeaks=1000, numthreads=8):
    """
    get pysam data
    ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
    for each counts, do a rolling average (or a convolving of the data) with numpy
    append to an array
    return array, normalized
    """
    loaded = {}
    center = [int((val['start'] + val['end']) / 2) for k, val in peaks.iterrows()]
    peaks['start'] = [c - window for c in center]
    peaks['end'] = [c + window - 1 for c in center]
    peaks[peaks.columns[:3]].sort_values(by=['chrom', 'start']).to_csv('temp/peaks.bed', sep='\t', index=False, header=False)
    bedpeaks = BedTool('temp/peaks.bed')

    fig, ax = plt.subplots(1, len(bams))
    peakset = peaks["foldchange"].values.argsort()[::-1][:numpeaks]
    for i, val in enumerate(bams):
        coverage = BedTool(folder + val).intersect(bedpeaks).genome_coverage(bga=True, split=True)\
            .intersect(bedpeaks).to_dataframe(names=['chrom', 'start', 'end', 'coverage'])
        cov = np.zeros((len(peaks), window * 2), dtype=int)
        j = 0
        pdb.set_trace()
        for i, (k, val) in enumerate(peaks.iterrows()):
            print(i / len(peaks), end='\r')
            while coverage.iloc[j].start > val.start:
                j -= 1
            while coverage.iloc[j].start < val.end:
                cov[i][coverage.iloc[j].start - val.start:coverage.iloc[j].end - val.start] =\
                    coverage.iloc[j].coverage
                j += 1
        sns.heatmap(coverage, ax=ax[i])
        ax[i].set_title(val.split('.')[0])
    fig.show()
    return None, fig


def computePeaksAt(peaks, bigwigs, folder='data/seqs/', window=1000, numpeaks=4000, numthreads=8,
                   width=5, length=10, name='temp/peaksat.png'):
    """
    get pysam data
    ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
    for each counts, do a rolling average (or a convolving of the data) with numpy
    append to an array
    return array, normalized
    """
    center = [int((val['start'] + val['relative_summit_pos'])) for k, val in peaks.iterrows()]
    pd.set_option('mode.chained_assignment', None)
    peaks['start'] = [c - window for c in center]
    peaks['end'] = [c + window for c in center]
    fig, ax = plt.subplots(1, len(bigwigs), figsize=[width, length])
    peaks = peaks.sort_values(by=["foldchange"], ascending=False)
    if numpeaks > len(peaks):
        numpeaks = len(peaks) - 1
    for num, bigwig in enumerate(bigwigs):
        bw = pyBigWig.open(folder + bigwig)
        cov = np.zeros((numpeaks, window * 2), dtype=int)
        for i, (k, val) in enumerate(peaks.iloc[:numpeaks].iterrows()):
            cov[i] = np.nan_to_num(bw.values(str(val.chrom), val.start, val.end), 0)
        sns.heatmap(cov, ax=ax[num], yticklabels=[], cmap=cmaps[num],
                    cbar=False)
        ax[num].set_title(bigwig.split('.')[0])
    fig.subplots_adjust(wspace=0.1)
    fig.show()
    fig.savefig(name)
    return cov, fig


"""
def computeMeanCov(bigwigFolder, meanOnly=True, ref="GRCh38", averageFragSize=150, outputfolder='/data/coverage/'):
    meancov = {}
    for val in os.listdir(bigwigFolder):
        bw = pyBigWig.open(folder + bigwig)
        if bw:
            meancov[val.split('.')[0]] = bw.
    return meancov


def substractPeaks(peaks1, to):
    peaks1 = peaks1.sort_values(by=['chrom', 'start'])
    peaks2 = to.sort_values(by=['chrom', 'start'])
    j = 0
    i = 0
    newpeaks = pd.DataFrame(columns=peaks2.columns)
    while i < len(peaks2):
        if peaks1.iloc[j].chrom == peaks2.iloc[j].chrom:
            else
        intersection = helper.intersection([peaks2.iloc[i]['start'], peaks2.iloc[i]['end']],
                                           [peaks1.iloc[j]['start'], peaks1.iloc[j]['start']])
        if intersection:

            j += 1
        while peaks2.iloc[i]['end'] > peaks1.iloc[j]['start']:
            if peaks2.iloc[i]['end'] > peaks1.iloc[j]['end']:
                newpeaks.append(peaks2.iloc[i], ignore_index=True)
            else:
                newpeaks.append({'start': peaks2.iloc[i]['start'],
                                 'end': peaks1.iloc[j]['start'],
                                 'chromomsome': peaks1.iloc[j]['chrom']}, ignore_index=True)
"""


def simpleMergedPeaks(peaks, window=0, totpeaknumber=0):
    peaks = peaks.sort_values(by=['chrom', 'start'])
    tfs = list(set(peaks['name'].tolist()))
    mergedpeaksdict = {}
    remove = []
    peaknumber = 0
    merged_bed = {
        "chrom": [peaks.iloc[0]['chrom']],
        "start": [peaks.iloc[0]['start']],
        "end": [],
        "peak_number": [peaknumber + totpeaknumber],
        "foldchange": [],
        "-log10pvalue": [],
        "-log10qvalue": [],
        "relative_summit_pos": []
    }
    foldchange = [peaks.iloc[0].get('foldchange', 1)]
    log10pvalue = [peaks.iloc[0].get('-log10pvalue', 0)]
    log10qvalue = [peaks.iloc[0].get('-log10qvalue', 0)]
    relative_summit_pos = peaks.iloc[1].get('relative_summit_pos', peaks.iloc[0]['start'])
    # computes overlap by extending a bit the window (100bp?) should be ~readsize
    prev_end = peaks.iloc[0]['end']
    prev_chrom = peaks.iloc[0]['chrom']
    tfmerged = {a: [0] for a in tfs}
    tfmerged[peaks.iloc[0]['name']][-1] = 1
    for i, (pos, peak) in enumerate(peaks.iloc[1:].iterrows()):
        print(str(i / len(peaks)), end="\r")
        if prev_end + window > peak['start'] and prev_chrom == peak['chrom']:
            # can be merged
            foldchange.append(peak.get('foldchange', 1))
            log10pvalue.append(peak.get('-log10pvalue', 0))
            log10qvalue.append(peak.get('-log10qvalue', 0))
            if peak.get('foldchange', 1) > max(foldchange):
                relative_summit_pos = peak.get('relative_summit_pos', peaks['start'])
        else:
            # newpeak
            for k, val in tfmerged.items():
                val.append(0)
            peaknumber += 1
            merged_bed['chrom'].append(peak['chrom'])
            merged_bed['start'].append(peak['start'])
            merged_bed['end'].append(prev_end)
            merged_bed['peak_number'].append(peaknumber + totpeaknumber)
            merged_bed['foldchange'].append(np.mean(foldchange))
            merged_bed['-log10pvalue'].append(max(log10pvalue))
            merged_bed['-log10qvalue'].append(max(log10qvalue))
            merged_bed['relative_summit_pos'].append(relative_summit_pos)
            foldchange = [peak.get('foldchange', 1)]
            log10pvalue = [peak.get('-log10pvalue', 0)]
            log10qvalue = [peak.get('-log10qvalue', 0)]
            relative_summit_pos = peak.get('relative_summit_pos', peak['start'])
        prev_end = peak['end']
        prev_chrom = peak['chrom']
        tfmerged[peak['name']][-1] = 1
    merged_bed['end'].append(prev_end)
    merged_bed['foldchange'].append(np.mean(foldchange))
    merged_bed['-log10pvalue'].append(max(log10pvalue))
    merged_bed['-log10qvalue'].append(max(log10qvalue))
    merged_bed['relative_summit_pos'].append(relative_summit_pos)

    merged_bed = pd.DataFrame(merged_bed)
    tfmerged = pd.DataFrame(tfmerged)
    return pd.concat([merged_bed, tfmerged], axis=1, sort=False)


def findpeakpath(folder, peakname):
    for val in os.listdir(folder):
        if peakname in val:
            return val
    raise ValueError('no peaks')


def mergeReplicatePeaks(peaks, reps, bigwigfolder, markedasbad=None, window=200,
                        sampling=10000, mincov=4, doPlot=True, cov={}, minKL=2, use='max',
                        MINOVERLAP=0.5, SIZECUTOFF=0.7, mergewindow=50):
    """
    /!/ should only be passed peaks with at least one good replicate
    for each TFpeaksets,
    1. find the replicate that have the most peaks
    2. correlate peaks and get in highest correlation order with the replicate found in 1
    3. find overlap of both and get size of second replicate
    4. if small(er)-> use only to increase statistics
      1. if a lot of uncalled peaks in replicate 2 at replicate 1 peaks (flag for mergebam)
    5. if similar size -> get only intersect
      2. add to intersect, find uncalled peaks in both replicates which are called in the other
    6. repeat for all replicates
    -------------------------
    if full overlap of one of the peak replicate, only use the overlapped one to increase confidence on peak
    if >80% average non overlap,
      print warning and percentage of overlap

    if <20% average non overlap,
      take the overlap and increase confidence and avg logfold

    if one is <20%:
      if other <40% average non overlap,
        take the overlap and increase confidence and avg logfold
      else
        take

    gets the max cov at the genomic window and if above some threshold, accepts the peak.

    extend peak by X bp if no TSS
    remove TSS from peaks

    create a new data frame containing merged peak size, reassembled peak data (p value etc..) and
    a boolean value for presence of each TF listed in previous df
    ------------------------------------

    args:
    ----
    peaks: df[bed-like] all the peaks into the sameBam with a column containing the 'name'
    being the filename of the tf replicate bam

    reps: dict(TFname:[file]) a dictionarry giving a list of bam filename for all replicates of the TF
    bamfolder: str, foldername
    avgCov: dict(filename:int) a dict where for each bam filename is given an averageCoverage
    window:


    returns:
    -------
    mergedpeaks: dict{df-peakslike}
    bamtomerge: [[bam1,bam2]]

    """
    print("/!/ should only be passed peaks with at least one good replicate")
    # for a df containing a set of peaks in bed format and an additional column of different TF
    tfs = list(set(peaks['name'].tolist()))
    totpeaknumber = 0
    mergedpeaksdict = {}
    remove = []
    tomergebam = []
    for tf, rep in reps.items():
        if len(rep) == 1:
            print("we only have one replicate for " + tf + " .. pass")
            mergedpeaksdict.update({tf: peaks[peaks['tf'] == tf]})
            continue
        print("merging " + str(len(rep)) + " " + tf + " peaks")
        merged = simpleMergedPeaks(peaks[peaks['tf'] == tf], window=mergewindow)
        merged_bed = merged[merged.columns[8:]]
        finalpeaks = merged[merged.columns[:8]]
        print('finish first overlaps lookup')
        # flag when  biggest is <1000 peaks
        if len(finalpeaks) < 1000:
            print('!TF has less than 1000 PEAKS!')
        # for each TF (replicates), compute number of peaks
        peakmatrix = merged_bed.values

        # compute overlap matrix (venn?)
        if peakmatrix.shape[1] < 7 and doPlot:
            presence = []
            for peakpres in peakmatrix.T:  # https://github.com/tctianchi/pyvenn
                presence.append(set([i for i, val in enumerate(peakpres) if val == 1]))
            Helper.venn(presence, merged_bed.columns)
        else:
            print('too many replicates for Venn')
        bigwigs = os.listdir(bigwigfolder)
        totpeak = np.sum(peakmatrix, 0)
        sort = np.argsort(totpeak)[::-1]
        biggest_ind = sort[0]
        peakmatrix = peakmatrix.T
        biggest = merged_bed.columns[biggest_ind]
        # starts with highest similarity and go descending
        for i, val in enumerate(sort[1:]):
            # if avg non overlap > 60%, and first, and none small flag TF as unreliable.
            overlap = len(presence[biggest_ind] & presence[val]) / len(presence[biggest_ind])
            peakname = merged_bed.columns[val]
            print('overlap: ' + str(overlap))
            if overlap < MINOVERLAP:
                smallsupport = len(presence[biggest_ind] & presence[val]) / len(presence[val])
                print('not enough overlap')
                if totpeak[val] > totpeak[biggest_ind] * SIZECUTOFF:
                    if i == 0:
                        print("Wrong TF")
                        remove.append(tf)
                        break
                    # if not first, throw the other replicate and continue
                    peakmatrix = np.delete(peakmatrix, val)
                # if small and small overlap more than 80% do findAdditionalPeaks only
                # on for small else throw small

                elif smallsupport > MINOVERLAP:
                    tolookfor = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val] == 0)
                    additionalpeaks = findAdditionalPeaks(finalpeaks, tolookfor,
                                                          bigwigfolder + findpeakpath(bigwigfolder, peakname),
                                                          sampling=sampling, mincov=mincov,
                                                          window=window, minKL=minKL, use=use)
                    peakmatrix[val] = np.logical_or(peakmatrix[val], additionalpeaks)
                    # if new set of peaks >MINOVERLAP size of big, flag for merge bams
                    recovered = np.sum(np.logical_and(peakmatrix[biggest_ind], peakmatrix[val])) / np.sum(
                        peakmatrix[biggest_ind])
                    print("we have recovered " + str(recovered))
                    if recovered > MINOVERLAP:
                        tomergebam.append([biggest, peakname])
                else:
                    print('not enough support: ' + str(smallsupport))
            else:
                print('enough overlap')
                tolookfor = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val] == 0)
                additionalpeaks = findAdditionalPeaks(finalpeaks, tolookfor,
                                                      bigwigfolder + findpeakpath(bigwigfolder, peakname),
                                                      sampling=sampling, mincov=mincov,
                                                      window=window, minKL=minKL, use=use)
                peakmatrix[val] = np.logical_or(peakmatrix[val], additionalpeaks)
                tolookfor = np.logical_and(peakmatrix[val], peakmatrix[biggest_ind] == 0)
                additionalpeaks = findAdditionalPeaks(finalpeaks, tolookfor,
                                                      bigwigfolder + findpeakpath(bigwigfolder, biggest),
                                                      sampling=sampling, mincov=mincov,
                                                      window=window, minKL=minKL, use=use)
                peakmatrix[biggest_ind] = np.logical_or(peakmatrix[biggest_ind], additionalpeaks)
                recovered = np.sum(np.logical_and(peakmatrix[biggest_ind], peakmatrix[val])) / np.sum(
                    peakmatrix[biggest_ind])
                print('we have recovered ' + str(recovered))
                tomergebam.append([biggest, peakname])
                # take the intersection of both peaksets
                peakmatrix[biggest_ind] = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val])
        if tf not in remove:
            finalpeaks = finalpeaks[peakmatrix[biggest_ind].astype(bool)]
            finalpeaks['name'] = tf
            mergedpeaksdict.update({tf: finalpeaks})
    return mergedpeaksdict, tomergebam, remove


def findAdditionalPeaks(peaks, tolookfor, filepath, sampling=10000, mincov=4,
                        window=200, cov={}, minKL=2, use='max'):
    """
    findAdditionalPeaks: for all peaks in A and/or B find in coverage file if zone has relative cov
    of more than thresh then add to peak
    if B is small and > 20% of peaks in A are found back, increase confidence and
    flag for mergeBams
    if < 20% don't flag for merge bam
    f B is big and now mean non overlap < 40%, take union and flag for mergeBam else, throw B.

    returns:
    -------
    peaks: np.array(bool) for each peaks in peakset, returns a binary
    """
    # def poisson(k, lamb, scale): return scale * (lamb**k / factorial(k)) * np.exp(-lamb)

    def KLpoisson(lamb1, lamb2): return lamb1 * np.log(lamb1 / lamb2) + lamb2 - lamb1
    bw = pyBigWig.open(filepath)
    res = np.zeros(len(peaks))
    prevchrom = ''
    lamb = {}
    for i, has in enumerate(tolookfor):
        if has:
            val = peaks.iloc[i]
            if val.chrom != prevchrom:
                if val.chrom not in cov:
                    cov[val.chrom] = bw.stats(str(val.chrom))[0]
                    prevchrom = val.chrom
                if use == 'poisson':
                    samples = np.zeroes(10 * sampling)
                    sampling = np.rand([sampling])
                    sampling = sampling * bw.chroms(str(val.chrom))
                    for j, sample in enumerate(sampling):
                        samples[j:j + 10] = np.nan_to_num(bw.values(str(val.chrom), sample, sample + 10), 0)
                    samples = np.bincount(samples)
                    # curve_fit(poisson,range(len(bins)))
                    lamb[val.chrom] = poisson.fit(samples, range(len(samples) - 1))

            start = max([val.start - window, 0])
            end = min(val.end + window, bw.chroms(str(val.chrom)))
            zone = np.nan_to_num(bw.values(str(val.chrom), start, end), 0)
            if use == 'max':
                if max(zone) / cov[val.chrom] > mincov or sum(zone) / (cov[val.chrom] * (end - start)) > mincov:
                    res[i] = 1
            elif use == 'poisson':
                lamb = poisson.fit(zone,)
                if KLpoisson(lamb, lamb[val.chrom]) > minKL:
                    res[i] = 1

    return res
    # BETTER WAY
    #
    # for each nucleotide in the sequence not in registered peaks,
    # create a count distribution accross nucleotide use bedgraph file output from
    # "bedtools genomecov -d -ibam youralignment.bam" command
    # number of times you see 0,1,2,3,4,5.. reads for each nucleotide. fit an expo distribution to it
    # if above a p-value of .1, mark as found, return the new marks with their p-value


def createCorrMatrix(peaks, bigwigs, window=100):
    """
    somewhat similar concept to computePeaksAt

    # get pysam data
    # ask for counts only at specific locus based on peak windows from mergedpeakset
    # append to an array
    # return array, normalized
    """
    res = np.zeros((len(bigwigs), len(peaks)), dtype=float)
    for i, bw in enumerate(bigwigs):
        print('doing file ' + str(bw))
        bw = pyBigWig.open(bw)
        for k, val in peaks.iterrows():
            start = max([val.start - window, 0])
            end = min(val.end + window, bw.chroms(str(val.chrom)))
            res[i][k] = bw.stats(str(val.chrom), start, end)[0]
    res = np.nan_to_num(res, 0)
    return (res.T / res.max(1)).T


def annotatePeaks():
    """
    get H3k27AC peaks and compute Rose
    for each TF peak
    assign super enhancer if within it and create a super enhancer TFgroup
    for each peaks
    apply similar to rose and merge TF together. (create new TFgroup)
    for each TF groups say
      if within super enhancer
      if within high h3k27ac marks
      its nearest gene (distance to it, cis/trans)


    FIND A WAY TO FILTER TO ONLY PLACES WITH H3K27AC marks??

    TAD points where most contacts on one side happened on this side.
    specific distance. zone of most concentration of contacts for a peak region.
    """
    # if way == 'Closest':

    # elif way == 'ClosestExpressed':

    # elif way == 'ActivityByContact':

    # else:
    #     raise ValueError('needs to be oneof Closest ClosestExpressed ActivityByContact')


# os.system()


def getCoLocalization():
    """
    for each annotations (super enhancer & TFgroups)
    for each TF, find highest peak/meanCov. if above thresh add to localization
    """


def refineGroupsWithHiC():
    """
    given HiC data, for each loops (should be less than X Mb distance. create a Xkb zone around both sides
    and find + merge each TF/TFgroup at this location)
    """


def getPeaksOverlap(peaks, isMerged=False, correlationMatrix=None, countMatrix=None, extend=1000,
                    onlyOn=[], quality={}):
    """
    generates a venn diagram of overlap between a replicate (quality anotation can be added)

    args:
    ----
      peakset: df of Bed Peak format
      onlyOn: list[string] of TF name
      quality: dict[string: {1,2,3}] a dict with a value between 1: good, 2: bad, 3: throw for each TF
    """


# def assignGene(peaks, bedFolder):
