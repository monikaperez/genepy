

def extractPairedSingleEndFrom(folder, pattern='R1/R2', sep='-', namepos=2):
  # given a folder, find using a specific patterns, single end and paired end files
  # return a list of single end files and a list of list of paired end files [[R1,R2]]
  single = []
  paired = {}
  for val in os.listdir(folder):
    if val.contains('R1'):
      name = val.split(sep)[namepos]
      paired[name]['R1'] = val
    elif val.contains('R2'):
      name = val.split(sep)[namepos]
      paired[name]['R2'] = val
    else:
      single.append(val)

  return single, pd.Dataframe(paired)


def findReplicates(folder, sep='-', namings='-r([0-9])', pos=2):


"""
creates a dict of name and replicate files
"""
  rep = {}
  for val in os.listdir(folder):
    number = re.search(namings, val):
    if number:
      name = val.split(sep)[namepos]
      if rep.has_key(name):
        rep[name].append(val)
      else:
        rep[name] = val

  return rep


def computeSingleEnd(singlend, folder="data/seqs/", numthreads=8, peaksFolder="peaks/",
                     ismapped=False, mappedFolder='mapped/', refFolder='data/reference/index'):
  # run the singleEnd pipeline
  # for alignment etc, one can use pysam ready made implementation of samtools
  for val in singlend:
    out1 = folder + mappedFolder + val.split('.')[0] + ".mapped.sam"
    if not ismapped:
      in1 = folder + val
      os.system("bowtie2 -x " + refFolder + " --threads " + str(numthreads) + " -t -k 1 --very-sensitive -U " + in1 + " -S " + out1)
    out2 = folder + peaksFolder + val.split('.')[0]
    print(out1)
    os.system("macs2 callpeak -f SAM -t " + out1 + " --outdir " + out2)
    # it can take many TB so better delete


def computePairedEnd(pairedend, folder="data/seqs/", numthreads=8, peaksFolder="peaks/",
                     ismapped=False, mappedFolder='mapped/', refFolder='data/reference/index'):
  # run the paired end pipeline
  for key, val in pairedend.iteritems():
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
  # run the bigwig command line for a set of bam files in a folder
  size = {"GRCh37": '2864785220',
          "GRCh38": '2913022398'}
  for i in bams:
    in1 = folder + i
    out1 = folder + "bigwig/" + i + split('.')[0] + '.bw'
    os.system("bamCoverage --effectiveGenomeSize " + size[genome] + " -p " + str(numthreads) + " -b " + in1 + "-of bigwig -o " + out1)


def mergeBams(rep):


"""
uses samtools to merge a set of replicates considered into one file
"""
for i, val in rep.iteritems():
  out1 = i + '.merged.bam'
  for bam in val:
    in1 += ' ' + bam
  os.system("samtools merge " + out1 + in1)


def loadNarrowPeaks():
  # for each data peak type file listed in a MACS2 way in a given MACS2 output folder, will merge them
  # all into one dataframe and output the dataframe
  bindings = pd.DataFrame()
  for folder in os.listdir("data/peaks"):
    if any(tf in folder for tf in CTF):
      binding = pd.read_csv("data/peaks/" + folder + "/NA_peaks.narrowPeak", sep='\t', header=None)
      binding['type'] = [folder] * binding.shape[0]
      bindings = bindings.append(binding)
  bindings = bindings.drop(5, 1)
  bindings = bindings.rename(columns={
      0: "chromosome",
      1: 'start',
      2: 'end',
      3: 'peak_number',
      4: 'size',
      6: "foldchange",
      7: "-log10pvalue",
      8: "-log10qvalue",
      9: 'relative_summit_pos'})
  bindings = bindings.sort_values(by=["chromosome", "start", "end"], axis=0)
  bindings.start = bindings.start.astype('int')
  bindings.end = bindings.end.astype('int')
  bindings.size = bindings.size.astype('int')
  bindings.relative_summit_pos = bindings.relative_summit_pos.astype('int')
  bindings.foldchange = bindings.foldchange.astype('float')
  bindings["-log10pvalue"] = bindings["-log10pvalue"].astype('float')
  bindings['-log10qvalue'] = bindings['-log10qvalue'].astype('float')
  return bindingss


def computePeaksAt(peaks, bams, folder='data/seqs', window=1000, numthreads=8):

  # get pysam data
  # ask for counts only at specific locus based on windows from center+-size from sorted MYC peaks
  # for each counts, do a rolling average (or a convolving of the data) with numpy
  # append to an array
  # return array, normalized
  loaded = {}
  res = {i: np.zeroes((len(peaks), window)) for i in bams}
  for val in bams:
    loaded.update({val: pysam.AlignmentFile(val, 'rb', threads=numthreads)})

  for val in peaks.iteritems():
    for k, bam in loaded.iteritems():
      center = int((val['Start'] + val['End']) / 2)
      for pileupcolumn in bam.pileup(contig=val['chromosome'], start=center - window,
                                     stop=center + window - 1, truncate=True):
        res[k][(center - window) - pileupcolumn.pos] = pileupcolumn.n
  return res


def computeCov(bamfolder, ref="data/reference/genome.fa"):
  for val in os.listdir(bamfolder):
    os.system("bedtools genomecov -ibam " + bamfolder + val)


def mergeReplicatePeaks(peaks, reps, bamfolder, avgCov, window=200, numthread=8):
  """
  /!/ should only be passed peaks with at least one good replicate

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

  args:
  ----
  peaks: df[bed-like] all the peaks into the sameBam with a column containing the 'type' 
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
  tfs = list(set(peaks['type'].tolist()))
  totpeaknumber = 0
  mergedpeaksdict = {}
  remove = []
  for k, rep in reps:
    print("doing peak " + k)
    peaknumber = 0
    merged_bed = {
        "chromosome": [tfpeaks.iloc[1]['chromosome']],
        "start": [tfpeaks.iloc[1]['start']],
        "end": [],
        "peak_number": [peaknumber + totpeaknumber],
        "size": [],
        "foldchange": [],
        "log10pvalue": [],
        "log10qvalue": [],
        "relative_summit_pos": [tfpeaks.iloc[1]['relative_summit_pos']]
    }
    foldchange = [tfpeaks.iloc[1]['foldchange']]
    log10pvalue = [tfpeaks.iloc[1]['log10pvalue']]
    log10qvalue = [tfpeaks.iloc[1]['log10qvalue']]
    relative_summit_pos = tfpeaks.iloc[1]['relative_summit_pos']
    # computes overlap by extending a bit the window (100bp?) should be ~readsize
    tfpeaks = peaks[peaks['type'] == tf]
    prev_end = -101
    prevchrom = 1
    tfmerged = {a: [0] for a in val}
    for i, pos, peak in enumerate(tfpeaks.iloc[1:].iterrows()):
      print(i / len(tfpeaks), end='\r')
      if prev_end + window > peak['start'] and prevchrom == peak['chromosome']:
        tfmerged[peak['type']] = 1
        foldchange.append(peak['foldchange'])
        log10pvalue.append(peak['log10pvalue'])
        log10qvalue.append(peak['log10qvalue'])
        if peak['foldchange'] > max(foldchange):
          relative_summit_pos = peak['relative_summit_pos']
      else:
        # newpeak
        for k, val in tfmerged.iteritems():
          val.append(0)
        peaknumber += 1
        merged_bed['chromosome'].append(peak['chromosome'])
        merged_bed['start'].append(peak['start'])
        merged_bed['end'].append(prev_end)
        merged_bed['peak_number'].append(peaknumber + totpeaknumber)
        merged_bed['size'].append(prev_end - merged_bed['start'][-2])
        merged_bed['foldchange'].append(mean(foldchange))
        merged_bed['log10pvalue'].append(max(log10pvalue))
        merged_bed['log10qvalue'].append(max(log10qvalue))
        merged_bed['relative_summit_pos'].append(relative_summit_pos)
        foldchange = [peak['foldchange']]
        log10pvalue = [peak['log10pvalue']]
        log10qvalue = [peak['log10qvalue']]
        relative_summit_pos = peak['relative_summit_pos']
      prev_end = peak['end']
      prev_chrom = peak['end']
    totpeaknumber += peaknumber
    merged_bed = pd.Dataframe(merged_bed)
    finalpeaks = pd.Dataframe(tfmerged)
    print('finish first overlaps lookup')
    # flag when biggest is <1000 peaks
    if len(finalpeaks) < 1000:
      print('!TF has less than 1000 PEAKS!')
    # for each TF (replicates), compute number of peaks
    peakmatrix = finalpeaks.values

    # compute overlap matrix (venn?)
    if len(peakmatrix) < 7 and doPlot:
      inp = []
      for peakpres in peakmatrix.T:  # https://github.com/tctianchi/pyvenn
        inp.append([i for i, val in enumerate(peakpres)])
      venn(inp, finalpeaks.columns)
    else print('too many replicates for Venn')

    # compute smilarity matrix (plot?)
    corrmat = np.corrcoef(peakmatrix)
    if doPlot:
      Helper.plotCorrelationMatrix(corrmat, names=finalpeaks.columns,
                                   title("Correlation matrix of " + k + " replicates"))
    totpeaks = np.sum(peakmatrix)
    biggest_ind = np.argsort(totpeaks)[0]
    corr2big = corrmat[biggest_ind]
    biggest = finalpeaks[biggest_ind]

    # starts with highest similarity and go descending
    for i, val in enumerate(np.argsort(corr2big)):
      # if avg non overlap > 60%, and first, and none small flag TF as unreliable.
      peakname = finalpeaks[val]
      if corr2big[val] < 0.7:
        if totpeak[val] > totpeak[biggest_ind] * 0.5
          if i == 0:
            print("Wrong TF")
            remove.append(k)
            break
          # if not first, throw the other replicate and continue
          finalpeaks.drop(peakname, 1)
        # if small and small overlaps more than 80% do findAdditionalPeaks only on for small else throw small
        elif np.sum(np.logical_and(peakmatrix[biggest_ind], peakmatrix[val])) / np.sum(peakmatrix[biggest_ind]) > 0.8:
          tolookfor = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val] == 0)
          additionalpeaks = findAdditionalPeaks(merged_bed, tolookfor, peakname, bamfolder,
                                                avgCov[peakname], window, numthread)

          # for new set of peaks >0.7 size of big, flag for merge bams
          if np.sum(np.logical_and(peakmatrix[biggest_ind], newpeakset)) / np.sum(peakmatrix[biggest_ind]) > 0.8
            tomergebam.append([biggest, peakname])
        else:
          finalpeaks.drop(peakname, 1)

      # else findAdditionalPeaks:
      else:
        tolookfor = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val] == 0)
        additionalpeaks = findAdditionalPeaks(merged_bed, tolookfor, peakname, bamfolder,
                                              avgCov[peakname], window, numthread)
        peakmatrix[val] = np.logical_or(peakmatrix[val], additionalpeaks)
        tolookfor = np.logical_and(peakmatrix[val], peakmatrix[biggest_ind] == 0)
        additionalpeaks = findAdditionalPeaks(merged_bed, tolookfor, biggest, bamfolder,
                                              avgCov[biggest], window, numthread)
        peakmatrix[biggest_ind] = np.logical_or(peakmatrix[biggest_ind], additionalpeaks)
        tomergebam.append([biggest, peakname])  # flag for merge bams
        # take the intersection of both peaksets
        peakmatrix[biggest_ind] = np.logical_and(peakmatrix[biggest_ind], peakmatrix[val])
        peakmatrix.drop(val, 1)

    mergedpeaksdict.update({k: finalpeaks})


def findAdditionalPeaks(peaks, tolookfor, filetofindin, folder, avgCov, window=200, numthreads=8,):
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
  bam = pysam.AlignmentFile(folder + filetofindin, 'rb', threads=numthreads)
  res = np.zeroes(len(peaks))
  for i, val in enumerate(tolookfor):
    if val:
      pos = peaks.iloc[i]
      center = int((pos['Start'] + pos['End']) / 2)
      zone = np.zeroes(window * 2)
      for pileupcolumn in bam.pileup(contig=pos['chromosome'], start=center - window,
                                     stop=center + window - 1, truncate=True):
        zone[(center - window) - pileupcolumn.pos] = pileupcolumn.n
      if max(zone) > 10 * avgCov:
        res[i] = 1

  return res
  # BETTER WAY
  #
  # for each nucleotide in the sequence not in registered peaks,
  # create a count distribution accross nucleotide use bedgraph file output from
  # "bedtools genomecov -d -ibam youralignment.bam" command
  # number of times you see 0,1,2,3,4,5.. reads for each nucleotide. fit an expo distribution to it
  # if above a p-value of .1, mark as found, return the new marks with their p-value


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
"""


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


def getPeaksOverlap(peaks, , isMerged=False, correlationMatrix=None, countMatrix=None, extend=1000, onlyOn=[], quality={}):

"""
generates a venn diagram of overlap between a replicate (quality anotation can be added)

args:
----
  peakset: df of Bed Peak format
  onlyOn: list[string] of TF name
  quality: dict[string: {1,2,3}] a dict with a value between 1: good, 2: bad, 3: throw for each TF
"""


def assignGene(peaks, bedFolder):
