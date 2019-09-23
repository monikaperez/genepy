# Jeremie Kalfon
# jkobject@gmail.com
# for BroadInstitute CDS
# in 2019

import time
import pandas as pd
from google.cloud import storage
import pdb


def waitForSubmission(wm, submissions, raise_errors=True):
  failed_submission = []
  timing = 0
  assert submissions is not None
  if type(submissions) is type(""):
    submissions = [submissions]
  for scount, submission_id in enumerate(submissions):
    finished = False
    while not finished:
      done = 0
      failed = 0
      finished = True
      submission = wm.get_submission(submission_id)["workflows"]
      for wcount, i in enumerate(submission):
        if i['status'] not in {'Done', 'Aborted', 'Failed', 'Succeeded'}:
          finished = False
        elif i["status"] in {'Failed', 'Aborted'}:
          failed += 1
          if i["workflowEntity"]["entityName"] not in failed_submission:
            print(i["workflowEntity"]["entityName"])
            failed_submission.append(i["workflowEntity"]["entityName"])
        elif i["status"] in {'Done', 'Succeeded'}:
          done += 1
      if not finished:
        time.sleep(40)
        print("status is: Done for " + str(done) + " jobs in submission " + str(scount) + ". " + str(timing) + ",5 mn elapsed.", end="\r")
        timing += 1
        time.sleep(20)
        print("status is: Failed for " + str(failed) + " jobs in submission " + str(scount) + ". " + str(timing) + " mn elapsed.", end="\r")
      else:
        print(str(done / (done + failed)) + " of jobs Succeeded in submission " + str(scount) + ".")
  if len(failed_submission) > 0 and raise_errors:
    raise RuntimeError(str(len(failed_submission)) + " failed submission")
  return failed_submission
  # print and return well formated data


def uploadFromFolder(gcpfolder, prefix, wm, sep='_', updating=False, fformat="fastq12", newsamples=None, samplesetname=None):
  """
  upload samples (virtually: only creates tsv file) from a google bucket to a terra workspace

  it also creates a sample set.

  gcpfolder
  prefix: the folder path
  wm: the workspace terra
  sep: the separator (only takes the first part of the name before the sep character)
  updating: if needs
  fformat bambai, fastq12, fastqR1R2
  newsamples
  samplesetname
  """
  print('please be sure you gave access to your terra email account access to this bucket')
  if samplesetname is None:
    samplesetname = 'from:' + gcpfolder + prefix
  files = list_blobs_with_prefix(gcpfolder, prefix, '/')
  if fformat == "bambai":
    if newsamples is None:
      data = {'sample_id': [], 'bam': [], 'bai': []}
      for file in files:
        if val.split('.')[-1] in ["bam", "bai"]:
          name = file.split('/')[-1].split('.')[0].split(sep)[0][:-2]
          if name in data['sample_id']:
            pos = data['sample_id'].index(name)
            if file[-4:] == ".bam":
              data['bam'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif file[-4:] == ".bai":
              data['bai'].insert(pos, 'gs://' + gcpfolder + '/' + file)
          else:
            data['sample_id'].append(name)
            if file[-4:] == ".bam":
              data['bam'].append('gs://' + gcpfolder + '/' + file)
            elif file[-4:] == ".bai":
              data['bai'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
        else:
          print("unrecognized file type : " + file)
      df = pd.DataFrame(data)
      df = df.set_index("sample_id")
      df["participant"] = pd.Series(data['sample_id'], index=data['sample_id'])
      wm.upload_samples(df)
      wm.update_sample_set(samplesetname, df.index.values.tolist())
    else:
      # TODO: check if each column exists and can be added, else don't add it
      for i, val in enumerate(newsample["file_path"]):
        if val.split('/')[-1].split('.')[1] != "WholeGenome" or val.split('/')[-2] != "bam":
          newsample = newsample.drop(i)
        elif val.split('/')[1] != 'gs:':
          newsample["file_path"][i] = gcpfolder + newsample["file_path"][i].split('/')[-1]
      newsample = newsample.reset_index(drop=True)
      newsample = newsample.rename(index=str, columns={"sample_name": "sample_id", "subject_name": "participant_id", "file_path": "WGS_bam"})
      currfile = ""
      bai = [''] * int(newsample.shape[0])
      # creating an array of bai and adding it to their coresponding bams
      for i in newsample.index:
        currfile = newsample["WGS_bam"][i]
        if currfile.split('/')[-1].split('.')[-1] == "bai":
          bai[int(newsample[newsample["WGS_bam"] == currfile[:-4]].index.values[0])] = currfile
      newsample["WGS_bam_index"] = pd.Series(bai, index=newsample.index)
      # removing original bai rows
      for i in newsample.index:
        currfile = newsample["WGS_bam"][i]
        if currfile.split('/')[-1].split('.')[-1] == "bai":
          newsample = newsample.drop(i)
      newsample = newsample.reset_index(drop=True)
      newsample["sample_set"] = pd.Series([samplesetname] * int(newsample.shape[0]), index=newsample.index)
      newsample.set_index("sample_id", inplace=True, drop=True)
      newsample = newsample[newsample.columns.tolist()[1:] + [newsample.columns.tolist()[0]]]
      newsample = newsample.loc[~newsample.index.duplicated(keep='first')]
      newsample.to_csv("temp/samples.bambai.tsv", sep="\t")
      wm.upload_samples(newsample)
      wm.update_sample_set(samplesetname, newsample.index)
  if fformat in {"fastq12", "fastqR1R2"}:
    data = {'sample_id': [], 'fastq1': [], 'fastq2': []}
    # print and return well formated data
    for file in files:
      if file[-9:] == ".fastq.gz" or file[-6:] == ".fq.gz":
        name = file.split('/')[-1].split('.')[0].split(sep)[0][:-2]
        if name in data['sample_id']:
          pos = data['sample_id'].index(name)
          if fformat == "fastqR1R2":
            if "R1" in file:
              data['fastq1'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif "R2" in file:
              data['fastq2'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
          else:
            if file.split('.')[-3][-1] == '1':
              data['fastq1'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            elif file.split('.')[-3][-1] == '2':
              data['falsstq2'].insert(pos, 'gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq 1/2 error", file)
        else:
          data['sample_id'].append(name)
          if fformat == "fastqR1R2":
            if "R1" in file:
              data['fastq1'].append('gs://' + gcpfolder + '/' + file)
            elif "R2" in file:
              data['fastq2'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
          else:
            if file.split('.')[-3][-1] == '1':
              data['fastq1'].append('gs://' + gcpfolder + '/' + file)
            elif file.split('.')[-3][-1] == '2':
              data['fastq2'].append('gs://' + gcpfolder + '/' + file)
            else:
              raise Exception("No fastq R1/R2 error", file)
      else:
        print("unrecognized file type : " + file)
    df = pd.DataFrame(data)
    df = df.set_index("sample_id")
    df["participant"] = pd.Series(data['sample_id'], index=data['sample_id'])
    wm.upload_samples(df)
    wm.update_sample_set(samplesetname, df.index.values.tolist())


def updateAllSampleSet(wm, newsample_setname, Allsample_setname='All_samples'):
  """
  update the previous All Sample sample_set with the new samples that have been added.

  It is especially useful for the aggregate task
  """
  prevsamples = list(wm.get_sample_sets().loc[Allsample_setname]['samples'])
  newsamples = list(wm.get_sample_sets().loc[newsample_setname]['samples'])
  prevsamples.extend(newsamples)
  wm.update_sample_set(Allsample_setname, prevsamples)


def addToSampleSet(wm, samplesetid, samples):
  prevsamples = wm.get_sample_sets()['samples'][samplesetid]
  samples.extend(prevsamples)
  wm.update_sample_set(samplesetid, samples)


def addToPairSet(wm, pairsetid, pairs):
  prevpairs = wm.get_pair_sets()[pairsetid].pairs.tolist()
  pairs.extend(prevpairs)
  wm.update_sample_set(pairsetid, list(set(pairs)))


def list_blobs_with_prefix(bucket_name, prefix, delimiter=None):
  """Lists all the blobs in the bucket that begin with the prefix.

  This can be used to list all blobs in a "folder", e.g. "public/".

  The delimiter argument can be used to restrict the results to only the
  "files" in the given "folder". Without the delimiter, the entire tree under
  the prefix is returned. For example, given these blobs:

      /a/1.txt
      /a/b/2.txt

  If you just specify prefix = '/a', you'll get back:

      /a/1.txt
      /a/b/2.txt

  However, if you specify prefix='/a' and delimiter='/', you'll get back:

      /a/1.txt

  """
  storage_client = storage.Client()
  bucket = storage_client.get_bucket(bucket_name)
  ret = []
  blobs = bucket.list_blobs(prefix=prefix, delimiter=delimiter)
  for blob in blobs:
    ret.append(blob.name)
  return(ret)


def saveOmicsOutput(wm, pathto_cnvpng='segmented_copy_ratio_img',
                    pathto_stats='sample_statistics',
                    specific_cohorts=[],
                    speicifc_celllines=[],
                    is_from_pairs=True,
                    pathto_snv='filtered_variants',
                    pathto_seg='cnv_calls',
                    datadir='gs://cclf_results/targeted/kim_sept/',
                    specific_samples=[]):
  if specific_cohorts:
    samples = wm.get_samples()
    samples = samples[samples.index.isin(specificlist)]
  if is_from_pairs:
    pairs = wm.get_pairs()
    pairs = pairs[pairs['case_sample'].isin(specificlist)]
  for i, val in samples.iterrows():
    os.system('gsutil cp ' + val[pathto_seg] + ' ' + datadir + i + '/')
    os.system('gsutil cp ' + val[pathto_cnvpng] + ' ' + datadir + i + '/')
    os.system('gsutil cp ' + val[pathto_stats] + ' ' + datadir + i + '/')
    if is_from_pairs:
      snvs = pairs[pairs["case_sample"] == i][pathto_snv]
      for snv in snvs:
        if snv is not np.nan:
          os.system('gsutil cp ' + snv + ' ' + datadir + i + '/')
          break
    else:
      os.system('gsutil cp ' + val[pathto_snv] + ' ' + datadir + i + '/')


# ["gs://fc-ed07d172-7980-475e-81de-108a694a3532/",
#  "gs://fc-c23078b3-05b3-4158-ba8f-2b1eeb1bfa16/",
#  "gs://fc-51050008-201e-4a40-8ec7-28b6fb2b1885/"]
# "gs://fc-secure-98816a9e-5207-4361-8bf0-f9e046966e62/"
def changeGSlocation(wmfrom, wmto=None, prevgslist=[], newgs='', index_func=None):
  data = {}
  try:
    a = wmfrom.get_participants()
    data.update({'participants': a})
  except:
    print('no participants')
  try:
    a = wmfrom.get_samples()
    data.update({'samples': a})
  except:
    print('no samples')
  try:
    a = wmfrom.get_pair_sets()
    data.update({'pair_sets': a})
  except:
    print('no pair_sets')
  try:
    a = wmfrom.get_pairs()
    data.update({'pairs': a})
  except:
    print('no pairs')
  try:
    a = wmfrom.get_sample_sets()
    data.update({'sample_sets': a})
  except:
    print('no sample_sets')
  # currently works only for sample, sample
  for k, entity in data.iteritems():
    for k, val in entity.iterrows():
      for i, v in enumerate(val):
        if type(v) is str:
          for prev in prevgslist:
            v = v.replace(prev, newgs)
          val[i] = v
        if type(v) is list:
          ind = []
          for j in v:
            for prev in prevgslist:
              j = j.replace(prev, newgs)
            ind.append(j)
          val[i] = ind
        entity.loc[k] = val
    if wmto is None:
      wmto = wmfrom
    if "participants" in data:
      wmto.upload_participants(data['participants'])
    if "samples" in data:
      wmto.upload_samples(data['samples'])
    if "pairs" in data:
      wmto.upload_pairs(data['pairs'])
    if "pair_set" in data:
      pairset = data['pair_set'].drop('pairs', 1)
      wmto.upload_entities('pair_set', pairset)
      for i, val in data['pair_set'].iterrows():
        wmto.update_pair_set(i, val.pairs)
    if "sample_set" in data:
      sampleset = data['sample_set'].drop('samples', 1)
      wmto.upload_entities('sample_set', sampleset)
      for i, val in data['sample_set'].iterrows():
        wmto.update_sample_set(i, val.samples)


def renametsvs(wmfrom, wmto=None, index_func=None):
  data = {}
  try:
    a = wmfrom.get_participants()
    data.update({'participants': a})
  except:
    print('no participants')
  try:
    a = wmfrom.get_samples()
    data.update({'samples': a})
  except:
    print('no samples')
  try:
    a = wmfrom.get_pair_sets()
    data.update({'pair_sets': a})
  except:
    print('no pair_sets')
  try:
    a = wmfrom.get_pairs()
    data.update({'pairs': a})
  except:
    print('no pairs')
  try:
    a = wmfrom.get_sample_sets()
    data.update({'sample_sets': a})
  except:
    print('no sample_sets')
  # currently works only for sample, sample
  for k, entity in data.iteritems():
    ind = []
    for i in entity.index:
      pos = val.find('-SM')
      if pos != -1:
        val = val[pos + 1:]
      pos = val.find('-SM')
      if pos != -1:
        val = val[:9] + val[pos + 1:]
      ind.append(val)
    entity.index = ind
    for k, val in entity.iterrows():
      for i, v in enumerate(val):
        if type(v) is list or type(v) is str:
          ind = []
          for j in v:
            pos = j.find('-SM')
            if pos != -1:
              j = j[pos + 1:]
            pos = j.find('-SM')
            if pos != -1:
              j = j[:9] + j[pos + 1:]
            ind.append(j)
          val[i] = ind
        entity.loc[k] = val
    if wmto is None:
      wmto = wmfrom
    if "participants" in data:
      wmto.upload_participants(data['participants'])
    if "samples" in data:
      wmto.upload_samples(data['samples'])
    if "pairs" in data:
      wmto.upload_pairs(data['pairs'])
    if "pair_set" in data:
      pairset = data['pair_set'].drop('pairs', 1)
      wmto.upload_entities('pair_set', pairset)
      for i, val in data['pair_set'].iterrows():
        wmto.update_pair_set(i, val.pairs)
    if "sample_set" in data:
      sampleset = data['sample_set'].drop('samples', 1)
      wmto.upload_entities('sample_set', sampleset)
      for i, val in data['sample_set'].iterrows():
        wmto.update_sample_set(i, val.samples)
