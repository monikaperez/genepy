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


def saveOmicsOutput(vm, pathto_cnvpng='segmented_copy_ratio_img',
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
