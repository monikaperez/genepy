import dalmatian as dm
import pandas as pd
import Helper
import os
import numpy as np


def getReport(workspace1="CCLF_TSCA_2_0_2", namespace1="nci-mimoun-bi-org",
              pathto_cnvpng='segmented_copy_ratio_img',
              pathto_stats='sample_statistics',
              pathto_snv='filtered_variants',
              pathto_seg='cnv_calls',
              workspacewes='CCLF_WES', namespacewes='nci-mimoun-bi-org',
              pathto_cnvpng_wes='segmented_copy_ratio_img',
              pathto_stats_wes='sample_statistics',
              pathto_snv_wes='mafLite',
              pathto_seg_wes='tumor_seg',
              is_from_pairs=True,
              datadir='gs://cclf_results/targeted/kim_sept/',
              tempdir='temp/cclfmerge/',
              specificlist=None
              ):
  """
  """
  # if we get many condition: merge and display the one that has been selected)
  # if we get WES data: add on the side as additional data and merge CNplot
  # if we get primary: add on the side as additional data and merge CNplot
  print('you need to be on macOS for now')
  wm1 = dm.WorkspaceManager(namespace1, workspace1)
  wm_wes = dm.WorkspaceManager(namespacewes, workspacewes)
  if type(specificlist) is str:
    # we consider it is a filename
    specificlist = pd.read_csv(specificlist).tolist()
  elif specificlist is None:
    print("taking all samples from the workspace")
    specificlist = wm.get_sample_sets().index

  sample = wm1.get_samples()
  sample_wes = wm_wes.get_samples()
  if is_from_pairs:
    pair = wm1.get_pairs()
  sample_part = sample.participant.tolist()
  samplewes_part = sample_wes.participant.tolist()
  for val in specificlist:
    found = False
    mutfile = pd.DataFrame()
    cnfile = pd.DataFrame()
    images = []

    if val in sample_part:
      sample = sample[sample.participant == val]
      for k, condition in sample.iterrows():
        if condition['sample_type'] == "Tumor":
          found = True
          cond_name = condition['media']
          outputloc = datadir + condition['primary_disease'].replace(' ', '_') + '/' + val + '/'
          os.system('gsutil cp ' + condition[pathto_seg] + ' ' + tempdir + 'copy_number.tsv')
          os.system('gsutil cp ' + condition[pathto_cnvpng] + ' ' + tempdir + cond_name + '_copy_number_map.png')
          images.append(tempdir + cond_name + '_copy_number_map.png')
          os.system('gsutil cp ' + condition[pathto_stats] + ' ' + outputloc + cond_name + '_sample_statistics.txt')
          if is_from_pairs:
            snvs = pair[pair["case_sample"] == k][pathto_snv]
            for snv in snvs:
              if snv is not np.nan:
                os.system('gsutil cp ' + snv + ' ' + tempdir + 'mutations.tsv')
                break
          else:
            os.system('gsutil cp ' + condition[pathto_snv_unpaired] + ' ' + tempdir + 'mutations.tsv')
          mut = pd.read_csv(tempdir + 'mutations.tsv', sep='\t')
          mut['condition'] = cond_name
          mutfile = mutfile.append(mut)
          cn = pd.read_csv(tempdir + 'copy_number.tsv', sep='\t')
          cn['condition'] = cond_name
          cnfile = cnfile.append(cn)
      if found:
        cnfile.to_csv(tempdir + 'cn.tsv', sep='\t')
        mutfile.to_csv(tempdir + 'mut.tsv', sep='\t')
        os.system('gsutil cp ' + tempdir + 'cn.tsv ' + outputloc + 'copy_number.tsv')
        os.system('gsutil cp ' + tempdir + 'mut.tsv ' + outputloc + 'mutation.tsv')
    if val in samplewes_part:
      sample_wes = sample_wes[sample_wes.participant.isin(specificlist)]
      for k, wes in sample_wes.iterrows():
        if wes['sample_type'] == "Tumor":
          primary = wes['primary_disease'].replace(' ', '_') if 'primary_disease' in wes else 'unknown'
          outputloc = datadir + primary + '/' + val + '/'
          found = True
          os.system('gsutil cp ' + wes[pathto_seg_wes] + ' ' + outputloc + 'wes_copy_number.tsv')

          os.system('gsutil cp ' + wes[pathto_cnvpng_wes] + ' ' + tempdir + 'wes_copy_number_map.pdf')
          os.system('sips -s format png ' + tempdir + 'wes_copy_number_map.pdf --out ' + tempdir + 'wes_copy_number_map.png')
          images.append(tempdir + 'wes_copy_number_map.png')

          if is_from_pairs:
            snvs = pair[pair["case_sample"] == k]
            if pathto_snv_wes in snvs:
              snvs = snvs[pathto_snv_wes]
            else:
              continue
            for snv in snvs:
              if snv is not np.nan:
                os.system('gsutil cp ' + snv + ' ' + outputloc + 'mutations.tsv')
                break
          else:
            os.system('gsutil cp ' + wes[pathto_snv_wes] + ' ' + outputloc + 'wes_mutations.tsv')

    if not found:
      print("we did not found any data for " + val)
    else:
      Helper.mergeImages(images, tempdir + 'merged.png')
      os.system('gsutil cp ' + tempdir + 'merged.png ' + outputloc + 'merged_copy_number_map.png')
