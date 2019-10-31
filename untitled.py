  for k, val in samp.iterrows():
    if val[bamcol] != 'NA' and val[baicol] != 'NA' and val[bamcol] is not None:
    # if bai has duplicate size
      new = val[bamcol]
      print(new)
      prev = new
      if newgs not in new:
        for prevgs in prevgslist:
          new = new.replace(prevgs, newgs)
        if flag_non_matching:
          if 'gs://' == prev[:5]:
            if new == prev:
              flaglist.append(prev)
              continue
      code = os.system('gsutil ls ' + new)
      if code==0:
        print('just changing name')
        samp.loc[k][bamcol]=new
        samp.loc[k][baicol]=new.replace('.bam','.bai') if os.system('gsutil ls ' + 
          new.replace('.bam','.bai')) == 0 else new.replace('.bam','.bam.bai')

      elif code == 256:
        print('not found')
        bai = new.replace('.bam','.bai') if os.system('gsutil ls ' + 
          new.replace('.bam','.bai')) == 0 else new.replace('.bam','.bam.bai')
        for va in names[sizes[bai]]:
          # for all duplicate size
          # if ls bam of bai duplicate size work
          # mv bam to bampath in folder
          if '.bam' in va:
            if os.system('gsutil ls ' + va.split('.bam.bai')[0] + ".bam") == 0:
              if va.split('.bam.bai')[0] + ".bam" in previous:
                print(str(k)+' was also found in wes... skipping '+str(va))
                continue
              print("refound")
              samp.loc[k,bamcol]=va.split('.bam.bai')[0] + ".bam"
              samp.loc[k,baicol]=va
              break
          elif os.system('gsutil ls ' + va.split('.bai')[0] + ".bam") == 0:
            if va.split('.bai')[0] + ".bam" in previous:
              print(str(k)+' was also found in wes... skipping '+str(va))
              continue
            print("refound")
            samp.loc[k,bamcol]=va.split('.bai')[0] + ".bam"
            samp.loc[k,baicol]=va
            break
      elif code == signal.SIGINT:
        print('Awakened')
        break
    else:
      print("no data for " + str(k))