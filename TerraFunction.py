# Jeremie Kalfon
# jkobject@gmail.com
# for BroadInstitute CDS
# in 2019

import time


def wait_for_submission(wm, submissions):
  failed_submission = []
  if type(submissions) is type(""):
    submissions = [submissions]
  for submission_id in submissions:
    if not wm.get_submission(submission_id)['workflows']:
      while wm.get_submission(submission_id)['status'] not in {'Done', 'Aborted'}:
        time.sleep(60)
    else:
      for i in wm.get_submission(submission_id)["workflows"]:
        while i['status'] not in {'Done', 'Aborted'}:
          time.sleep(60)
        if i["status"] == 'Failed':
          print(newsample.loc[i["workflowEntity"]["entityName"]]["attr_sex"])
          print(i["workflowEntity"]["entityName"])
          failed_submission.append(i["workflowEntity"]["entityName"])
  # print and return well formated data
