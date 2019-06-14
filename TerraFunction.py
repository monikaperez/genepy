# Jeremie Kalfon
# jkobject@gmail.com
# for BroadInstitute CDS
# in 2019

import time


def wait_for_submission(wm, submissions):
	time = 0
  failed_submission = []
  if type(submissions) is type(""):
    submissions = [submissions]
  for submission_id in submissions:
  	submissions = wm.get_submission(submission_id)['workflows']
    	print("processing "+len(submissions["workflows"])+" submissions")
      for i in submissions["workflows"]:
        while i['status'] not in ['Done', 'Aborted', 'Failed', 'Succeeded']:
          time.sleep(60)
          time+=1
          print(time + " mn elapsed")
        if i["status"] in ['Failed', 'Aborted']:
          failed_submission.append(i["workflowEntity"]["entityName"])
  return failed_submission


def retry_failed(wm, submission_id):
	newset = wm.build_retry_set(submission_id)
	submission = wm.get_submission(submission_id)
	subid = wm.create_submission(submission['methodConfigurationName'], newset['name'],newset['type'], expression=newset['expression'])
	prevfailed = wait_for_submission(wm, submission_id)
	failed = wait_for_submission(wm, subid)
	notfailed = list(set(prevfailed)-set(failed))
	if len(notfailed) == 0:
		print("all failed")
		return failed
	else:
		samples = pd.Dataframe(wm.get_sample_sets().iloc[newset['name']])
		# add values to old failed set
		samples = wm.get_sample_sets
		# remove values of notfailed samples from newset with some still failed
		# upload them back again
		# return still failed ones
