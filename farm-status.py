#!/usr/bin/env python
import subprocess
import sys
import re

jobid = sys.argv[1]
# incomplete script
max_status_checks = 30
job_state = None
i = 1
while is.None(job_state):
    try:
        job_info = sp.check_output(shlex.split("scontrol -o show jobid {}".format(jobid)))
        job_state = re.search("JobState=(\w+)", job_info.decode())
        # do I need a .str() for states like RUNNING,FAILED
        break
    except sp.CalledProcessError as e:
        logger.error("scontrol process error for show jobid")
        logger.error(e)
        if i >= max_status_checks:
            job_state = "JOBSTATE COULD NOT BE DETERMINED"
            exit(0)
        else:
            i = i + 1
            time.sleep(1)

output = job_state #job_state.strip()
#output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if output == "COMPLETED":
  print("success")
elif any(r in output for r in running_status):
  print("running")
else:
  print("failed")
