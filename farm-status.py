#!/usr/bin/env python
import subprocess
import sys
import re

#jobid = sys.argv[1]
farm_message = re.search("Submitted batch job (\d+)", sys.argv[1])
jobid = farm_message.group(1)

output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
if "COMPLETED" in output:
  print("success")
elif any(r in output for r in running_status):
  print("running")
else:
  print("failed")
