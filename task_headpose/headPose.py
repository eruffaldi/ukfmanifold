import csv
import livestat
import sys

a = csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "headPose.txt","rb"))
o = csv.DictWriter(open(sys.argv[2] if len(sys.argv) > 2 else "headPoseout.csv","wb"),fieldnames=[x.strip() for x in a.fieldnames])
o.writeheader()
w = None
for x in a:
	allfields = a.fieldnames
	asarray = [x[y] for y in allfields]
	# find bad field
	idx = allfields.index(" pose_Rx")
	# rebuild the dictionary
	asarray[idx] = asarray[idx].strip().split(" ")
	asarray = asarray[0:idx] + asarray[idx] + asarray[idx+1:]
	d = dict([(allfields[i].strip(),float(asarray[i])) for i in range(0,len(allfields))])
	# create the statistics verifyier
	if w is None:
		w = dict([(k.strip(),livestat.LiveStat(k.strip())) for k in d.keys()])

	# update the statistics
	for k in d.keys():
		w[k].append(float(d[k]))

	print d["pose_Tx"],d["pose_Ty"],d["pose_Tz"],d["pose_Rx"],d["pose_Ry"],d["pose_Rz"]

	o.writerow(d)

for k in sorted(w.keys()):
	print w[k]