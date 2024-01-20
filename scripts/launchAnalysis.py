import ROOT
import os,sys
from subprocess import Popen, PIPE

if len(sys.argv)<4:
    print("python scripts/launchAnalysis.py input output tag [extraOpts]")
    exit(-1)

input=sys.argv[1]
output=sys.argv[2]
tag=sys.argv[3]
cmssw=os.environ['CMSSW_BASE']
extraOpts=' '.join(sys.argv[4:]) if len(sys.argv)>4 else ''

#prepare output
os.system('mkdir -p %s'%output)

#create condor jobs
print('Condor submission script is condor_%s.sub'%tag)
with open('condor_%s.sub'%tag,'w') as c:
    c.write('executable  = %s/src/HeavyIonsAnalysis/topskim/scripts/wrapAnalysis.sh\n'%cmssw)
    c.write('output      = condor_%s.out\n'%tag)
    c.write('error       = condor_%s.err\n'%tag)
    c.write('log         = condor_%s.log\n'%tag)
    #c.write('requirements = (OpSysAndVer =?= "SLCern6")\n')
    c.write('+AccountingGroup = "group_u_CMST3.all"\n')
    c.write('+MaxRuntime = 10800\n')

    chunks=os.listdir(input)
    for i in range(len(chunks)):
        c.write('arguments   = {0} {1}/{2} {3}/{4}_{5}.root {6}\n'.format(cmssw,input,chunks[i],output,tag,i,extraOpts))
        c.write('queue 1\n')

#submit to condor
os.system('condor_submit condor_%s.sub'%tag)
