#!/usr/bin/env python

import os
import sys
import ROOT

def checkPacked(args):
    sub=args

    toRun=[]
    nGood=0
    with open(sub,'r') as f:
        for line in f:
            if not 'arguments' in line : continue
            jobArgs=line.split()[2:]
            outF=jobArgs[2]
            if not os.path.isfile(outF):
                toRun.append((jobArgs,'missing'))
            else:
                try:
                    fIn=ROOT.TFile.Open(outF)
                    n=fIn.GetListOfKeys().GetSize()
                    if n<=1:
                        raise ValueError
                    fIn.Close()
                    nGood+=1
                except Exception as e:
                    toRun.append((jobArgs,'corrupted/no keys'))    

    print(sub,'good jobs:',nGood,'re-run these jobs:',len(toRun))
    for x in toRun:
        print('sh ${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/scripts/wrapAnalysis.sh %s'%(' '.join(x[0])))
        os.system('sh ${CMSSW_BASE}/src/HeavyIonsAnalysis/topskim/scripts/wrapAnalysis.sh %s'%(' '.join(x[0])))


task_list=[]
condor_dir=sys.argv[1]
for f in os.listdir(condor_dir):
    if not '.sub' in f : continue

    #in practice can't run in parallel because the local outputs are always named tt2l.root...
    checkPacked(os.path.join(condor_dir,f))
