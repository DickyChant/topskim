import os, sys, ROOT


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
                    if not os.path.getsize(outF):
                        print('file has size 0', outF)
                        raise ValueError
                    fIn=ROOT.TFile.Open(outF)
                    if not fIn:
                        print('file is a null pointer', outF)
                        raise ValueError
                    if (fIn.IsZombie()):
                        print('file is zombie', outF)
                        raise ValueError
                    n=fIn.GetListOfKeys().GetSize()
                    if n<=1:
                        print('no keys', outF)
                        raise ValueError
                    fIn.Close()
                    nGood+=1
                except ValueError as e:
                    toRun.append((jobArgs,'corrupted/no keys'))    

    print(sub,'good jobs:',nGood,'re-run these jobs:',len(toRun))
    return toRun
    ##for x in toRun:
    ##    os.system('sh scripts/wrapAnalysis.sh %s'%(' '.join(x[0])))

resubArgs = []

task_list=[]
condor_dir=sys.argv[1]
for f in os.listdir(condor_dir):
    if not '.sub' in f : continue
    if     '.swp' in f : continue
    if 'resubmit' in f : continue
    resubArgs += checkPacked(f)

if len(resubArgs):
    print('Condor resubmission')
    print('will resubmit this number of jobs', len(resubArgs))

    with open('condor_resubmit.sub','w') as c:
        c.write('executable  = %s/src/HeavyIonsAnalysis/topskim/scripts/wrapAnalysis.sh\n'%resubArgs[0][0][0])
        c.write('output      = condor_resubmit_$(ProcId).out\n')
        c.write('error       = condor_resubmit_$(ProcId).err\n')
        c.write('log         = condor_resubmit_$(ProcId).log\n')
        #c.write('requirements = (OpSysAndVer =?= "SLCern6")\n')
        c.write('+AccountingGroup = "group_u_CMST3.all"\n')
        c.write('+MaxRuntime = 14400\n')
    
        for res in resubArgs:
            c.write('\n')
            c.write('arguments   = {args} \n'.format(args= ' '.join(res[0])))
            c.write('queue 1\n')

    resub = input("do you want to resub these jobs to the batch? [y/n]")

    if resub in ['y','yes','Y','YES', 'Yes']:
        os.system('condor_submit condor_resubmit.sub')
    else:
        print('did not resubmit. the file is there if you want to do it by hand')

else:
    print(' all good. all files there!')
