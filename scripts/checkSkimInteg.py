import os,sys
import ROOT

work='%s/src/HeavyIonsAnalysis/topskim'%os.environ['CMSSW_BASE']
script='%s/scripts/skimScript.sh'%work
condor_file=sys.argv[1]
with open(condor_file,'r') as condor:
    lines=condor.readlines()
    args=[l.split('=')[1].strip() for l in lines if 'arguments' in l]
    header=[l.strip() for l in lines if l.find('arguments')!=0 and l.find('queue')!=0]
print(header)
input()
    
redo=[]
for ia in args:
    a=ia.split()
    
    #check that output was produced
    outF=os.path.join(a[-1],a[-2])
    print(outF, end=' ')
    if not os.path.isfile(outF):
        redo.append(a)
        print('not found')
        continue

    #check it was correctly closed
    outF=outF.replace('/eos/cms/','root://eoscms//')
    try:
        inF=ROOT.TFile.Open(outF)
        keys=inF.GetListOfKeys()
        evts=inF.Get('Events')
        nentries=evts.GetEntriesFast()
        if inF.IsZombie() or keys.GetSize()==0 or nentries==0:
            raise Exception('corrupted file')
        inF.Close()
    except:
        redo.append(a)
        print('is corrupted')
        continue

    print('is good with',nentries,'events')

with open(condor_file.replace('.condor','_resub.condor'),'w') as condor:
    for l in header:
        if len(l)==0 : continue
        condor.write(l+'\n')
    for ia in redo:
        condor.write('arguments = %s\n'%' '.join(ia))
        condor.write('queue 1\n')

#local run
#for ia in redo:
#    print 'Reprocessing locally',ia
#    cfg=ia[1]
#    os.system("sed -i 's/file:\/eos\/cms\/store/\/store/g' %s"%cfg)
#    os.system("sh %s %s"%(script,' '.join(ia)))
