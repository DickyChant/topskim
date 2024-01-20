from DataFormats.FWLite import Handle, Runs
import ROOT

lheruninfo=Handle('LHERunInfoProduct')

runs=Runs('root://cmsxrootd.fnal.gov///store/himc/HINPbPbAutumn18DR/TTJets_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/AODSIM/mva98_103X_upgrade2018_realistic_HI_v11-v1/260000/B9597CA0-3FCA-F841-BCE0-63BF3AD54B31.root') ## hin aMC@NLO
#runs=Runs('root://cmsxrootd.fnal.gov///store/himc/HINPbPbAutumn18DR/TT_TuneCP5_HydjetDrumMB_5p02TeV-powheg-pythia8/AODSIM/mva98_103X_upgrade2018_realistic_HI_v11-v1/270000/C860F47F-A00D-FC46-8F03-74E8674DAB62.root') ## hin powheg
#store/himc/RunIIpp5Spring18MiniAOD/TT_TuneCP5_5p02TeV-powheg-pythia8/MINIAODSIM/94X_mc2017_realistic_forppRef5TeV-v1/60000/04907A40-C3B7-E811-BA32-7CD30ACDE0CC.root ## pp
for r in runs:
    r.getByLabel('externalLHEProducer',lheruninfo)
    it=lheruninfo.product().headers_begin()
    while it!=lheruninfo.product().headers_end():
        lines=it.lines()
        allowPrint=False
        wgtCtr=0
        for i in range(0,lines.size()):
            linestr=lines.at(i)
            if '<weightgroup' in linestr : allowPrint=True
            if '</weightgroup' in linestr : allowPrint=False
            if not allowPrint : continue
            if 'weightgroup' in linestr :
                print('*'*50)
                print(linestr)
                print('*'*50)
            else:
                if not 'weight' in linestr : continue
                print(wgtCtr,linestr)
                wgtCtr+=1
        next(it)
