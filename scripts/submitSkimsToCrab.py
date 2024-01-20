#!/usr/bin/env python

import os

def createJob(dataset,pset,lumiMask,lfnDirBase,removePrevCrab=True):
        
    """outputs a crab cfg file to submit the job"""

    request=dataset.split('/')[1]

    if removePrevCrab:
        os.system('rm -rf grid/crab_%s'%request) 

    cfg_file='grid/%s_cfg.py'%request
    cfg=open(cfg_file,'w')
    cfg.write('from WMCore.Configuration import Configuration\n')
    cfg.write('import os\n')
    cfg.write('config = Configuration()\n')
    cfg.write('\n')
    cfg.write('config.section_("General")\n')
    cfg.write('config.General.requestName = "%s"\n' % request)
    cfg.write('config.General.workArea = "grid"\n')
    cfg.write('config.General.transferOutputs=True\n')
    cfg.write('\n')
    cfg.write('config.section_("JobType")\n')
    cfg.write('config.JobType.pluginName = "Analysis"\n')
    cfg.write('config.JobType.psetName = "%s"\n'%pset)
    cfg.write('config.JobType.disableAutomaticOutputCollection = False\n')
    #cfg.write('config.JobType.pyCfgParams = [\'isPP=False\',\'maxEvents=-1\',\'outputFile=HiForest.root\']\n')
    cfg.write('config.JobType.outputFiles = [\'HiForestAOD.root\']\n')
    cfg.write('config.JobType.allowUndistributedCMSSW = True\n')
    cfg.write('\n')
    cfg.write('config.section_("Data")\n')
    cfg.write('config.Data.inputDataset = "%s"\n' % dataset)
    cfg.write('config.Data.inputDBS = "global"\n')
    cfg.write('config.Data.splitting = "LumiBased"\n')
    cfg.write('config.Data.unitsPerJob = 10\n')
    if lumiMask:
        cfg.write('config.Data.lumiMask = \'%s\'\n' %lumiMask)
    cfg.write('config.Data.ignoreLocality = False\n')    
    cfg.write('config.Data.publication = False\n')
    cfg.write('config.Data.outLFNDirBase = \"%s\"\n' % lfnDirBase)
    cfg.write('\n')
    cfg.write('config.section_("Site")\n')
    cfg.write('config.Site.storageSite = "T2_CH_CERN"\n')
    cfg.close()
    
    return cfg_file


#prepare output
os.system('mkdir -p grid')
cmssw=os.environ['CMSSW_BASE']
pset={False:"%s/src/HeavyIonsAnalysis/JetAnalysis/test/runForestAOD_pponAA_MIX_103X.py"%cmssw,
      True:"%s/src/HeavyIonsAnalysis/JetAnalysis/test/runForestAOD_pponAA_DATA_103X.py"%cmssw}
lumiMask={False:None,
          True:"/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/HI/PromptReco/Cert_326381-327564_HI_PromptReco_Collisions18_JSON.txt"}
lfnDirBase="/store/group/cmst3/group/hintt/grid_17Jul/"
submit=True

if submit:
    print('Will submit crab jobs. Make sure you have sourced crab.sh')
    print('source /cvmfs/cms.cern.ch/crab3/crab.sh')

for dset,isData in [
        #('/WJetsToLNu_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/TT_TuneCP5_HydjetDrumMB_5p02TeV-powheg-pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/TTJets_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_HydjetDrumMB_5p02TeV-powheg-pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/QCD_Pt_80to170_bcToE_TuneCP5_HydjetDrumMB_5p02TeV_pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v2/AODSIM',False),
        #('/QCD_Pt_20to30_bcToE_TuneCP5_HydjetDrumMB_5p02TeV_pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/QCD_Pt-20to30_EMEnriched_TuneCP5_HydjetDrumMB_5p02TeV_pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/QCD_Pt-30to50_EMEnriched_TuneCP5_HydjetDrumMB_5p02TeV_pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/QCD_Pt-50to80_EMEnriched_TuneCP5_HydjetDrumMB_5p02TeV_pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        #('/QCD_Pt-80to120_EMEnriched_TuneCP5_HydjetDrumMB_5p02TeV_pythia8/HINPbPbAutumn18DR-mva98_103X_upgrade2018_realistic_HI_v11-v1/AODSIM',False),
        ('/HIHardProbes/HIRun2018A-04Apr2019-v1/AOD',True),
        ('/HISingleMuon/HIRun2018A-04Apr2019-v1/AOD',True)
    ]:

    cfg=createJob(dataset=dset,pset=pset[isData],lumiMask=lumiMask[isData],lfnDirBase=lfnDirBase)
    if submit : 
        print('Submitting',cfg)
        os.system('alias crab=\'/cvmfs/cms.cern.ch/crab3/crab-env-bootstrap.sh\' && crab submit -c %s' % cfg)

print('Config files are stored in grid/')
