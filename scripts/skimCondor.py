import os, optparse

## usage:
## python skimCondor.py -d /eos/cms/store/hidata/HIRun2018A/HIHardProbes/AOD/04Apr2019-v1/ -n 100 -o /eos/cms/store/cmst3/group/hintt/SkimAOD_HardProbes_Aug05/ --dryRun


if __name__ == '__main__':                                                                                                                                                                    

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('-d', '--directory', type=str, default='',    help='directory with the output directories [%default]')
    parser.add_option(      '--dryRun',              default=False, help='dry run (do not submit to condor) [%default]', action='store_true')
    parser.add_option('-n', '--nfiles',  type=int, default=100 ,  help='skim nfiles files in each job.')
    parser.add_option('-o', '--targetdir', type=str, default='',    help='target directory for the output [%default]')
    (options, args) = parser.parse_args()

    if not options.targetdir:
        print('no target directory given... exiting')
        sys.exit(0)
    if not options.directory:
        print('no source directory given... exiting')
        sys.exit(0)

    cmsenv = os.environ['CMSSW_BASE']

    filelist = []

    os.system('mkdir -p '+options.targetdir)

    t = 'SingleMuon' if 'SingleMuon' in options.directory else 'HardProbes' if 'HardProbes' in options.directory else ''

    ## get the list of valid files. the list should be readable for everyone
    validFileListBase = '/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_10_3_3_patch1/src/HeavyIonsAnalysis/topskim/scripts/'

    if t == 'SingleMuon':
        fn_validFileList = validFileListBase+'/validFiles_SingleMuon.txt'
    elif t == 'HardProbes':
        fn_validFileList = validFileListBase+'/validFiles_HardProbes.txt'

    f_validFileList = open(fn_validFileList, 'r')

    validFileList = [i.strip().split('/')[-1] for i in f_validFileList.readlines() ]

    ## done getting the valid file list

    nInvalid = 0

    for p,d,f in os.walk(options.directory) : #'/eos/cms/store/hidata/HIRun2018A/HISingleMuon/AOD/04Apr2019-v1/'):
        for fn in f:
            if not '.root' in fn: continue
            if not fn in validFileList:
                nInvalid += 1
                continue
            filelist.append(os.path.join(p,fn))

    print('there are {n} valid files in the directory {d} and its subdirectories'.format(n=len(filelist), d=options.directory))

    print('grouping them into bunches of', str(options.nfiles))


    bunches = []
    nbunches = 0
    tmp_str = ''
    for fi,f in enumerate(filelist):
        if fi and ( fi%options.nfiles == 0 or fi+1 == len(filelist) ) :
            bunches.append(tmp_str)
            tmp_str = ''
        else:
            tmp_str += '"file:{f}",\n'.format(f=f)

    #print bunches

    dummy = open('skim_dummy_cfg.py','r')
    dummylines = dummy.readlines()

    arguments = []

    os.system('mkdir -p skimCondor'+t)
    for ib,bunch in enumerate(bunches):
        nf = open('skimCondor{t}/skimCfg_{t}_cfg_{n}.py'.format(t=t,n=ib), 'w')
        tmpfn = 'skim'+t+'_'+str(ib)+'.root'
        for line in dummylines:
            nf.write(line.replace('XXXXX',bunch).replace('YYYYY',tmpfn))

        arguments.append('{c} {f} {fn} {tar}'.format(c=os.path.abspath('.'), f=os.path.abspath('.')+'/'+nf.name, fn=tmpfn, tar=options.targetdir))

        nf.close()
        

    tmp_condor_filename = 'condor_skim{t}.condor'.format(t=t)
    tmp_condor = open(tmp_condor_filename,'w')
    tmp_condor.write('''Executable = {here}/skimScript.sh
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = skimCondor{t}/skim_{t}_$(ProcId).log
Output     = skimCondor{t}/skim_{t}_$(ProcId).out
Error      = skimCondor{t}/skim_{t}_$(ProcId).error
getenv      = True

environment = "LS_SUBCWD={here}"
{ag}
+MaxRuntime = 12000\n\n'''.format(t=t, here=os.environ['PWD'],ag = '+AccountingGroup = "group_u_CMST3.all"' if os.environ['USER'] in ['mdunser', 'psilva'] else ''))

    for arg in arguments:
        tmp_condor.write('arguments = '+arg+'\n')
        tmp_condor.write('queue 1\n\n')

    tmp_condor.close()

    
    if options.dryRun:
        print('this was just a dry run...')
        print('wrote file', tmp_condor_filename)
     
    else:
        print('submitting to condor...')
        os.system('condor_submit '+ tmp_condor_filename)


