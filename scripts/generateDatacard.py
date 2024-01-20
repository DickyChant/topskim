import ROOT
import sys
import os

def generateDatacardFor(cat):

    with open('datacard_%s.dat'%cat,'w') as dc:
        dc.write('imax *\n')
        dc.write('jmax *\n')
        dc.write('kmax *\n')
        dc.write('---------------\n')
        dc.write('shapes * * plotter_count.root {0}_acopl/acopl_$PROCESS {0}_acopl/acopl_$PROCESS_$SYSTEMATIC\n'.format(cat))
        dc.write('shapes data_obs * plotter_count.root {0}_acopl/acopl_Data\n'.format(cat))
        dc.write('---------------\n')
        dc.write('bin {0}\n'.format(cat))
        dc.write('observation -1\n')
        dc.write('------------------------------\n')
        dc.write('bin             {0}        {0}   {0}  {0}              {0}\n'.format(cat))
        dc.write('process         tbart_tbart tW_tW  VV_VV Combdata_Combdata Zgamma_Zgamma\n')
        dc.write('process         0           1     2       3                4\n')
        dc.write('rate            -1          -1     -1    -1                -1\n')
        dc.write('------------------------------\n')
        dc.write('lumi     lnN    1.05        1.05   1.05  -                -\n')
        dc.write('combnorm{0} lnN    -           -      -     1.3               -\n'.format(cat))
        dc.write('dynorm{0}   lnN    -           -      -     -                 1.3\n'.format(cat))
        dc.write('twnorm   lnN    -           1.3      -     -               -\n')
        dc.write('vvnorm   lnN    -           -      1.3     -               -\n')
        dc.write('{0} autoMCStats 0.0 1\n'.format(cat)) 

def getFitresult(url):

    res=[]
    fIn=ROOT.TFile.Open(url)
    t=fIn.Get('limit')
    for i in range(3):
        t.GetEntry(i)
        res.append( t.r )
        if i==0 : continue
        res[-1]=res[-1]-res[0]
    fIn.Close()
    return res

fit=sys.argv[1]
cats=[]
if fit=='inc':
    cats=['em','ee','mm']
elif fit=='inchpur':
    cats=['emhpur','eehpur','mmhpur']
elif fit=='incbtag':
    cats=['em0pfb','emgeq1pfb','mm0pfb','mmgeq1pfb','ee0pfb','eegeq1pfb']
elif fit=='hpurbtag':
    cats=['emhpur0pfb','emhpurgeq1pfb','mmhpur0pfb','mmhpurgeq1pfb','eehpur0pfb','eehpurgeq1pfb']
else:
    print('Unknown fit:',fit)
    sys.exit(-1)

opts=''
for c in cats: 
    generateDatacardFor(c)
    opts='{0} {1}=datacard_{1}.dat'.format(opts,c)

os.system('combineCards.py {0} > finaldatacard.dat'.format(opts))
os.system('combine -M MultiDimFit finaldatacard.dat -t -1 --expectSignal=1 --saveFitResult --robustFit=1 --algo=cross --cl=0.68')
resTot=getFitresult('higgsCombineTest.MultiDimFit.mH120.root')
#os.system('combine -M MultiDimFit finaldatacard.dat -t -1 --expectSignal=1 --saveFitResult --robustFit=1 --algo=cross --cl=0.68 -S0')
resStat=getFitresult('higgsCombineTest.MultiDimFit.mH120.root')

systHi=ROOT.TMath.Sqrt(resTot[1]**2-resStat[1]**2)
systLo=ROOT.TMath.Sqrt(resTot[2]**2-resStat[2]**2)
print('%3.3f +%3.3f-%3.3f (syst) +%3.3f-%3.3f (stat)'%(resTot[0],systHi,systLo,resStat[1],resStat[2]))
print(resTot)
print(resStat)
