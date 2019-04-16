import ROOT
import sys
import os
import pickle
import random
from HeavyIonsAnalysis.topskim.LeptonObject import *
from HeavyIonsAnalysis.topskim.HistoTool import *

url=sys.argv[1]
plotter=sys.argv[2] if len(sys.argv)>2 else None

#parse the data files
if not '.pck' in url:

    #build the chain
    t=ROOT.TChain('tree')
    for f in os.listdir(url):
        if not 'Skim' in f or not '.root' in f: continue
        t.Add(os.path.join(url,f))

    #loop over events
    print 'Analysing',t.GetEntries(),'events'
    dilCollection={11*11:[],11*13:[],13*13:[]}
    for iev in range(t.GetEntries()):
        t.GetEntry(iev)
        try:
            dil=getDilepton(t,[13,11])
            dilCollection[dil.flavour].append(dil)
        except Exception as e:
            print e
            pass

    #save dict in a cache
    url='dilepton_summary.pck'
    with open(url,'w') as cache:
        pickle.dump( dilCollection,cache,pickle.HIGHEST_PROTOCOL)

#load dilepton summary from pickle
with open(url,'r') as cache:
    dilCollection=pickle.load(cache)

#load previously plotted distributions
mixWeights={}
flavCats={13*13:'mm',11*11:'ee',13*11:'em'}
weightDist='ptll'
if plotter:
    fIn=ROOT.TFile.Open(plotter)    
    for flav in flavCats:
        cat=flavCats[flav]
        ratio=fIn.Get('{0}_{1}/{1}_{0}_SS'.format(cat,weightDist))
        ratio.Divide(fIn.Get('{0}_{1}/{1}_{0}mix_mix'.format(cat,weightDist)))
        mixWeights[flav]=ROOT.TGraph(ratio)
    fIn.Close()

ht=HistoTool()
ht.add( ROOT.TH1F("lpt",   ";Lepton transverse momentum [GeV];Leptons",10,20,200) )
ht.add( ROOT.TH1F('leta',  ";Lepton pseudo-rapidity;Leptons",10,0,2.5) )
ht.add( ROOT.TH1F('liso',  ";Relative isolation;Leptons",10,-1,4) )
ht.add( ROOT.TH1F('acopl', ";1-#Delta#phi(l,l')/#pi;Events",10,0,1.0) )
ht.add( ROOT.TH1F('ptll',  ";Dilepton transverse momentum [GeV];Events",10,0,200) )
ht.add( ROOT.TH1F('mll',   ";Dilepton invariant mass [GeV];Events",20,0,200) )

for key in dilCollection:
    print 'processing',flavCats[key]

    for dil in dilCollection[key]:
        
        if dil.isZ : continue

        toPlot=[]
        leptons=[dil.l1,dil.l2]

        #mix with another event of the same flavour
        for i in range(1,3):
            a=getattr(dil,'l%d'%i)
            while True:
                otherDil=random.choice(dilCollection[key])
                if dil.isZ : continue
                b=otherDil.l1 if otherDil.l1.pdgId==a.pdgId else otherDil.l2
                toPlot.append( (Dilepton(a,b,dil.isOF,dil.isSS,dil.isZ,dil.isIso,True),'mix', 1) )                
                if key in mixWeights:
                    obsVal=(a.p4+b.p4).M() if weightDist=='mll' else (a.p4+b.p4).Pt()
                    wgt=mixWeights[key].Eval(obsVal)
                    toPlot.append( (Dilepton(a,b,dil.isOF,dil.isSS,dil.isZ,dil.isIso,True),'mixwgt', wgt) )
                break
            
        if dil.isSS:
            toPlot.append( (dil,'',1) )

        #fill histograms
        for obj,tag,wgt in toPlot:
            cats=[flavCats[key]+tag]
            if obj.isIso : cats += [flavCats[key]+tag+'iso']
            ht.fill( (obj.l1.p4.Pt(),wgt),           'lpt',  cats, 'lead')
            ht.fill( (abs(obj.l1.p4.Eta()),wgt),     'leta', cats, 'lead')
            ht.fill( (obj.l1.getIsolation(),wgt),    'liso', cats, 'lead')
            ht.fill( (obj.l2.p4.Pt(),wgt),           'lpt',  cats, 'sublead')
            ht.fill( (abs(obj.l2.p4.Eta()),wgt),     'leta', cats, 'sublead')
            ht.fill( (obj.l2.getIsolation(),wgt),    'liso', cats, 'sublead')
            acopl=1-abs(leptons[0].p4.DeltaPhi(dil.l2.p4))/ROOT.TMath.Pi()
            ht.fill( (acopl,wgt),       'acopl', cats)
            ht.fill( (obj.p4.M(),wgt),  'mll',   cats)
            ht.fill( (obj.p4.Pt(),wgt), 'ptll' , cats)


url='combbackground_plots.root'
ht.writeToFile(url)
print 'Plots can be found in',url
