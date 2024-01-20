import ROOT
import sys
import numpy as np
from optimizeLeptonIsolation import canvasHeader

def getEfficiency(url,flavList,ptRange):

    t=ROOT.TChain('tree')
    t.AddFile(url)

    baseHistos={'csv'     : ROOT.TH1F('csv',';CSV;Jets',50,0,1),
                'csvdisc' : ROOT.TH1F('csv',';CSV;Jets',20,0,1),
                'pt'      : ROOT.TH1F('pt',';Transverse momentum [GeV];Jets',20,30,250),
                'eta'     : ROOT.TH1F('eta',';Pseudo-rapidity;Jets',20,0,2.0),
                'cent'    : ROOT.TH1F('centrality',';Centrality;Jets',10,0,100) }
    histos={}
    for key in baseHistos:
        for pfix in ['den','num','numtight']:
            histos[key+'_'+pfix]=baseHistos[key].Clone(key+'_'+pfix)
            histos[key+'_'+pfix].SetDirectory(0)
            histos[key+'_'+pfix].Sumw2()

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        ncoll=t.ncollWgt        
        for j in range(t.nbjet):

            genpt=t.bjet_genpt[j]
            geneta=t.bjet_geneta[j]
            if(abs(geneta)>2.0) : continue
            if genpt<ptRange[0] : continue
            if genpt>ptRange[1] : continue

            bflav=t.bjet_flavorB[j]
            if not abs(bflav) in flavList : continue
            csvVal=min(1.,max(0.,t.bjet_csvv2[j]))            

            #for efficiency extraction
            xbin=histos['csv_num'].GetXaxis().FindBin(csvVal)
            for ix in range(histos['csv_num'].GetNbinsX()+1):
                xcen=histos['csv_num'].GetXaxis().GetBinCenter(ix+1)
                if ix<xbin:
                    histos['csv_num'].Fill(xcen,ncoll)
                histos['csv_den'].Fill(xcen,ncoll)
                
            histos['csvdisc_den'].Fill(csvVal,ncoll)
            histos['pt_den'].Fill(genpt,ncoll)
            histos['eta_den'].Fill(abs(geneta),ncoll)
            histos['cent_den'].Fill(t.cenbin,ncoll)
            if csvVal>0.81:
                histos['csvdisc_num'].Fill(csvVal,ncoll)
                histos['pt_num'].Fill(genpt,ncoll)
                histos['eta_num'].Fill(abs(geneta),ncoll)
                histos['cent_num'].Fill(t.cenbin,ncoll)
            if csvVal>0.91:
                histos['csvdisc_numtight'].Fill(csvVal,ncoll)
                histos['pt_numtight'].Fill(genpt,ncoll)
                histos['eta_numtight'].Fill(abs(geneta),ncoll)
                histos['cent_numtight'].Fill(t.cenbin,ncoll)


    gr_eff=ROOT.TGraphAsymmErrors()
    gr_eff.Divide(histos['csv_num'],histos['csv_den'])    
    return gr_eff,histos


def showEfficiencyCurves(grColl,name):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    c.SetGridy()
    c.SetLogy()
    mg=ROOT.TMultiGraph()
    for g in grColl: mg.Add(g,'p')
    mg.Draw('ap')
    mg.GetXaxis().SetTitle('CSV')
    mg.GetYaxis().SetTitle('Efficiency')
    mg.GetXaxis().SetRangeUser(0,1)
    mg.GetYaxis().SetRangeUser(1e-4,1)
    leg=c.BuildLegend(0.15,0.15,0.4,0.3)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    l=ROOT.TLine()
    l.SetLineStyle(0)
    l.SetLineColor(ROOT.kGray+1)
    for y in [0.001,0.01,0.1]:
        l.DrawLine(0,y,1,y)
    canvasHeader()
    for ext in ['png','pdf']:
        c.SaveAs('%s.%s'%(name,ext))


def showHistos(histos,flav):

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    c.SetGridy()
    c.SetLogy()
    for key in histos:
        if '_num' in key: continue

        frame=histos[key].Clone('frame')
        frame.Reset('ICE')
        frame.GetYaxis().SetTitle('Efficiency or PDF')
        frame.GetYaxis().SetRangeUser(1e-3,1)
        frame.Draw()

        gr_eff=ROOT.TGraphAsymmErrors()
        gr_eff.SetMarkerStyle(20)
        gr_eff.Divide(histos[key.replace('_den','_num')],histos[key])

        gr_efftight=ROOT.TGraphAsymmErrors()
        gr_efftight.SetMarkerStyle(24)
        gr_efftight.SetMarkerColor(ROOT.kGray+1)
        gr_efftight.SetLineColor(ROOT.kGray+1)
        gr_efftight.Divide(histos[key.replace('_den','_numtight')],histos[key])

        histos[key].Scale(1./histos[key].Integral())        
        histos[key].Draw('histsame')
        histos[key].SetFillStyle(1001)
        histos[key].SetFillColor(ROOT.kGray)
        histos[key].SetLineColor(1)

        gr_eff.Draw('p')
        gr_efftight.Draw('p')
        
        if flav!='b':
            leg=ROOT.TLegend(0.65,0.93,0.95,0.75)
        else:
            leg=ROOT.TLegend(0.65,0.75,0.95,0.5)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.AddEntry(gr_eff,'Loose w.p.','ep')
        leg.AddEntry(gr_efftight,'Tight w.p.','ep')
        leg.AddEntry(histos[key],'Distribution','l')
        leg.Draw()

        canvasHeader(extraTxt=[flav])
        c.RedrawAxis()
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('%s_%s.%s'%(key,flav,ext))
        frame.Delete()
        gr_eff.Delete()


def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    baseDir=sys.argv[1]
    plottag=''
    if len(sys.argv)>2:
        plottag=sys.argv[2]
    mixSig='TTJets_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8.root'
    ppSample='/eos/cms/store/cmst3/group/hintt/PbPb2018_skim27Apr/TT_TuneCP5_5p02TeV-powheg-pythia8.root'
    ptRange=[30,500]

    #get the efficiency curves
    csv={}
    for tag,url,flavList,ms,ci in [ ('b',         baseDir+mixSig, [5],         20, 1),
                                    ('udsg',      baseDir+mixSig, [1,2,3,21],  24, 1),
                                    ('unmatched', baseDir+mixSig, [0],         21, 1),
                                ]:
        
        csv[tag],histos=getEfficiency(url,flavList,ptRange)
        csv[tag].SetTitle(tag)
        csv[tag].SetMarkerStyle(ms)
        csv[tag].SetLineColor(ci)
        csv[tag].SetMarkerColor(ci)

        showHistos(histos,tag)


    #tune the working point
    wpEff=0.05
    wpEff=0.01
    bestCut=0.0
    bestEff=csv['udsg'].Eval(bestCut)
    for x in np.arange(bestCut,1,0.01):
        eff=max(csv['udsg'].Eval(x),csv['unmatched'].Eval(x))
        if abs(eff-wpEff)>abs(bestEff-wpEff): continue
        bestCut=x
        bestEff=eff
    print('<'*50)
    print('csv>',bestCut)
    for tag in csv:
        print(tag,'eff=',csv[tag].Eval(bestCut))
    print('<'*50)

    #compare the curves
    for name,grNames in [ ('csveff',       ['b',  'udsg',  'unmatched']), 
                          #('csveff',       ['b',  'b (pp)', 'udsg',  'udsg (pp)', 'unmatched', 'unmatched (pp)']), 
                          ('beff',         ['b',    'b (0-30)',    'b (30-100)']),
                          ('udsgeff',      ['udsg', 'udsg (0-30)', 'udsg (30-100)']),
                          ('unmatchedeff', ['unmatched', 'unmatched (0-30)', 'unmatched (30-100)']),
                          ]:
        showEfficiencyCurves(grColl=[csv[x].Clone() for x in grNames],name=name+plottag)


if __name__ == "__main__":
    main()
