import ROOT
import sys
import os
import numpy as np
import pickle
from array import array
from HeavyIonsAnalysis.topskim.EventReader import *
from prepareCombinatorialBackgroundTree import prepareDileptonCollection

ISOTITLES=['I_{raw}','I(R=0.3)','I(R=0.2)','Mini isolation']
REQ_EFFB=0.4
MCTAG='DYJetsToLL_MLL-50_TuneCP5_HydjetDrumMB_5p02TeV-amcatnloFXFX-pythia8' #'TTJets_TuneCP5_HydjetDrumMB-amcatnloFXFX'

def canvasHeader(extraTxt=[]):
    txt=ROOT.TLatex()
    txt.SetNDC(True)
    txt.SetTextFont(42)
    txt.SetTextSize(0.04)
    txt.SetTextAlign(12)
    txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
    txt.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
    txt.DrawLatex(0.95,0.97,'#scale[0.8]{1752 #mub^{-1} (#sqrt{s_{NN}}=5.02 TeV)}')
    txt.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
    txt.SetTextSize(0.035)
    for i in range(len(extraTxt)):
        txt.DrawLatex(0.15,0.9-i*0.04,extraTxt[i])


def drawIsolationProfile(hmain,hextra=[],extraTxt=[],name='isoprof',doLegend=False,printBinContent=False,logX=False):

    """ draws the isolation profile """
   
    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.1)
    c.SetLogx(logX)
    if hmain.InheritsFrom('TH2'):
        c.SetRightMargin(0.12)
        hmain.Draw('colztexte' if printBinContent else 'colz')
    elif hmain.InheritsFrom('TH1'):
        c.SetRightMargin(0.03)
        hmain.Draw('hist')
    for h,drawOpt in hextra:
        h.Draw(drawOpt)
    if doLegend:
        leg=c.BuildLegend(0.6,0.93,0.95,0.93-0.05*(1+len(hextra)))
        leg.SetTextFont(42)
        leg.SetTextSize(0.03)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
    canvasHeader(extraTxt)
    c.Modified()
    c.Update()
    c.RedrawAxis()
    for ext in ['png','pdf']: #,'root']:
        c.SaveAs('{0}.{1}'.format(name,ext))


def drawROCSummary(grColl,name):

    colors=[ROOT.kBlack, ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kMagenta-9,ROOT.kRed+1,ROOT.kAzure+7, ROOT.kBlue-7]
    ncols=len(colors)
    lines=[1,2,3]

    c=ROOT.TCanvas('c','c',500,500)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.03)
    c.SetBottomMargin(0.1)
    c.SetGridy()

    mg=ROOT.TMultiGraph()
    for i in range(len(grColl)):
        grColl[i].SetLineColor( colors[i%ncols] )
        grColl[i].SetLineStyle( lines[i/ncols] )
        grColl[i].SetLineWidth(2)
        mg.Add(grColl[i],'l')
    mg.Draw('al')
    mg.GetXaxis().SetTitle('Signal efficiency')
    mg.GetYaxis().SetTitle('Background efficiency')
    mg.GetXaxis().SetRangeUser(0.5,1)
    mg.GetYaxis().SetRangeUser(0.2,1)
    leg=c.BuildLegend(0.15,0.94,0.4,0.94-0.03*len(grColl))
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)   
    l=ROOT.TLine(0.5,REQ_EFFB,1,REQ_EFFB)
    l.SetLineColor(ROOT.kGray+2)
    l.Draw()    
    canvasHeader()
    c.Modified()
    c.Update()
    c.RedrawAxis()
    for ext in ['png','pdf']:
        c.SaveAs("%s.%s"%(name,ext))


def getVariables(dilColl):

    """ prune the variables needed for the optimization of the isolation """
    kin=[]
    isoEstimators=[]
    globalEvent=[]
    for ll in dilColl:
        zWindow=True if abs(ll.p4.M()-91)<15 else False
        for i in range(2):
            l=getattr(ll,'l%d'%(i+1))
            kin.append( [l.pt,l.eta,l.phi] )

            uncorIso=l.chiso+l.nhiso+l.phiso
            if abs(l.pdgId)==13 : uncorIso=l.isofull30

            isoEstimators.append( [uncorIso,uncorIso,l.isofull20,l.miniiso*l.pt] )
            globalEvent.append([l.rho ,l.cenbin,l.ncollWgt,zWindow])

    return kin,isoEstimators,globalEvent


def getROC(sig, bkg, cutAtEffB):

    """build the efficiency and ROC curves"""

    #start the graphs
    roc_gr=ROOT.TGraph()
    roc_gr.SetName('roc')
    effs_gr=roc_gr.Clone('effs')
    effb_gr=roc_gr.Clone('effb')

    #loop over bins
    tot_sig=sig.Integral()
    tot_bkg=bkg.Integral()
    nbins=sig.GetNbinsX()
    bestEff,bestCut=None,None
    for xbin in range(nbins+1):
        cut=sig.GetXaxis().GetBinUpEdge(xbin)
        
        effsig=sig.Integral(0,xbin+1)/tot_sig if tot_sig>0 else 0
        effbkg=bkg.Integral(0,xbin+1)/tot_bkg if tot_bkg>0 else 0

        #if not bestEff or abs(effsig-cutAtEffS)<abs(bestEff-cutAtEffS):
        #    bestEff=effsig
        #    bestCut=cut

        if not bestEff or abs(effbkg-cutAtEffB)<abs(bestEff-cutAtEffB):
            bestEff=effbkg
            bestCut=cut

        roc_gr.SetPoint(xbin-1,effsig,effbkg)
        effs_gr.SetPoint(xbin-1,cut,effsig)
        effb_gr.SetPoint(xbin-1,cut,effbkg)

    return roc_gr,effs_gr,effb_gr, bestEff, bestCut

def tuneIsolation(mixFile,ch,matchedll=None):
    
    #readout Z and same-sign candidate events
    with open(mixFile,'r') as cache:
        allDileptons=pickle.load(cache)
        zll=allDileptons[(ch,True)]
        ssll=[ ll for ll in allDileptons[(ch,False)] if ll.isSS ] 

    kin,isoEstimators,globalEvent=getVariables(zll)
    ss_kin,ss_isoEstimators,ss_globalEvent=getVariables(ssll)
    mc_kin,mc_isoEstimators,mc_globalEvent=None,None,None
    if matchedll:
        mc_kin,mc_isoEstimators,mc_globalEvent=getVariables(matchedll)

    qkin=np.percentile(kin,list(range(0,120,30)),axis=0)
    qiso=np.percentile(isoEstimators,list(range(0,110,10)),axis=0)
    qGlobalEvent=np.percentile(globalEvent,list(range(0,110,10)),axis=0)

    rocCurves=[]
    sfisoHistos=[]
    for i in range(len(ISOTITLES)):

        #parametrize the profile versus rho

        #pol2
        #isoProf=ROOT.TF1('isoprof','[0]*pow(x+[1],2)+[2]*(x+[1])',0,1000)

        #piecewise
        isoFormula  = 'x<[0] ? '
        isoFormula += '[1]*pow(x,2)+[2]*x+[3]: '
        isoFormula += '[4]*pow(log(x),2)+[5]*log(x)'
        isoFormula += '+[3]+[1]*pow([0],2)-[4]*pow(log([0]),2)+[2]*[0]-[5]*log([0])'
        isoProf=ROOT.TF1('isoprof',isoFormula,0,1000)
        isoProf.SetParLimits(0,50,120)
        isoProf.SetParLimits(1,0,100)
        isoProf.SetParLimits(3,0,30) #pp typical range
        isoProf.SetParLimits(4,0,100)

        isovsrho=ROOT.TH2F('iso%dvsrho'%i, ';#rho;%s [GeV]'%ISOTITLES[i], 10,array('d',qGlobalEvent[:,0]),10,array('d',qiso[:,i]) )
        for k in range(len(isoEstimators)): 
            isovsrho.Fill(globalEvent[k][0],isoEstimators[k][i])
        
        #subtract same-sign contribution in the Z-window
        ss_isovsrho=isovsrho.Clone('ssiso%dvsrho'%i)
        ss_isovsrho.Reset('ICE')
        for k in range(len(ss_isoEstimators)): 
            if ss_globalEvent[k][3]: 
                ss_isovsrho.Fill(ss_globalEvent[k][0],ss_isoEstimators[k][i])
        isovsrho.Add(ss_isovsrho,-1)
        ss_isovsrho.Delete()

        #profile and fit
        isovsrho_prof=isovsrho.ProfileX()
        isovsrho_prof.Fit(isoProf,'MR+')
        isovsrho_prof.SetMarkerStyle(20)
        drawIsolationProfile(hmain=isovsrho,
                             hextra=[(isovsrho_prof,'e1same')],
                             extraTxt=['#bf{%s}'%ISOTITLES[i],
                                       'UE=#rho<#rho_{0} ? a#rho^{2}+b#rho+c : dlog(#rho)^{2}+elog(#rho)+f', 
                                       '#rho_{0}=%f'%isoProf.GetParameter(0),
                                       'a=%f'%isoProf.GetParameter(1),
                                       'b=%f'%isoProf.GetParameter(2),
                                       'c=%f'%isoProf.GetParameter(3),
                                       'd=%f'%isoProf.GetParameter(4),
                                       'e=%f'%isoProf.GetParameter(5)],
                             name='iso%dprof_%d'%(i,ch))


        #correct the isolation estimator
        coriso=ROOT.TH1F('coriso%d'%i, ';[%s-UE(#rho)]/p_{T};'%ISOTITLES[i],50,-3,4)
        corisovscen=ROOT.TH2F('coriso%dvscen'%i, ';Centrality (%%);[%s-UE(#rho)]/p_{T}'%ISOTITLES[i], 10,array('d',qGlobalEvent[:,1]),50,-3,4)
        corIsoEstimators=[]
        for k in range(len(isoEstimators)):            
            if i>0:
                corIso=(isoEstimators[k][i]-isoProf.Eval(globalEvent[k][0]))/kin[k][0]
            else:
                corIso=(isoEstimators[k][i])/kin[k][0]
            corIsoEstimators.append(corIso)
            corisovscen.Fill(globalEvent[k][1],corIso)
            coriso.Fill(corIso)

        #repeat for the same-sign events and subtract from the previous
        ss_corisovscen={}
        ss_coriso={}
        for tag in ['z','inc']:
            ss_corisovscen[tag]=corisovscen.Clone('%ssscoriso%dvscen'%(tag,i))
            ss_corisovscen[tag].Reset('ICE')
            ss_coriso[tag]=coriso.Clone('%ssscoriso%d'%(tag,i))
            ss_coriso[tag].Reset('ICE')
        ss_corIsoEstimators=[]
        for k in range(len(ss_isoEstimators)):
            if i>0:
                corIso=(ss_isoEstimators[k][i]-isoProf.Eval(ss_globalEvent[k][0]))/ss_kin[k][0]
            else:
                corIso=(ss_isoEstimators[k][i])/ss_kin[k][0]
            tag='z' if ss_globalEvent[k][3] else 'inc'
            ss_corIsoEstimators.append(corIso)
            ss_corisovscen[tag].Fill(ss_globalEvent[k][1],corIso)
            ss_coriso[tag].Fill(corIso)            
        corisovscen.Add(ss_corisovscen['z'],-1)
        coriso.Add(ss_coriso['z'],-1)        
        coriso.Scale(1./coriso.Integral())
        coriso.SetLineWidth(2)
        coriso.SetLineColor(1)
        coriso.GetYaxis().SetTitle('PDF')
        
        ss_coriso['inc'].Scale(1./ss_coriso['inc'].Integral())
        ss_coriso['inc'].SetLineColor(ROOT.kGray+2)
        ss_coriso['inc'].SetFillStyle(3001)
        ss_coriso['inc'].SetFillColor(ROOT.kGray)
            
        #repeat the same for MC truth if available
        mc_coriso=None
        mc_corIsoEstimators=[]
        if mc_isoEstimators:
            mc_coriso=coriso.Clone('mccoriso%d'%i)
            mc_coriso.Reset('ICE')
            mc_coriso.SetLineWidth(2)
            mc_coriso.SetLineColor(ROOT.kCyan+1)
            mc_coriso.SetLineStyle(1)
            for k in range(len(mc_isoEstimators)):
                if i>0:
                    corIso=(mc_isoEstimators[k][i]-isoProf.Eval(mc_globalEvent[k][0]))/mc_kin[k][0]
                else:
                    corIso=(mc_isoEstimators[k][i])/mc_kin[k][0]
                mc_corIsoEstimators.append(corIso)
                ncoll=mc_globalEvent[k][2]
                mc_coriso.Fill(corIso,ncoll)
            mc_coriso.Scale(1./mc_coriso.Integral())


        #ROC curve
        result=getROC(coriso,ss_coriso['inc'],REQ_EFFB)
        rocCurves.append(result[0])
        isoTitle='#bf{[%s-UE(#rho)]/p_{T}}'%ISOTITLES[i]
        rocCurves[-1].SetTitle(isoTitle)
        result[1].SetLineColor(1)
        result[2].SetLineColor(ROOT.kRed)
        cut=result[-1]

        #compare distributions
        corisovscen_prof=corisovscen.ProfileX()
        corisovscen_prof.SetMarkerStyle(20)
        drawIsolationProfile(hmain=corisovscen,
                             hextra=[(corisovscen_prof,'e1same')],
                             extraTxt=['#bf{[%s-UE(#rho)]/p_{T}}'%ISOTITLES[i]],
                             name='coriso%dprof_%d'%(i,ch))

        coriso.SetTitle('Z#rightarrow ll data')        
        ss_coriso['inc'].SetTitle('SS data')
        hextra=[(ss_coriso['inc'],'histsame')]
        if mc_coriso: 
            mc_coriso.SetTitle('Z#rightarrow ll MC')
            hextra.append( (mc_coriso,'histsame') )
        drawIsolationProfile(hmain=coriso,
                             hextra=hextra,
                             extraTxt=[isoTitle,'I_{rel}<%f'%cut],
                             name='coriso%d_%d'%(i,ch),
                             doLegend=True)
        

        #data/MC scale factors
        if mc_coriso:
            for tag in ['cen','periph','inc']:
                ptBins=[20,40,100]
                etaBins=[0,1.4442,1.5,2.5] #EE-EB
                if ch==169: 
                    etaBins=[0,0.8,2.5] #DT-RPC+CSC
                isoeff=ROOT.TH2F('effvsptvseta_iso%d'%i, 
                                 ';Transverse momentum [GeV];Pseudo-rapidity;Efficiency', 
                                 len(ptBins)-1,array('d',ptBins),len(etaBins)-1,array('d',etaBins))
                isoeff.Sumw2()
                isoeff_den=isoeff.Clone('isoeffden')               
                for k in range(len(kin)):
                    pt,eta,_=kin[k]
                    cenbin=globalEvent[k][1]
                    if tag=='cen' and cenbin>30 : continue
                    elif tag=='periph' and cenbin<=30 : continue
                    isoeff_den.Fill(pt,abs(eta))
                    if corIsoEstimators[k]>cut : continue
                    isoeff.Fill(min(pt,ptBins[-1]-0.1),min(abs(eta),etaBins[-1]-0.1))
                isoeff.Divide(isoeff_den)
                isoeff_den.Delete()

                mcisoeff=isoeff.Clone('mceffvsptvseta_iso%d'%i)
                mcisoeff.Reset('ICE')
                mcisoeff_den=mcisoeff.Clone('mcisoeffden')
                for k in range(len(mc_kin)):
                    pt,eta,_=mc_kin[k]
                    cenbin=mc_globalEvent[k][1]
                    if tag=='cen' and cenbin>30 : continue
                    elif tag=='periph' and cenbin<30 : continue
                    ncoll=mc_globalEvent[k][2]
                    ncoll=1
                    mcisoeff_den.Fill(pt,abs(eta),ncoll)
                    if mc_corIsoEstimators[k]>cut : continue
                    mcisoeff.Fill(pt,abs(eta),ncoll)

                mcisoeff.Divide(mcisoeff_den)
                mcisoeff_den.Delete()
                sfisoeff=isoeff.Clone('sfiso')
                sfisoeff.Divide(mcisoeff)
                mcisoeff.Delete()

                isoeff.GetXaxis().SetMoreLogLabels()
                isoeff.GetZaxis().SetRangeUser(0,1)
                drawIsolationProfile(hmain=isoeff,
                                     hextra=[],
                                     extraTxt=[isoTitle,'I_{rel}<%f'%cut],                                       
                                     name='%siso%deff_%d'%(tag,i,ch),
                                     printBinContent=True,
                                     logX=True
                                     )

                sfisoeff.GetXaxis().SetMoreLogLabels()
                sfisoeff.GetZaxis().SetTitle('SF = #varepsilon(data)/#varepsilon(MC)')
                sfisoeff.GetZaxis().SetRangeUser(0.8,1.2)
                #save a copy before drawing
                sfisoHistos.append( sfisoeff.Clone('sfiso%deff_%s_%d'%(i,tag,ch) ) )
                sfisoHistos[-1].SetDirectory(0)
                drawIsolationProfile(hmain=sfisoeff,
                                     hextra=[],
                                     extraTxt=[isoTitle,'I_{rel}<%f'%cut],                                       
                                     name='sfiso%deff_%s_%d'%(i,tag,ch),
                                     printBinContent=True,
                                     logX=True)

        continue
            
        #look for hotspots in same-sign data
        ssHotSpots=ROOT.TH2F('ssetavsphi',';Pseudo-rapidity;Azimuthal angle [rad];Candidates',25,-2.5,2.5,25,-3.15,3.15)
        isossHotSpots=ssHotSpots.Clone('isossetavsphi')
        if i==5:
            sshs_iso=ROOT.TH1F('sshsiso%d'%i,';%s;Leptons'%ISOTITLES[i],10,array('d',qiso[:,i]))
        else:
            sshs_iso=ROOT.TH1F('sshsiso%d'%i,';%s;Leptons'%ISOTITLES[i],50,-3,4)
        ssnhs_iso=sshs_iso.Clone('ssnhsiso%d'%i)
        ziso=sshs_iso.Clone('ziso%d'%i)
        for k in range(len(ss_kin)):
            
            pt,eta,phi=ss_kin[k]
            HEM1516failure=True if eta<-1.39 and -1.6<phi<-0.9 else False

            isoVal=ss_isoEstimators[k][i] if i==5 else ss_corIsoEstimators[k]
            if HEM1516failure:
                sshs_iso.Fill(isoVal)
            else:
                ssnhs_iso.Fill(isoVal)

            ssHotSpots.Fill(eta,phi)
            if ss_corIsoEstimators[k]>cut : continue
            isossHotSpots.Fill(eta,phi)

        for k in range(len(kin)):
            isoVal=isoEstimators[k][i] if i==5 else corIsoEstimators[k]
            ziso.Fill(isoVal)

        if i==0: drawIsolationProfile(hmain=ssHotSpots,hextra=[],extraTxt=[],name='sshotspots_%d'%ch)
        drawIsolationProfile(hmain=isossHotSpots,hextra=[],extraTxt=[isoTitle,'I_{rel}<%f'%cut],name='iso%dsshotspots_%d'%(i,ch))
        ssHotSpots.Delete()
        isossHotSpots.Delete()

        if ch==121:
            for h in [ssnhs_iso,sshs_iso,ziso]:
                h.Scale(1./h.Integral())
                h.SetLineWidth(2)
                h.GetYaxis().SetTitle('PDF')
                if i!=5: continue
                for xbin in range(ssnhs_iso.GetNbinsX()):
                    wid=ssnhs_iso.GetXaxis().GetBinWidth(xbin+1)                
                    h.SetBinContent(xbin+1,h.GetBinContent(xbin+1)/wid)
                    h.SetBinError(xbin+1,h.GetBinError(xbin+1)/wid)
            ziso.SetLineColor(ROOT.kRed)
            sshs_iso.SetLineColor(ROOT.kGray+2)
            sshs_iso.SetFillStyle(3001)
            sshs_iso.SetFillColor(ROOT.kGray)
            drawIsolationProfile(hmain=ziso,
                                 hextra=[(sshs_iso,'histsame'),(ssnhs_iso,'histsame')],
                                 extraTxt=[ISOTITLES[i] if i==5 else isoTitle],
                                 name='iso%d_hotspotcomp_%d'%(i,ch))

    #overall performance summary
    drawROCSummary(rocCurves,'rocsummary_%d'%ch)
    return sfisoHistos

def main():
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kLightTemperature)        
    ROOT.gStyle.SetPaintTextFormat('4.3f')

    matchedll={11*11:None,13*13:None}
    if len(sys.argv)>1:
       
        print('Processing MC truth from',MCTAG)
        prepareDileptonCollection(sys.argv[1],MCTAG,maxEvents=100000)

        #readout matched leptons
        print('Opening pickle file')
        with open('dilepton_summary_%s.pck'%MCTAG,'r') as cache:
            allDileptons=pickle.load(cache)
            for ch in matchedll:
                #matchedll[ch]=[ll for ll in allDileptons[(ch,False)] if ll.l1.matched and ll.l2.matched]
                matchedll[ch]=[ll for ll in allDileptons[(ch,True)] ]
                print(len(matchedll[ch]),ch)

    sfisoHistos=[]
    for ch in matchedll:
        sfisoHistos += tuneIsolation('dilepton_summary.pck',ch,matchedll[ch])


    fOut=ROOT.TFile.Open('data/isolation_sf.root','RECREATE')
    for h in sfisoHistos:
        h.SetDirectory(fOut)
        h.Write()
    fOut.Close()


if __name__ == "__main__":
    main()
