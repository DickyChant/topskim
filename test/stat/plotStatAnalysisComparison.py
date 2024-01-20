import ROOT
import os
import sys
import optparse
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.ndimage import filters
import numpy as np

from statAnalysisTools import *
from plotPostFitDistributions import *

basedir = '/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/'

def main():

    fitres={
        ## 'sameFlavorsphericity'        : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix/sphericity/sameFlavor.card',
        ##                                  'title':'$S,ee,\mu\mu$'},
        ## 'sameFlavorsphericity25'      : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/sphericity/sameFlavor.card',
        ##                                  'title':'$S,ee,\mu\mu,25/25$'},
        ## 'sameFlavorsphericitybtags'   : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix_jetAnalysis/sphericity/sameFlavor.card',
        ##                                  'title':'$S,ee,\mu\mu+b$'},
        ## 'sameFlavorsphericitybtags25' : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/sphericity/sameFlavor.card',
        ##                                  'title':'$S,ee,\mu\mu+b,25/25$'},
        ## 'sameFlavorbdtcomb'           : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix/bdtcomb/sameFlavor.card',
        ##                                  'title':'BDT,$ee,\mu\mu$'},
        ## 'sameFlavorbdtcomb25'         : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/bdtcomb/sameFlavor.card',
        ##                                  'title':'$BDT,ee,\mu\mu,25/25$'},
        ## 'sameFlavorbdtcombbtags'      : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix_jetAnalysis/bdtcomb/sameFlavor.card',
        ##                                  'title':'$BDT,ee,\mu\mu+b$'},
        ## 'sameFlavorbdtcombbtags25'    : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/bdtcomb/sameFlavor.card',
        ##                                  'title':'$BDT,ee,\mu\mu+b,25/25$'},
        ## 'emsphericity'                : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix/sphericity/em.card',
        ##                                  'title':'$S,e\mu$'},
        ## 'emsphericity25'              : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/sphericity/em.card',
        ##                                  'title':'$S,e\mu,25/25$'},
        ## 'emsphericitybtags'           : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix_jetAnalysis/sphericity/em.card',
        ##                                  'title':'$S,e\mu+b$'},
        ## 'emsphericitybtags25'         : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/sphericity/em.card',
        ##                                  'title':'$S,e\mu+b,25/25$'},
        ## 'embdtcomb'                   : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix/bdtcomb/em.card',
        ##                                  'title':'$BDT,e\mu$'},
        ## 'embdtcomb25'                 : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/bdtcomb/em.card',
        ##                                  'title':'$BDT,e\mu,25/25$'},
        ## 'embdtcombbtags'              : {'dir':'datacards_2019-09-19_elPt25muPt20_ttbarBRFix_jetAnalysis/bdtcomb/em.card',
        ##                                  'title':'$BDT,e\mu+b$'},
        ## 'embdtcombbtags25'            : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/bdtcomb/em.card',
        ##                                  'title':'$BDT,e\mu+b,25/25$'},
        'allFlavorssphericity'        : {'dir':basedir+'datacards_2019-10-21_forApproval/sphericity/combinedCard/',
                                         'title':'$S,ee,\mu\mu,e\mu$'},
        ## 'allFlavorssphericity25'      : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/sphericity/allFlavors.card',
        ##                                  'title':'$S,ee,\mu\mu,e\mu,25/25$'},
        'allFlavorssphericitybtags'   : {'dir':basedir+'datacards_2019-10-21_forApproval_jetAnalysis/sphericity/combinedCard/',
                                         'title':'$S,ee,\mu\mu,e\mu+b$'},
        ## 'allFlavorssphericitybtags25' : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/sphericity/allFlavors.card',
        ##                                  'title':'$S,ee,\mu\mu,e\mu+b,25/25$'},
        'allFlavorsbdtcomb'           : {'dir':basedir+'datacards_2019-10-21_forApproval/bdtcomb/combinedCard/',
                                         'title':'$BDT,ee,\mu\mu,e\mu$'},
        ## 'allFlavorsbdtcomb25'         : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/bdtcomb/allFlavors.card',
        ##                                  'title':'$BDT,ee,\mu\mu,e\mu,25/25$'},
        'allFlavorsbdtcombbtags'      : {'dir':basedir+'datacards_2019-10-21_forApproval_jetAnalysis/bdtcomb/combinedCard/',
                                         'title':'$BDT,ee,\mu\mu,e\mu+b$'},
        ## 'allFlavorsbdtcombbtags25'    : {'dir':'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/bdtcomb/allFlavors.card',
        ##                                  'title':'$BDT,ee,\mu\mu,e\mu+b,25/25$'},
    }
    toCompare=[]
    for ch in ['allFlavors']: ##'sameFlavor','em'
        for d in ['sphericity','bdtcomb']:
            for a in ['','btags']: ##'25', 'btags25'
                toCompare.append(ch+d+a)

    """
    basedir = '/afs/cern.ch/work/m/mdunser/public/cmssw/heavyIons/CMSSW_9_4_6_patch1/src/CMGTools/TTHAnalysis/python/plotter/hin-ttbar/'
    fitres={
        'sphericity2520_bjet': {'dir': basedir+'datacards_2019-09-19_elPt25muPt20_ttbarBRFix_jetAnalysis/sphericity/combinedCard/', 'title': 'sphericity (25/20) + b'},
        'sphericity2520_incl': {'dir': basedir+'datacards_2019-09-19_elPt25muPt20_ttbarBRFix/sphericity/allFlavors/'              , 'title': 'sphericity (25/20)    '},
        'sphericity2525_bjet': {'dir': basedir+'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/sphericity/combinedCard/', 'title': 'sphericity (25/25) + b'},
        'sphericity2525_incl': {'dir': basedir+'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/sphericity/allFlavors/'              , 'title': 'sphericity (25/20)    '},
        'bdtcomb2520_bjet'   : {'dir': basedir+'datacards_2019-09-19_elPt25muPt20_ttbarBRFix_jetAnalysis/bdtcomb/combinedCard/'   , 'title': 'bdt (25/20) + b'},
        'bdtcomb2520_incl'   : {'dir': basedir+'datacards_2019-09-19_elPt25muPt20_ttbarBRFix/bdtcomb/allFlavors/'                 , 'title': 'bdt (25/20)    '},
        'bdtcomb2525_bjet'   : {'dir': basedir+'datacards_2019-09-20_elPt25muPt25_ttbarBRFix_jetAnalysis/bdtcomb/combinedCard/'   , 'title': 'bdt (25/25) + b'},
        'bdtcomb2525_incl'   : {'dir': basedir+'datacards_2019-09-20_elPt25muPt25_ttbarBRFix/bdtcomb/allFlavors/'                 , 'title': 'bdt (25/25)    '},
    }
    toCompare=[]
    for d in ['sphericity','bdtcomb']:
        for pt in ['2520','2525']:
            for bj in ['bjet', 'incl']:
                toCompare.append(d+pt+'_'+bj)
    """

    muPlots={'obs':ROOT.TGraphAsymmErrors(),'exp':ROOT.TGraphAsymmErrors()}
    sigPlots={'obs':ROOT.TGraphAsymmErrors(),'exp':ROOT.TGraphAsymmErrors()}
    frame=ROOT.TH2F('frame','frame',10,0,10,len(toCompare),0,len(toCompare))

    for ipt in range(len(toCompare)):

        r=toCompare[ipt]
        if not r in fitres : continue

        url=fitres[r]['dir']
        title=fitres[r]['title']
        print('Parsing results from',url)
                  
        fitres[r]['sigtoys'] = getResultsFromToys(url+'/higgsCombinesig_toys.Significance.mH120.123456.root','limit')
        sigobs=getResultsFromToys(url+'/higgsCombinesig_obs.Significance.mH120.root','limit')
        fitres[r]['sigobs']  = sigobs[0] 
        fitres[r]['rtoys']   = getResultsFromToys(url+'/higgsCombineexpsinglestoys.MultiDimFit.mH172.5.123456.root','r')
        for fit in ['exp','obs']:
            fitres[r]['r'+fit]=getNLLResults('%s/higgsCombine%sscan.MultiDimFit.mH172.5.root'%(url,fit),
                                               '%s/higgsCombine%ssingles.MultiDimFit.mH172.5.root'%(url,fit))
            if os.path.isfile(url+'/higgsCombine%sscan2d.MultiDimFit.mH172.5.root'%fit):
                fitres[r]['ll2d'+fit] = getResultsFromToys(url+'/higgsCombine%sscan2d.MultiDimFit.mH172.5.root'%fit,['r','muzg','deltaNLL'])


        #show postfit distributions
        doPostFitPlot(url+'/fitDiagnosticsobs.root')

        #show nll scan
        exp=fitres[r]['rexp']        
        obs=fitres[r]['robs']
        showNLLScan(exp=exp, obs=obs, name='%s_nllscan'%r,title=title)
        if 'll2dexp' in fitres[r]:
            showNLLContour(exp=fitres[r]['ll2dexp'], obs=fitres[r]['ll2dobs'],name='%s_nll2dscan'%r,title=title)

        #add to fit summary
        rtitle=title.replace('\\','#')
        rtitle=rtitle.replace('$','')
        frame.GetYaxis().SetBinLabel(ipt+1,rtitle)

        muPlots['exp'].SetPoint(ipt,exp[1][0],ipt+0.6)
        muPlots['exp'].SetPointError(ipt,abs(exp[1][1]),exp[1][2],0,0)
        muPlots['obs'].SetPoint(ipt,obs[1][0],ipt+0.4)
        muPlots['obs'].SetPointError(ipt,abs(obs[1][1]),obs[1][2],0,0)
        qtl=[np.percentile(fitres[r]['sigtoys'],pval) for pval in [16,50,84] ]
        sigPlots['exp'].SetPoint(ipt,qtl[0],ipt+0.6)
        sigPlots['exp'].SetPointError(ipt,qtl[1]-qtl[0],qtl[2]-qtl[1],0,0)
        sigPlots['obs'].SetPoint(ipt,fitres[r]['sigobs'],ipt+0.4)


        #compare toys with observations
        for var,xtitle,xbins,xran,yran in [ ('r',   '$\mu=\sigma/\sigma_{exp}$', 25, (0.25,1.7), (0,120)),
                                            ('sig', 'Significance',              50, (0.,10.),    (0,120))]:
            toys=fitres[r][var+'toys']
            obs=fitres[r][var+'obs'] 
            if var=='r' : obs=obs[1][0]
            showToys(toys=toys,
                     obs=obs,
                     name='%s_%s_toys'%(r,var),
                     title=title,
                     xtitle=xtitle,
                     xbins=xbins,
                     xran=xran,
                     yran=yran)

        #copy plots in sub-directories locally
        urlplots=[fname for fname in os.listdir(url) if '.png' in fname or '.pdf' in fname]
        for fname in urlplots:
            os.system('cp -v %s %s_%s'%(os.path.join(url,fname),r,fname))


    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    c=ROOT.TCanvas('c','c',500,500)
    c.SetLeftMargin(0.32)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetBottomMargin(0.1)
    c.SetGridy()
    for xvar in ['r','sig']:
        
        plots=muPlots if xvar=='r' else sigPlots

        frame.Draw()
        if xvar=='r' : frame.GetXaxis().SetRangeUser(0,2)
        else : frame.GetXaxis().SetRangeUser(0,8)
        frame.GetYaxis().SetLabelSize(0.05)
        frame.GetXaxis().SetTitleSize(0.05)
        frame.GetXaxis().SetTitleOffset(0.9)
        frame.GetXaxis().SetTitle('#mu=#sigma/#sigma_{th}' if xvar=='r' else 'Significance')

        plots['exp'].SetTitle('Exp.')
        plots['exp'].SetMarkerStyle(24)
        plots['exp'].SetLineColor(ROOT.kGray)
        plots['exp'].SetMarkerColor(ROOT.kGray)
        plots['exp'].Draw('p')

        plots['obs'].SetTitle('Obs.')
        plots['obs'].SetLineWidth(2)
        plots['obs'].SetMarkerStyle(20)
        plots['obs'].Draw('p')

        leg=ROOT.TLegend(0.75,0.95,0.98,0.82)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.045)
        leg.SetTextFont(42)
        leg.AddEntry(plots['obs'],'Obs.','lp')
        leg.AddEntry(plots['exp'],'Exp.','lp')
        leg.Draw()

        tex=ROOT.TLatex()
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.SetNDC()
        tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
        tex.DrawLatex(0.32,0.96,'#bf{CMS} #it{preliminary}')
        tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
        tex.DrawLatex(0.95,0.96,'1.7 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)')
    
        c.Modified()
        c.Update()
        for ext in ['png','pdf']:
            c.SaveAs('fitsummary_%s.%s'%(xvar,ext))

    #put all under StatAnaSummary
    os.system('mkdir -p StatAnaSummary')
    os.system('mv *.{png,pdf} StatAnaSummary')
    print('Plots are available in StatAnaSummary')

if __name__ == "__main__":
    sys.exit(main())
