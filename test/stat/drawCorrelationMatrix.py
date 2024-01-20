import ROOT
import sys

ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPaintTextFormat('4.2f')

fIn=ROOT.TFile.Open(sys.argv[1]+'/fitDiagnosticsobs.root')
h=fIn.Get('fit_s').correlationHist()

c=ROOT.TCanvas('c','c',500,500)
c.SetTopMargin(0.05)
c.SetLeftMargin(0.12)
c.SetRightMargin(0.2)
c.SetBottomMargin(0.1)

h.Draw('colz')
h.GetZaxis().SetTitle('Post-fit correlation')

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.035)
tex.SetNDC()
tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
tex.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
tex.DrawLatex(0.95,0.97,'1.8 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)')

c.SaveAs('correlation.png')

corr=[]
for rbin in range(1,h.GetNbinsX()+1):
    name=h.GetXaxis().GetBinLabel(rbin)
    if name!='r': continue

    for ybin in range(1,h.GetNbinsY()+1):
        corr.append( (h.GetXaxis().GetBinLabel(rbin),
                      h.GetYaxis().GetBinLabel(ybin),
                      abs(h.GetBinContent(rbin,ybin)),
                      h.GetBinContent(rbin,ybin),
                      ybin) )
    break

corr.sort(key=lambda x: x[2],reverse=True) 
redcorr=[x for x in corr if x[2]>0.1]
#redcorr=[x for x in corr]
for a,b,_,corr,ybin in redcorr:
    print(a,b,corr,ybin)

#by value 
binsOfInterest=[x[4] for x in redcorr]

#forced selection
#forcedSelection=['r','lumi','tW_lnN','data_comb_lnN','ptZ','muzg',
#                 'quenching','btagB','btagUDSG']
#binsOfInterest=[x[4] for x in redcorr if x[1] in forcedSelection]



nred=len(binsOfInterest)
hred=ROOT.TH2F('corr','',nred,0,nred,nred,0,nred)
for iy in range(nred):

    ybin=binsOfInterest[iy]
    label=h.GetYaxis().GetBinLabel(ybin)
    hred.GetYaxis().SetBinLabel(iy+1,label)
    hred.SetBinContent(iy,iy,1.)
        
    for ix in range(nred):
        xbin=h.GetNbinsX()-(binsOfInterest[ix]-1)
        label=h.GetXaxis().GetBinLabel(xbin)
        hred.GetXaxis().SetBinLabel(ix+1,label)
        hred.SetBinContent(ix+1,iy+1,h.GetBinContent(xbin,ybin))
        hred.SetBinContent(iy+1,ix+1,h.GetBinContent(xbin,ybin))

hred.Draw('colz text')
hred.GetZaxis().SetTitle('Post-fit correlation')

tex=ROOT.TLatex()
tex.SetTextFont(42)
tex.SetTextSize(0.035)
tex.SetNDC()
tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
tex.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
tex.DrawLatex(0.95,0.97,'1.8 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)')

c.SaveAs('redcorrelation.png')
