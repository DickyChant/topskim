import ROOT
import sys

baseDir=sys.argv[1]
sig=ROOT.TFile.Open('%s/TTJets_TuneCP5_HydjetDrumMB-amcatnloFXFX.root'%baseDir)
sigT=sig.Get('tree')
data=ROOT.TFile.Open('%s/SkimMuons_04Apr2019-v1.root'%baseDir)
dataT=data.Get('tree')

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

cnv=ROOT.TCanvas('c','c',500,500)
cnv.SetLeftMargin(0.12)
cnv.SetTopMargin(0.05)
cnv.SetBottomMargin(0.1)
cnv.SetRightMargin(0.03)
sigT.SetLineWidth(2)
sigT.SetLineColor(ROOT.kGray)
sigT.Draw('cenbin','','histnorm')
sigT.SetLineColor(1)
sigT.Draw('cenbin','ncoll','histnormsame')
dataT.SetMarkerStyle(20)
dataT.Draw('cenbin','','e1normsame')

mc=cnv.GetListOfPrimitives().At(0)
mc.GetXaxis().SetTitle('Centrality bin')
mc.GetYaxis().SetTitle('PDF')
mc.GetYaxis().SetRangeUser(0,0.1)
mcw=cnv.GetListOfPrimitives().At(1)
data=cnv.GetListOfPrimitives().At(2)

leg=ROOT.TLegend(0.15,0.94,0.4,0.6)
leg.SetBorderSize(0)
leg.SetTextFont(42)
leg.SetTextSize(0.05)
leg.SetFillStyle(0)
leg.AddEntry(data,'Data','ep')
leg.AddEntry(mc,'MC','l')
leg.AddEntry(mcw,'MC weighted','l')
leg.Draw()

txt=ROOT.TLatex()
txt.SetNDC(True)
txt.SetTextFont(42)
txt.SetTextSize(0.045)
txt.SetTextAlign(12)
txt.DrawLatex(0.12,0.97,'#bf{CMS} #it{preliminary}')
cnv.Modified()
cnv.Update()
cnv.SaveAs('cenweight.pdf')
input()


