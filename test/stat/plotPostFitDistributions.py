import ROOT
import sys
import os

mcStack=[('total_background', 'Background', ROOT.kGray),
         ('total_signal',     'Signal',     ROOT.kRed+2)]

def doPostFitPlot(url):

   ROOT.gStyle.SetOptTitle(0)
   ROOT.gStyle.SetOptStat(0)
   ROOT.gROOT.SetBatch(True)

   inF = ROOT.TFile(url)
   xtitle='BDT' if 'bdtcomb' in url else 'Sphericity'
   xtitle += ' bin'
   for ch in inF.Get('shapes_fit_s').GetListOfKeys():
      chName=ch.GetName()

      chTitle=chName.replace('mm','#mu#mu')
      chTitle=chTitle.replace('em','e#mu')
      if url.find('_jetAnalysis')>0 : chTitle += '+b-tags'

      chDir=ch.ReadObj()
      plotsPrefit={}
      plotsPostfit={}
      for proc in chDir.GetListOfKeys():
         pname=proc.GetName()
         plotsPostfit[pname]=proc.ReadObj()
         try:
            plotsPostfit[pname].SetDirectory(0)
         except:
            pass
         plotsPrefit[pname]=inF.Get('shapes_prefit/%s/%s'%(chName,pname))
         try:
            plotsPrefit[pname].SetDirectory(0)
         except:
            pass

      compareFitResult(plotsPostfit=plotsPostfit,plotsPrefit=plotsPrefit,
                       plotName=os.path.dirname(url)+'/postfit_'+chName,
                       xtitle=xtitle,extraTxt=[chTitle])
    
def compareFitResult(plotsPrefit,plotsPostfit,plotName,xtitle,extraTxt=[]):

   c = ROOT.TCanvas("c","c",600,600)
   c.SetTopMargin(0)
   c.SetLeftMargin(0)
   c.SetRightMargin(0)
   c.SetBottomMargin(0)

   #data/MC
   p1 = ROOT.TPad("p1", "p1", 0., 0.5, 1., 1.)
   p1.SetTopMargin(0.14)
   p1.SetRightMargin(0.03)
   p1.SetLeftMargin(0.12)
   p1.SetBottomMargin(0.01)
   p1.Draw()
   p1.cd()
   frame=plotsPostfit['total'].Clone('frame')
   frame.Reset('ICE')
   frame.GetYaxis().SetTitle('Events')
   frame.GetYaxis().SetTitleOffset(0.75)
   frame.GetYaxis().SetNdivisions(5)
   frame.GetYaxis().SetTitleSize(0.08)
   frame.GetYaxis().SetLabelSize(0.08)
   frame.GetXaxis().SetTitleSize(0)
   frame.GetXaxis().SetLabelSize(0)
   frame.GetYaxis().SetRangeUser(1e-1,1.5*plotsPostfit['total'].GetMaximum())
   frame.Draw()

   leg = ROOT.TLegend(0.7,0.55,0.97,0.83)
   leg.SetFillStyle(0)
   leg.SetBorderSize(0)
   leg.SetTextSize(0.08)
   leg.AddEntry(plotsPostfit['data'],'Data','ep')

   stack = ROOT.THStack("stack","stack")
   ratios=[]
   for proc,title,ci in mcStack:
      plotsPostfit[proc].SetTitle(title)
      plotsPostfit[proc].SetLineColor(1)
      plotsPostfit[proc].SetFillStyle(1001)
      plotsPostfit[proc].SetFillColor(ci)
      stack.Add(plotsPostfit[proc])
      leg.AddEntry(plotsPostfit[proc],title,'f')

      ratios.append( plotsPostfit[proc].Clone(proc+'_2prefit') )
      ratios[-1].SetDirectory(0)
      ratios[-1].Divide( plotsPrefit[proc] )
      ratios[-1].SetFillStyle(0)
      ratios[-1].SetFillColor(0)
      ratios[-1].SetLineColor(ci)
      ratios[-1].SetLineWidth(3)

   stack.Draw('histsame')

   totalUnc=ROOT.TGraphErrors(plotsPostfit['total'])
   totalUnc.SetFillStyle(3444)
   totalUnc.SetFillColor(1)
   totalUnc.SetMarkerStyle(1)
   totalUnc.Draw('e2')

   plotsPostfit['data'].SetMarkerStyle(20)
   for i in range(plotsPostfit['data'].GetN()):
      plotsPostfit['data'].SetPointEXhigh(i,0)
      plotsPostfit['data'].SetPointEXlow(i,0)
   plotsPostfit['data'].Draw('p')

   leg.Draw('same')

   tex=ROOT.TLatex()
   tex.SetTextFont(42)
   tex.SetTextSize(0.08)
   tex.SetNDC()
   tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
   tex.DrawLatex(0.12,0.93,'#bf{CMS} #it{preliminary}')
   tex.SetTextAlign(ROOT.kHAlignRight+ROOT.kVAlignCenter)
   tex.DrawLatex(0.95,0.93,'1.7 nb^{-1} (#sqrt{s_{NN}}=5.02 TeV)')
   tex.SetTextAlign(ROOT.kHAlignLeft+ROOT.kVAlignCenter)
   for it in range(len(extraTxt)):
      tex.DrawLatex(0.2,0.8+it*0.05,extraTxt[it])

   p1.RedrawAxis()

   #prefit/postfit
   c.cd()

   p2 = ROOT.TPad("p2", "p1", 0., 0.3, 1., 0.5)
   p2.SetTopMargin(0.01)
   p2.SetRightMargin(0.03)
   p2.SetLeftMargin(0.12)
   p2.SetBottomMargin(0.01)
   p2.SetGridy()
   p2.Draw()
   p2.cd()

   frame2=frame.Clone()
   frame2.GetYaxis().SetRangeUser(0.5,1.5)
   frame2.GetYaxis().CenterTitle()
   frame2.GetYaxis().SetTitle('Postfit/Prefit')
   frame2.GetYaxis().SetTitleOffset(0.3)
   frame2.GetYaxis().SetTitleSize(0.2)
   frame2.GetYaxis().SetLabelSize(0.2)
   frame2.Draw()

   leg2 = ROOT.TLegend(0.15,0.77,0.65,0.93)
   leg2.SetNColumns(len(ratios))
   leg2.SetFillStyle(0)
   leg2.SetBorderSize(0)
   leg2.SetTextSize(0.2)
   for r in ratios:
      r.Draw('histsame')
      leg2.AddEntry(r,r.GetTitle(),'l')
   leg2.Draw()

   p2.RedrawAxis()

   #data/MC
   c.cd()
   p3 = ROOT.TPad("p3", "p3", 0., 0.0, 1., 0.3)
   p3.SetTopMargin(0.01)
   p3.SetRightMargin(0.03)
   p3.SetLeftMargin(0.12)
   p3.SetBottomMargin(0.3)
   p3.SetGridy()
   p3.Draw()
   p3.cd()

   frame3=frame2.Clone()
   frame3.GetYaxis().SetTitle('Obs./Fit')
   frame3.GetYaxis().SetTitleOffset(0.45)
   frame3.GetYaxis().SetTitleSize(0.13)
   frame3.GetYaxis().SetLabelSize(0.13)
   frame3.GetXaxis().SetTitleSize(0.13)
   frame3.GetXaxis().SetLabelSize(0.13)
   frame3.GetXaxis().SetTitle(xtitle)
   frame3.Draw()

   def getRelUnc(plotColl,name,ci,fill):
      totalNoUnc=plotColl[name].Clone(name+'_nounc')
      for i in range(totalNoUnc.GetNbinsX()):
         totalNoUnc.SetBinError(i+1,0)
      relUnc=plotColl[name].Clone(name+'_relUnc')
      relUnc.Divide(totalNoUnc)
      relUncGr=ROOT.TGraphErrors(relUnc)
      relUncGr.SetFillStyle(fill)
      relUncGr.SetFillColor(ci)
      relUncGr.SetMarkerStyle(1)
      relUnc.Delete()
      totalNoUnc.Delete()
      return relUncGr

   relPreUncGr=getRelUnc(plotsPrefit,'total',ROOT.kGray,1001)
   relPreUncGr.Draw('e2')
   relFitUncGr=getRelUnc(plotsPostfit,'total',1,3344)
   relFitUncGr.Draw('e2')   
   data2fitGr=plotsPostfit['data'].Clone('data2fit')
   x,y=ROOT.Double(0),ROOT.Double(0)
   for i in range(data2fitGr.GetN()):
      den=plotsPostfit['total'].GetBinContent(i+1)
      if float(den)==0 : 
         data2fitGr.SetPointEYhigh(i,0)
         data2fitGr.SetPointEYlow(i,0)
      else:
         denUnc=plotsPostfit['total'].GetBinError(i+1)
         data2fitGr.GetPoint(i,x,y)
         data2fitGr.SetPoint(i,x,y/den)
         eyhi=data2fitGr.GetErrorYhigh(i)
         eylo=data2fitGr.GetErrorYlow(i)
         data2fitGr.SetPointEYhigh(i,eyhi/den)
         data2fitGr.SetPointEYlow(i,eylo/den)
         data2fitGr.Draw('p0')

   leg3 = ROOT.TLegend(0.15,0.85,0.8,0.95)
   leg3.SetNColumns(3)
   leg3.SetFillStyle(0)
   leg3.SetBorderSize(0)
   leg3.SetTextSize(0.13)
   leg3.AddEntry(data2fitGr,'Data','ep')
   leg3.AddEntry(relPreUncGr,'Prefit unc.','f')
   leg3.AddEntry(relFitUncGr,'Postfit unc.','f')
   leg3.Draw()

   p3.RedrawAxis()
   
   c.cd()
   c.Modified()
   c.Update()
   for ext in ['png','pdf']:
      c.SaveAs('%s.%s'%(plotName,ext))

   p1.Delete()
   p2.Delete()
   p3.Delete()


def main():
   url=sys.argv[1]
   if not os.path.isfile(url):
      print(url,'is not a file')
      sys.exit()
   doPostFitPlot(url)
   
if __name__ == "__main__":
    sys.exit(main())
