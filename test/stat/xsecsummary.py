import ROOT
from math import sqrt

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

c=ROOT.TCanvas('c','c',900,700)
c.SetBottomMargin(0.12)
c.SetLeftMargin(0.05)
c.SetRightMargin(0.05)
c.SetTopMargin(0.03)

A=208
P2PPBSCALE=1.#A*1e-3
E8p16TO8=1.#266.24/252.90
C90=1.645#in case the input not already scaled down

datasetleg={'PbPb': '#scale[1.1]{PbPb, 1.7 nb^{-1}, (#sqrt{s_{NN}}=5.02 TeV)}',
            'pp' : '#scale[1.1]{pp, 27.4 pb^{-1}, (#sqrt{s}=5.02 TeV)}'}

theo={
    'PbPb' : [
        #('#left(#splitline{CT14+EPPS16}{NLO MCFM}#right) #upoint K_{NNLO+NNLL}(Top++)',     59.0, 5.3, 5.3, 2.0, 1.6),
        #('#left(#splitline{CT10+EPS09}{NLO MCFM}#right) #upoint K_{NNLO+NNLL}(Top++)',      57.5, 3.3, 4.3, 2.0, 1.5)
        ('#splitline{NNPDF30 NNLO}{NNLO+NNLL Top++}',   68.9, 2.5, 2.7, 2.3, 1.9),  
        ('#splitline{CT14 NNLO}{NNLO+NNLL Top++}',      70.2, 3.2, 4.0, 2.4, 1.9) 
        
        ],
    'pp'  : [
        ('#splitline{NNPDF30 NNLO}{NNLO+NNLL Top++}',  68.9, 2.5, 2.7, 2.3, 1.9),
        ('#splitline{CT14 NNLO}{NNLO+NNLL Top++}',     70.2, 3.2, 4.0, 2.4, 1.9),
        ]
    }

refXsec=68.9
exp={'PbPb':
         [ 
             ('2l_{OS}',           refXsec*0.76,  0, 0, refXsec*0.24, refXsec*0.22),
             ('2l_{OS}+b-tags',    refXsec*0.59,  0, 0, refXsec*0.21, refXsec*0.18),
         ],
     'pp':
        [
            ('2l_{OS}/l+jets', 69.5, 6.1, 6.1, 5.8, 5.8),
        ]
     }
print(exp['PbPb'])
#datasetleg={'pPb': 'pPb, 180 nb^{-1}, (#sqrt{s_{NN}}=8.16 TeV)',
#            'pp' : 'pp, 19.6 fb^{-1}, (#sqrt{s}=8 TeV)'}
#
#theo={
#    'pPb' : [
#        ('#left(#splitline{CT14+EPPS16}{#scale[0.8]{NLO MCFM}}#right) #upoint K_{NNLO}(Top++)',     55.78, 4.91,                  5.02,                   1.95, 1.45),
#        ('#left(#splitline{CT10+EPS09}{#scale[0.8]{NLO MCFM}}#right) #upoint K_{NNLO}(Top++)',      57.56, sqrt(3.10**2+1.21**2), sqrt(3.79**2+2.13**2),  2.01, 1.50),
#        ('#left(#splitline{NNPDF3.0+EPS09}{#scale[0.8]{NLO MCFM}}#right) #upoint K_{NNLO}(Top++)',  56.63, sqrt(1.31**2+1.19**2), sqrt(1.31**2+2.1**2),   1.98, 1.47)
#        ],
#    'pp'  : [
#        ('#splitline{CT14}{#scale[0.8]{NNLO Top++}}',     266.24, 14.95/C90, 16.92/C90, 9.31,  6.84),
#        ('#splitline{MMHT14}{#scale[0.8]{NNLO Top++}}',   265.42, 7.38,      5.41,      9.23,  6.81),
#        ('#splitline{NNPDF3.0}{#scale[0.8]{NNLO Top++}}', 263.77, 6.09,      6.09,      9.23,  6.81)
#        ]
#    }
#
#exp={'pPb':
#         [ 
#        ('#mu+jets',  48, 11, 11, 13.8, 13.8),
#        ('e+jets',    66, 18, 18, 20.8, 20.8),
#        ('l+jets',    54, 11, 11, 13.3, 13.3)
#        ],
#     'pp':
#         [
#        ('e#mu #scale[0.6]{JHEP 1608 (2016) 029}', 244.90, 1.4, 1.4, 9.089, 8.554),
#        ('l+jets #scale[0.6]{EPJC 77 (2017) 15}',  228.5,  3.8, 3.8, 15.43, 15.43)
#        ]
#     }
#
colors  = {'PbPb':1,          'pp':1} #ROOT.TColor.GetColor('#fc8d59')}
fill    = {'PbPb':ROOT.kAzure-2, 'pp':ROOT.kGreen-5} 
markers = {'PbPb':20,         'pp':20}

dy=3.75
for key in exp:
    dy+=len(exp[key])
    if key in theo:
        dy+=len(theo[key])-1

frame=ROOT.TH2F('frame',';#sigma #times A^{2} [pb]',1,0,115,1,-1.25,dy)
frame.GetXaxis().SetTitleOffset(0.85)
frame.GetXaxis().SetTitleSize(0.06)
frame.GetXaxis().SetLabelSize(0.04)
frame.GetYaxis().SetNdivisions(0)
frame.Draw()

cms=ROOT.TLatex()
cms.SetTextFont(42)
cms.SetTextSize(0.06)
cms.SetNDC()
cms.DrawLatex(0.82,0.9,'#bf{CMS}') # #it{Preliminary}')

labels=ROOT.TLatex()
labels.SetTextFont(42)
labels.SetTextSize(0.05)

theoPDFGr,theoTotGr=[],[]
expStatGr,expTotalGr=[],[]
imeas=0

dy=0
for key in exp:

    if imeas>0: dy+=1
    expdy=len(exp[key])
    theody=0

    #check if there are theory references acompanying
    if key in theo:

        theody=len(theo[key])

        for i in range(0,len(theo[key])):
            title,mu,pdfDn,pdfUp,scaleDn,scaleUp=theo[key][i]
            if key=='pp':
                mu      *= P2PPBSCALE
                pdfUp   *= P2PPBSCALE
                pdfDn   *= P2PPBSCALE
                scaleUp *= P2PPBSCALE
                scaleDn *= P2PPBSCALE
            totDn=ROOT.TMath.Sqrt(pdfDn**2+scaleDn**2)
            totUp=ROOT.TMath.Sqrt(pdfUp**2+scaleUp**2)
            print(key,mu,totDn,totUp)
            theoPDFGr.append( ROOT.TGraph() )
            theoPDFGr[-1].SetTitle(title)
            theoPDFGr[-1].SetName('theopdf_%s_%d'%(key,i))
            #theoPDFGr[-1].SetFillColor(fill[key]-i)
            theoPDFGr[-1].SetFillColorAlpha(fill[key]-i, 0.5);
            theoPDFGr[-1].SetFillStyle(1001)

            theoTotGr.append( theoPDFGr[-1].Clone('theotot_%s_%d'%(key,i)) )
            #theoTotGr[-1].SetFillColor(fill[key]-i)
            theoTotGr[-1].SetFillColorAlpha(fill[key]-i, 0.5);
            theoTotGr[-1].SetFillStyle(3001)

            ymin,ymax=dy,dy+expdy+1
            if i>0: ymin,ymax=dy+expdy+i,dy+expdy+1+i

            theoPDFGr[-1].SetPoint(0,mu-pdfDn,ymin)
            theoPDFGr[-1].SetPoint(1,mu-pdfDn,ymax)
            theoPDFGr[-1].SetPoint(2,mu+pdfUp,ymax)
            theoPDFGr[-1].SetPoint(3,mu+pdfUp,ymin)
            theoPDFGr[-1].SetPoint(4,mu-pdfDn,ymin)

            theoTotGr[-1].SetPoint(0,mu-totDn,ymin)
            theoTotGr[-1].SetPoint(1,mu-totDn,ymax)
            theoTotGr[-1].SetPoint(2,mu+totUp,ymax)
            theoTotGr[-1].SetPoint(3,mu+totUp,ymin)
            theoTotGr[-1].SetPoint(4,mu-totDn,ymin)
                
            theoTotGr[-1].Draw('f')
            theoPDFGr[-1].Draw('f')

            labels.DrawLatex(88,ymax-0.75,'#scale[0.6]{%s}'%title)

    dy+=expdy+theody
    labels.DrawLatex(5,dy-0.5,'#scale[0.7]{#bf{%s}}'%datasetleg[key])
    #if key=='pp':
    #    labels.DrawLatex(5,dy-1.2,'#scale[0.5]{#it{Data scaled by A #upoint #frac{#sigma_{NNLO+NNLL}(8.16 TeV)}{#sigma_{NNLO+NNLL}(8 TeV)}}}')

    ROOT.gStyle.SetEndErrorSize(5)
    for i in range(0,len(exp[key])):
        expStatGr.append( ROOT.TGraphAsymmErrors() )
        expStatGr[-1].SetMarkerStyle(0)
        expStatGr[-1].SetLineWidth(3)
        expStatGr[-1].SetLineColor(colors[key])
        expStatGr[-1].SetFillStyle(0)

        expTotalGr.append( ROOT.TGraphAsymmErrors() )
        expTotalGr[-1].SetLineWidth(3)
        expTotalGr[-1].SetMarkerStyle(markers[key])
        expTotalGr[-1].SetMarkerColor(colors[key])
        expTotalGr[-1].SetLineColor(colors[key])
        expTotalGr[-1].SetFillStyle(0)

        title,xsec,statUp,statDn,totalUp,totalDn=exp[key][i]
        if key=='pp':
            xsec    *= P2PPBSCALE*E8p16TO8
            statUp  *= P2PPBSCALE*E8p16TO8
            statDn  *= P2PPBSCALE*E8p16TO8
            totalUp *= P2PPBSCALE*E8p16TO8
            totalDn *= P2PPBSCALE*E8p16TO8

        yanchor=imeas+i+1 if imeas==0 else imeas+i

        expStatGr[-1].SetPoint(0,xsec,yanchor)
        expStatGr[-1].SetPointError(0,statUp,statDn,0,0)
        expStatGr[-1].Draw('p[]')

        expTotalGr[-1].SetPoint(0,xsec,yanchor)
        expTotalGr[-1].SetPointError(0,totalUp,totalDn,0,0)
        expTotalGr[-1].Draw('p')
        #if "ell" not in title:
        #    print title
        #    labels.DrawLatex(5,yanchor,'#scale[0.6]{%s}'%title)
        #else:
        #    print 'else ', title
        labels.DrawLatex(5,yanchor,title)

    imeas += dy+2        

grlegstat=ROOT.TGraphErrors()
grlegstat.SetMarkerStyle(0)
grlegstat.SetLineWidth(2)
grlegstat.SetLineColor(1) #ROOT.kGray)
grlegstat.SetFillStyle(0)
grlegstat.SetPoint(0,97.5,0.5)
grlegstat.SetPointError(0,2.5,0)

grlegtot=ROOT.TGraphErrors()
grlegtot.SetLineWidth(2)
grlegtot.SetLineColor(1) #ROOT.kGray)
grlegtot.SetFillStyle(0)
grlegtot.SetPoint(0,97.5,0.5)
grlegtot.SetPointError(0,8,0)
grlegtot.SetMarkerStyle(20)
grlegtot.SetMarkerColor(1)

thlegpdf=ROOT.TGraphErrors()
thlegpdf.SetMarkerStyle(0)
#thlegpdf.SetFillColor(ROOT.kGray)
#thlegpdf.SetFillStyle(1001)
thlegpdf.SetFillColorAlpha(ROOT.kGray,0.5)
thlegpdf.SetFillStyle(1001)
thlegpdf.SetPoint(0,97.5,-0.5)
thlegpdf.SetPointError(0,2.5,0.25)

thlegtot=ROOT.TGraphErrors()
thlegtot.SetMarkerStyle(0)
thlegtot.SetFillColorAlpha(ROOT.kGray,0.5)
thlegtot.SetFillStyle(3001)
#thlegtot.SetFillColor(ROOT.kGray)
#thlegtot.SetFillStyle(3001)
thlegtot.SetPoint(0,97.5,-0.5)
thlegtot.SetPointError(0,8,0.25)


grlegstat.Draw('p[]')
grlegtot.Draw('p')
labels.SetTextSize(0.045)
#labels.DrawLatex(84,0.75,'#scale[0.6]{#color[15]{#it{Exp. unc.: stat  stat#oplussyst}}}')
labels.DrawLatex(84,0.75,'#scale[0.6]{#it{Exp. unc.: stat  stat#oplussyst}}')

thlegpdf.Draw('e2')
thlegtot.Draw('e2')
labels.SetTextSize(0.045)
#labels.DrawLatex(85.5,-0.15,'#scale[0.6]{#color[15]{#it{Th. unc.: pdf  pdf#oplusscales}}}')
labels.DrawLatex(85.5,-0.15,'#scale[0.6]{#it{Th. unc.: pdf  pdf#oplusscales}}')

#labels.DrawLatex(5,-0.6,'#scale[0.7]{#it{All theory predictions use: #mu_{F}=#mu_{R}=m_{t}=172.5 GeV, #alpha_{s} = 0.1180}}')

c.RedrawAxis()
c.Modified()
c.Update()
c.SaveAs('xsecsummary.pdf')
c.SaveAs('xsecsummary.png')


