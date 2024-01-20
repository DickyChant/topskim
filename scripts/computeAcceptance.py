import ROOT
import sys
import numpy as np

wgtList={
    'qcdscale':list(range(1,7)),
    'nnpdf':list(range(8,110)),
    'alpS':list(range(110,112)),
    'epps16ct14':list(range(9,9+96))
}

url=sys.argv[1]
fIn=ROOT.TFile.Open(url)
isPbPb=True if 'Drum' in url else False 

h=fIn.Get('gen_fidcounter')
h.SetDirectory(0)

hgen=h.ProjectionY("hgen",1,1)
hstep=[]
for xbin in range(1,h.GetNbinsX()+1):
    hstep.append( h.ProjectionY("hstep%d"%xbin,xbin,xbin) )
    hstep[-1].Divide(hgen)
    hstep[-1].SetDirectory(0)
    hstep[-1].Scale(100.)

    nomVal=hstep[-1].GetBinContent(1)
    print('%d %20s %3.2f'%(xbin,h.GetXaxis().GetBinLabel(xbin),nomVal), end=' ')

    qcdDiff=[abs(hstep[-1].GetBinContent(i+1)-nomVal) for i in wgtList['qcdscale']]
    print('+/-%3.2f (scales)'%max(qcdDiff), end=' ')

    if isPbPb:
        npdfDiff=[hstep[-1].GetBinContent(i+1)-nomVal for i in wgtList['epps16ct14']]
        deltaO_pVec=[ npdfDiff[i] for i in range(0,len(npdfDiff),2) ]
        deltaO_mVec=[ npdfDiff[i] for i in range(1,len(npdfDiff),2) ]
        deltaO_p,deltaO_m=0.,0.
        for i in range(len(deltaO_pVec)):
            deltaO_p += max(deltaO_pVec[i],deltaO_mVec[i])**2
            deltaO_m += min(deltaO_pVec[i],deltaO_mVec[i])**2
        deltaO_p=np.sqrt(deltaO_p)
        deltaO_m=np.sqrt(deltaO_m)
        print('+%3.2f-%3.2f (nPDFs)'%(deltaO_p,deltaO_m), end=' ')
    else:
        pdfDiff=sum([(hstep[-1].GetBinContent(i+1)-nomVal)**2 for i in wgtList['nnpdf']])
        alpsDiff=max([(hstep[-1].GetBinContent(i+1)-nomVal)**2 for i in wgtList['alpS']])         
        print('%3.2f (PDFs)'%(np.sqrt(pdfDiff)), end=' ')


#    for w in wgtList:
#        continue
#        diff=
#        
#        if w=='qcdscale':
#            print '+',max(diff),min(diff),'(scale)',
#        if w=='pdf':
#            print '+/-',np.std([abs(x) for x in diff]),'(pdf)',

    print('')
    
#
#
#hiso=fIn.Get('lep_fidcounter')
#hiso.SetDirectory(0)
#hfakeiso=fIn.Get('fakelep_fidcounter')
#hfakeiso.SetDirectory(0)
#fIn.Close()
#
#hiso.Add(hfakeiso)
#hfakeiso.Divide(hiso)
#fakefrac=hfakeiso.ProjectionY('fakefrac',3,3)
#fakefrac.Scale(100.)
#
#nomVal=fakefrac.GetBinContent(1)
#print 'fake fraction',nomVal,
#for w in wgtList:
#    diff=[fakefrac.GetBinContent(i+1)-nomVal for i in wgtList[w]] 
#    if w=='qcdscale':
#        print '+',max(diff),min(diff),'(scale)',
#    if w=='pdf':
#        print '+/-',np.std([abs(x) for x in diff]),'(pdf)',
#print ''
