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

def getNLLResults(urlscan,urlsingles):

    """gets the results of the scan and fit"""

    fitres=[]
    scan=[]

    #likelihood scan
    fIn=ROOT.TFile.Open(urlscan,'READ')
    t=fIn.Get('limit')
    for k in range(t.GetEntriesFast()):
        t.GetEntry(k)
        if t.quantileExpected<0: continue
        if t.deltaNLL<0 : continue
        scan.append( [t.r,t.deltaNLL] )
    scan=np.array(scan)
    fIn.Close()
    
    #fit result
    fIn=ROOT.TFile.Open(urlsingles,'READ')
    t=fIn.Get('limit')
    for i in range(3):
        t.GetEntry(i)
        fitres.append(t.r)
    fitres[1]=fitres[1]-fitres[0]
    fitres[2]=fitres[2]-fitres[0]
    fIn.Close()
    

    return scan,fitres

def getResultsFromToys(url,branch='limit'):

    """ reads out the results of the toys """

    sigToys=[]

    try:
        fIn=ROOT.TFile.Open(url,'READ')
        t=fIn.Get('limit')
        for k in range(t.GetEntriesFast()):
            t.GetEntry(k)
            if isinstance(branch, str):
                sigToys.append(getattr(t,branch))
            else:
                sigToys.append( [getattr(t,b) for b in branch] )
    except:
        print('Check',url)
        pass

    return sigToys



def finalizePlot(plt,ax,name,title='',banner='1.7 nb$^{-1}$ ($\sqrt{s_{NN}}$=5.02 TeV)'):

    """ final cosmetics on the plot """
    
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=22,weight='bold')
    ax.text(1.0,1.02,r'%s'%banner, transform=ax.transAxes,horizontalalignment='right',fontsize=20)
    if len(title)>0:
        ax.text(0.05,0.95,r'[%s]'%title, transform=ax.transAxes,horizontalalignment='left', fontsize=20)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for ext in ['png','pdf']:
        plt.savefig('%s.%s'%(name,ext))


def showToys(toys,obs,name,title,xtitle,xbins,xran,yran):

    """compares the results of the toys with the observation"""

    plt.clf()
    fig, ax = plt.subplots()
    fig.set_size_inches(10,10)

    ax.hist(toys,range=xran, bins=xbins, histtype='step', label=r'Exp. (median) %3.2f'%(np.median(toys)),linewidth=2)
    pval=0.
    for istoy in toys: 
        if istoy>obs: continue
        pval+=1.0
    pval /= len(toys)
    if pval>0.5 : pval=1-pval
    print(obs,yran,pval)
    plt.plot([obs,obs], [yran[0],yran[1]*0.8], ls='-', linewidth=3, label=r'Obs. %3.2f (p-val=%3.2f)'%(obs,pval))

    ax.set_ylim(yran[0],yran[1])
    plt.xlabel(r'%s'%xtitle, fontsize=22)
    plt.ylabel(r'Toys',fontsize=22)
    ax.legend(framealpha=0.0, fontsize=20, loc='upper right', numpoints=1)

    finalizePlot(plt,ax,name,title)


def showNLLScan(exp,obs,name,title,xran=[0,2],yran=[0,25]):

    """displays the result of the likelihood scan"""
    
    plt.clf()
    fig, ax = plt.subplots()
    fig.set_size_inches(10,10)
    ax.set_ylim(yran[0],yran[1])

    plt.plot(exp[0][:,0], 2*exp[0][:,1], ls='--', label=r'Exp. $%3.2f^{+%3.2f}_{%3.2f}$'%(exp[1][0],exp[1][2],exp[1][1]),linewidth=2)
    plt.plot(obs[0][:,0], 2*obs[0][:,1], ls='-', label=r'Obs. $%3.2f^{+%3.2f}_{%3.2f}$'%(obs[1][0],obs[1][2],obs[1][1]),linewidth=3)
    
    plt.xlim(xran)
    for cl,clname in [(1, r'68 %CL'),(3.84, r'95 %CL')]:
        plt.plot(xran, [cl,cl], ':', color='m')
        ax.text(xran[1]-0.05,cl+0.05,clname,horizontalalignment='right', fontsize=16)

    ax.legend(framealpha=0.0, fontsize=20, loc='upper right', numpoints=1)

    plt.xlabel(r'$\mu=\sigma/\sigma_{th}$', fontsize=26)
    plt.ylabel(r'$-2\log(\lambda)$',fontsize=26)

    finalizePlot(plt,ax,name,title)


def showNLLContour(exp,obs,name,title,method='linear',levels=[2.30,4.61,9.21],levelLabels=['68.3%','90%','99%']):

    """ interpolates the grid to obtain the likelihood contour 2 parameter fit levels (see PDG Statistics Table 38.2) """

    plt.clf()
    fig, ax = plt.subplots() 
    fig.set_size_inches(10,10)

    cntr=[]
    for data,ls,c,label in [(exp,'--','gray','Exp.'),(obs,'-','k','Obs.')]:

        data=np.array(data)
        x=data[:,0]
        y=data[:,1]
        z=data[:,2]

        #interpolate in a regular grid and draw the contour
        xi = np.linspace(min(x),max(x),1000)
        yi = np.linspace(min(y),max(y),1000)
        zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method=method)
        cntr.append( ax.contour(xi, yi, 2*zi, levels=levels, linewidths=1.0, linestyles=ls, colors=c) )
        cntr[-1].collections[0].set_label(r'%s'%label)

        #add contour level names
        fmt = {}
        for l, s in zip(levels,levelLabels):
            fmt[l] = s
        ax.clabel(cntr[-1], cntr[-1].levels[:], inline=True, fmt=fmt, fontsize=26)

    plt.xlabel(r'$\mu_{t\bar{t}}$', fontsize=28)
    plt.ylabel(r'$\mu_{DY}$',       fontsize=28)
    plt.ylim(0.6,1.4)
    plt.xlim(0,2)
    ax.legend(framealpha=0.0, fontsize=20, loc='upper right', numpoints=1)
    finalizePlot(plt,ax,name,title)


""" deprecated? """

def showSignificanceScan(scan,pfix='',lumi=0):
    plt.clf()
    fig, ax = plt.subplots()
    fig.set_size_inches(10,10)
    x_ran=[200,2000]
    ax.set_xlim(x_ran)
    ax.set_ylim((1,6))
    ithick=2.5
    for _,tit,d in scan:
        if ithick>1.0:ithick-=0.5
        plt.plot(d[:,0],d[:,1],label=r'%s'%tit,linewidth=ithick)
    ax.legend(framealpha=0.0, fontsize=20, loc='upper left', numpoints=1)
    for cl,clname in [(3, r'3$\sigma$'),(5, r'5$\sigma$')]:
        plt.plot(x_ran, [cl,cl], ':', color='m')
        ax.text(x_ran[1]-0.05,cl+0.05,clname,horizontalalignment='right', fontsize=18)
    plt.plot([458,458], [0,3], '--', color='m')
    plt.plot([1752,1752],   [0,5], '--', color='m')

    plt.xlabel(r'Integrated luminosity nb$^{-1}$', fontsize=22)
    plt.ylabel(r'Significance (asymptotic)',fontsize=22)
    finalizePlot(plt,ax,'significancescan_r'+pfix,lumi)


def getSignificanceScan(url):
    sigScan=[('higgsCombinesig_0.5.Significance.mH172.5.root',0.5),
             ('higgsCombinesig_1.Significance.mH172.5.root',1.0),
             ('higgsCombinesig_1.5.Significance.mH172.5.root',1.5),
             ('higgsCombinesig_2.Significance.mH172.5.root',2.0),
             ('higgsCombinesig_2.5.Significance.mH172.5.root',2.5),
             ('higgsCombinesig_3.Significance.mH172.5.root',3.0),
             ('higgsCombinesig_3.5.Significance.mH172.5.root',3.5),
             ('higgsCombinesig_4.Significance.mH172.5.root',4.0),
             ('higgsCombinesig_4.5.Significance.mH172.5.root',4.5),
             ('higgsCombinesig_5.Significance.mH172.5.root',5.0)]
    sigScanVals=[]
    for f,lumi in sigScan:
        fIn=ROOT.TFile.Open(url+'/'+f,'READ')
        t=fIn.Get('limit')
        t.GetEntry(0)
        sigScanVals.append( [lumi*458,t.limit] )
        fIn.Close()
    return np.array(sigScanVals)
