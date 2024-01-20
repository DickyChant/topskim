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


def doNLLScan(data,pList,pTitle,sigma=5):

    """ simple 1 parameter scan """

    plt.clf()
    fig, ax = plt.subplots()
    fig.set_size_inches(10,10)
    ax.set_ylim(0,10)

    x_ran=[0,1]
    for d in data:
        #raw values
        x=data[d][:,0]
        z=data[d][:,1]
        z=2*z

        #interpolate to generate equally spaced grid
        #and apply a gaussian filter
        minX,maxX=min(x),max(x)
        x_unif = np.arange(minX,maxX,0.001*(maxX-minX))
        z_spline = interp1d(x,z,kind='cubic',fill_value='extrapolate')
        z_spline_val=z_spline(x_unif)
        z_filt = filters.gaussian_filter1d(z_spline_val,sigma=sigma)

        #determine the maximum for which the likelihood passes the value of 10
        for i in range(len(x_unif)):
            if z_filt[i]<10:
                x_ran[1]=max([x_ran[1],x_unif[i]])
        
        #plt.plot(x,  z,                   'o')
        ls='-' if 'obs' in d or 'Obs' in d else '--'
        plt.plot(x_unif, z_filt, ls, label=d,linewidth=2)

    plt.xlim(x_ran)
    for cl,clname in [(1, r'68 %CL'),(3.84, r'95 %CL')]:
        plt.plot(x_ran, [cl,cl], ':', color='m')
        ax.text(x_ran[0]+0.12,cl+0.05,clname, fontsize=14)

    ax.legend(framealpha=0.0, fontsize=14, loc='upper right', numpoints=1)

    plt.xlabel(r'%s'%pTitle[0], fontsize=18)
    plt.ylabel(r'$-2\log(\lambda)$',fontsize=18)
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=18)
    ax.text(1.0,1.02,r'1.6 nb$^{-1}$ ($\sqrt{s_{NN}}$=5.02 TeV)', transform=ax.transAxes,horizontalalignment='right',fontsize=16)

    for ext in ['png','pdf']:
        plt.savefig('nllcontour1d_%s.%s'%('_'.join(pList),ext))


def doContour(data,    
              pList,
              pTitles,
              method='linear',
              levels=[2.30,4.61,9.21],
              levelLabels=['68.3%','90%','99%']):

    """ interpolates the grid to obtain the likelihood contour 
    2 parameter fit levels (see PDG Statistics Table 38.2) """

    #add contour level names
    fmt = {}
    for l, s in zip(levels,levelLabels) : fmt[l] = s

    fig, ax = plt.subplots()

    cntr=[]
    for d in data:

        #raw values
        x=data[d][:,0]
        y=data[d][:,1]
        z=data[d][:,2]

        if len(z)<10:
            print('This scan has less than 10 points...')
            return

        #interpolate and find minimum
        xi = np.linspace(min(x),max(x),100)
        yi = np.linspace(min(y),max(y),100)
        zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method=method, fill_value=99999.)

        minz = zi.min()
        bestFitIdx = np.where(zi==minz)
        try:
            bestX=xi[ bestFitIdx[1][0] ]
            bestY=yi[ bestFitIdx[0][0] ]
        except:
            bestX=xi[0]
            bestY=yi[0]
            print('Failed to find best fit idx',d)
        zi=(zi-minz)*2
        z=(z-minz)*2
        
        linestyles=['-' if 'obs' in d or 'Obs' in d else '--']*3
        cntr.append( plt.contour(xi, yi, zi, levels=levels, linewidths=1.0, linestyles=linestyles) )
        cntr[-1].collections[0].set_label(d) #'$-2\Delta ln(\lambda)$')
        plt.clabel(cntr[-1], cntr[-1].levels[:], inline=True, fmt=fmt, fontsize=10)

        #add best-fit point
        if 'obs' in d or 'Obs' in d :
            plt.plot([bestX], [bestY], '+', mew=4, markersize=12, color='k',label='Best fit')

    plt.xlabel(r'%s'%pTitles[0], fontsize=18)
    plt.ylabel(r'%s'%pTitles[1], fontsize=18)
    ax.text(0,1.02,'CMS preliminary', transform=ax.transAxes, fontsize=18)
    ax.text(1.0,1.02,r'1.6 nb$^{-1}$ ($\sqrt{s_{NN}}$=5.02 TeV)', transform=ax.transAxes,horizontalalignment='right',fontsize=16)
    ax.legend(framealpha=0.0, fontsize=14, loc='lower left', numpoints=1)

    for ext in ['png','pdf']:
        plt.savefig('nllcontour2d_%s.%s'%('_'.join(pList),ext))


def main():

    usage = 'usage: %prog title2:file2 title2:file2 [options]'
    parser = optparse.OptionParser(usage)
    parser.add_option('-p', '--params',          
                      dest='params',       
                      help='parameter list (CSV) [%default]',  
                      default='sigma,eb',
                      type='string')
    parser.add_option('--titles',
                      dest='ptitles',       
                      help='titles for the parameters (CSV) [%default]',  
                      default='$\sigma$ [$\mu b$],$\epsilon_{b}$',
                      type='string')
    (opt, args) = parser.parse_args()

    
    TITLES={'sigma':'$\sigma$ [$\mu b$]',
            'eb':'$\epsilon_{b}$',
            'delta':'$\delta_{medium}$'}

    pList=opt.params.split(',')
    pTitles=[TITLES[x] for x in pList]

    #read the likelihood
    fitres={}
    for ia in range(len(args)):
        tit,url=args[ia].split('=')
        print('Reading',url)
        fIn=ROOT.TFile.Open(url,'READ')
        t=fIn.Get('limit')
        fitres[tit]=[]
        for i in range(t.GetEntriesFast()):
            t.GetEntry(i)
            if t.deltaNLL<0 : continue
            fitres[tit].append( [getattr(t,p) for p in pList]+[t.deltaNLL] )
        fitres[tit]=np.array(fitres[tit])
        print(len(fitres[tit]),'scan points found for',tit)

    #plot the contour interpolating the available points
    if len(pList)==2: 
        doContour(fitres,pList,pTitles)
    else:
        doNLLScan(fitres,pList,pTitles)

if __name__ == "__main__":
    sys.exit(main())
