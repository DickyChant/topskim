import ROOT

def doRinRout(nin,nout,ndyin,ndyout,ch1,ch2):

    """
    implementation of RinRout a la TOP-16-005
    expect dicts with counts+/-unc  where keys are ee,mm,em
    ch1 is the channel to predict
    ch2 is the channel used as control
    """
    
    rinll    = nin[ch1][0]/nin[ch2][0]
    rinllUnc = rinll*ROOT.TMath.Sqrt(1./nin[ch1][0]+1./nin[ch2][0])

    kll    = ROOT.TMath.Sqrt(rinll)
    kllUnc = 0.5*kll*rinllUnc/rinll
    
    routin    = ndyout[ch1][0]/ndyin[ch1][0]
    routinUnc = routin*ROOT.TMath.Sqrt( (ndyout[ch1][1]/ndyout[ch1][0])**2 
                                        + (ndyin[ch1][1]/ndyin[ch1][0])**2 )

    nllout    = routin*(nin[ch1][0]-0.5*nin['em'][0]*kll)
    nlloutUnc = ROOT.TMath.Sqrt( (routinUnc*(nin[ch1][0]-0.5*nin['em'][0]*kll))**2
                                 + (routin**2)*nin[ch1][0]
                                 + nin['em'][0]*(0.5*routin*kll)**2
                                 + (0.5*routin*nin['em'][0]*kllUnc)**2 )
    
    sf    = nllout/ndyout[ch1][0]
    sfUnc = sf*ROOT.TMath.Sqrt( (nlloutUnc/nllout)**2+(ndyout[ch1][1]/ndyout[ch1][0])**2 )

    #print report
    print('-'*20,ch1,'-'*20)
    print('R_{out/nin}(MC)=%f +/- %f'%(routin,routinUnc))
    print('k_{\\ell,\\ell}=%f +/- %f'%(kll,kllUnc))
    print('N_{\\rm in}(data)=',nin[ch1][0])
    print('N_{\\rm out}=%f +/- %f'%(nllout,nlloutUnc))
    print('SF_{\\rm DY}=%f +/- %f'%(sf,sfUnc))
    print('-'*50)

    return (sf,sfUnc)
