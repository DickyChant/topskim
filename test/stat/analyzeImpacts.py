import json
import sys
import math

with open(sys.argv[1]) as json_file:
    impacts = json.load(json_file)


uncs={}
for param in impacts['params']:
    name=param['name']
    r0=param['r'][1]
    dr=max([abs(param['r'][0]-r0),abs(param['r'][2]-r0)])

    systName=name.replace('_lnN','')
    if 'prop' in systName:
        systName='Background and \\ttbar signal distribution'
    if 'pdf' in systName or 'muR' in systName or 'muF' in systName:
        systName='nPDF, $\\mu_\\text{R}/\\mu_\\text_F}$ scales, and $\\alpS(m_{\\cPZ}$'
    if systName in ['data_comb','tW','VV','zg']:
        systName='Background normalization'
    if systName in ['ptZ','ptTop']:
        systName='Top quark and Z boson \\pt modelling'
    if systName in ['lepSF','lepIsoSF']:
        systName='Lepton selection efficiency'
    if systName in ['JER','JEC']:
        systName='Jet energy scale and resolution'
    if 'btag' in systName or 'quenching' in systName:
        systName='b tagging efficiency'
    if not systName in uncs: uncs[systName]=0    
    uncs[systName] += dr**2

for s in uncs:
    print('%100s & %3.3f'%(s,math.sqrt(uncs[s])))
