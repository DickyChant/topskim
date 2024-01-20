import os
import sys

td = 'treeProducerHIN'
url=sys.argv[1]

for sd in os.listdir(url):
    if not '.root' in sd: continue
    dname = sd.replace('.root','')
    cmd_mkdir = 'mkdir -p {url}/{dn}/{td}'.format(url=url,dn=dname, td=td)
    cmd_mv    = 'cp {url}/{sd} {url}/{dn}/{td}/tree.root'.format(url=url,sd=sd,dn=dname,td=td)

    print('==============================')
    print(cmd_mkdir)
    os.system(cmd_mkdir)
    print(cmd_mv)
    os.system(cmd_mv)


combDir='{url}/Combinatorial'.format(url=url)
for sd in os.listdir(combDir):
    dname=os.path.splitext(sd)[0]
    cmd_mkdir = 'mkdir -p {url}/{dn}_Combinatorial/{td}'.format(url=url,dn=dname, td=td)
    cmd_mv    = 'cp {url}/Combinatorial/{sd} {url}/{dn}_Combinatorial/{td}/tree.root'.format(url=url,sd=sd,dn=dname,td=td)
    print('==============================')
    print(cmd_mkdir)
    os.system(cmd_mkdir)
    print(cmd_mv)
    os.system(cmd_mv)
