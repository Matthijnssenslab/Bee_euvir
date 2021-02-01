import pandas as pd
import pickle
import argparse

parser =argparse.ArgumentParser(description="parse tabfiles / pickled lists containing accesionID,taxID")
parser.add_argument('-i', help='give input', required=True)
parser.add_argument('-t', help='tab of pickle', choices=['tab','pickle'], required=True)
parser.add_argument('--nodes', help='location of nodesRED.dmp', required=True)
parser.add_argument('-o', help='output filename', required=True)

args = parser.parse_args()
if args.t == 'tab':
    nodetax = []
    with open(args.i) as f:
        for line in f:
            if not line.startswith('#'):
                temLis = [line.strip().split()[0],line.strip().split()[1]]
                nodetax.append(temLis)
elif args.t == 'pickle':
    with open(args.i, 'rb') as f:
        nodetax = pickle.load(f)

#with open("TotalBP_lowkmer_9780.tab") as f:
#    for line in f:
#        if not line.startswith("#"):
#            lijstje = [line.strip().split()[0],line.strip().split()[1]]
#            nodetax.append(lijstje)
taxpathlist = []
totlen = len(nodetax)
count = 0
nodes = pd.read_csv(args.nodes,delimiter='\t',encoding='utf-8')
keepranks = ['superkingdom','species','genus','family','subfamily','order', 'root']
for i in nodetax:
    name = i[0]
    tax = i[1]
    templis = [name]
    taxlis = [tax]
    while tax != 1:
        newnod = nodes.loc[nodes['tax'] == int(tax)]
        try:
            tax = newnod.iat[0,1]
        except:
            tax = 1
        taxlis.append(tax)
    taxlis = taxlis[::-1]
    finlis = templis + taxlis
    taxpathlist.append(finlis)
    count += 1
    if count % 200 == 0:
        print(str(count) + ' out of: ' + str(totlen))
print("Parsed through {}".format(args.i))
with open(args.o, 'wb') as f:
    pickle.dump(taxpathlist, f)
