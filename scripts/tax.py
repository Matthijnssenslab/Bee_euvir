import pandas as pd
import pickle
nodetax = []
with open("TotalBP_lowkmer_9780.tab") as f:
    for line in f:
        if not line.startswith("#"):
            lijstje = [line.strip().split()[0],line.strip().split()[1]]
            nodetax.append(lijstje)
taxpathlist = []
totlen = len(nodetax)
count = 0
nodes = pd.read_csv("nodesRED.dmp",delimiter='\t',encoding='utf-8')
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
    print(str(count) + ' out of: ' + str(totlen))
with open('totalscaffolds_size500_9780.pkl', 'wb') as f:
    pickle.dump(taxpathlist, f)
