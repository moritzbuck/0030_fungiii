from ete3 import NCBITaxa

ncbi = NCBITaxa()

with open("data/metadat.csv") as handle:
    name2file = { l.split("\t")[0] : l.split("\t")[1][:-1]  for l in handle}

species = list(set([" ".join(n.split()[0:2]).replace(" sp.","") for n in name2file.keys()]))
f2species = { "_".join(v.split("_")[:-2]) : " ".join(k.split()[0:2]).replace(" sp.","") for k, v in name2file.items()}

spec2tax = {k : dict([ (ncbi.get_rank([t])[t], ncbi.get_taxid_translator([t])[t]) for t in ncbi.get_lineage(v[0])]) for k,v in ncbi.get_name_translator(species).items()}
levels = ['superkingdom','phylum','class', 'order','family','genus','species']
levels_long = ['superkingdom','phylum','class', 'order','family','genus','species']
f2taxstr = {k : ",".join([spec2tax[v].get(l,"") for l in levels]) for k,v in f2species.items()}

with open("data/fungiii.tax","w") as handle:
    handle.write(",".join(["ID"] + levels) + "\n")
    handle.writelines([k + ".fasta," + v + "\n" for  k,v in  f2taxstr.items()])

    
