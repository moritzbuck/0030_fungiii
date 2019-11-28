import os, sys
from tqdm import tqdm
from os.path import join as pjoin
from ete3 import NCBITaxa
import gzip
from Bio import SeqIO

# Some variables
path = "/home/moritz/people/0023_anoxicencyclo"
ncbi = NCBITaxa()
fungi_tax_id = 4751
oomycetes = 4762

# To save taxa I checked already, the get_lineage thing from ete3.NVBITaxa is not very fast
id_dict = {}

# A function to convert a tax idea to a lineage
def get_taxa(i):
    try:
        lineage = ncbi.get_lineage(i)
    except:
        return None
    names = ncbi.get_taxid_translator(lineage)
    names = {k : v for k , v in names.items() if v not in ['root', 'cellular organisms']}

    ranks = ncbi.get_rank(lineage)
    return [names[l] for l in lineage if names.get(l)]


# a function that checks if a tax ID is a fungus or an Oomycete
def is_fungi(id):
    if id not in id_dict:
        lineage = ncbi.get_lineage(id)
        if len(lineage) > 4 and (lineage[4] == fungi_tax_id or lineage[4] == oomycetes):
            id_dict[id] = True
        else :
            id_dict[id] = False
    return id_dict[id]

# main loop
# for every file in that folder
for s in os.listdir(pjoin(path, "1000_processed_reads/")):
    print("Doing sample",s)
    reads_path = pjoin(path, "1000_processed_reads/", s, "reads")
    kaiju_path = pjoin(reads_path, "kaiju", 'nr_euk_mycocosm')
    # check to make sure it has some kaiju results and the output folder does not exist yet (e.g. I do not want to rerun stuff)
    if os.path.exists(kaiju_path) and not os.path.exists(pjoin(reads_path, "fungi")):
        # open the kaiju results file
        with open(pjoin(kaiju_path, s + ".kaiju")) as handle:
            # the actual parsing returns a dictionary with the read ID and the taxon ID it has, only for the lines that starts with "C" and "is_fungi"
            read2fungus = {l.split()[1] : l[:-1].split()[2]  for l in tqdm(handle) if l[0] == "C" and is_fungi(l[:-1].split()[2])}

        #Get full taxonomy names for all taxIDs found, I only keep the 6 first levels
        taxdict = {f : "_".join(get_taxa(f)[0:6]) for f in set(read2fungus.values())}


        #openning fwd and rev fastqs, I am using the gzip library to open them directly, it is not as fast as external gzip or pigzip, but is more lisible
        fwd_io = gzip.open(pjoin(reads_path, "fwd.fastq.gz"), "rt")
        fwd_parser = SeqIO.parse(fwd_io, "fastq")
        rev_io = gzip.open(pjoin(reads_path, "rev.fastq.gz"), "rt")
        rev_parser = SeqIO.parse(rev_io, "fastq")

        # initializing read dictionaries, I will make one file per taxa found
        fwds = {k : [] for k in set(taxdict.values())}
        revs = {k : [] for k in set(taxdict.values())}
        # looping through the fastqs
        for f, r in tqdm(zip(fwd_parser, rev_parser)):
            # just making sure I did not screw up the fastqs
            assert f.id == r.id
            # if the id is in our read2fungus dict, put it in the right read lists
            if f.id in read2fungus:
                tax = taxdict[read2fungus[f.id]]
                fwds[tax] += [f]
                revs[tax] += [r]

        # write the fastqs
        for f in set(taxdict.values()):
            os.makedirs(pjoin(reads_path, "fungi", f), exist_ok=True)
            with open(pjoin(reads_path, "fungi", f, "fwd.fastq"), "w") as handle:
                SeqIO.write(fwds[f], handle, "fastq")
            with open(pjoin(reads_path, "fungi", f, "rev.fastq"), "w") as handle:
                SeqIO.write(revs[f], handle, "fastq")


#something else (e.g. parsing my blast of the ITS sequences, not finished)
def parse_itss():
    head =['query', 'subject', 'identity', 'length', 'mismatch', 'gaps', 'qstart', 'qend', 'sstart','send','evalue','bitscore']
    raw_data = DataFrame.from_csv("Phytopthora-ID-ITS.vs_all_genomes.blastn", header = -1, index_col = False)

    with open("Phytopthora-ID-ITS.fasta") as handle:
        db = [s for s in SeqIO.parse(handle, "fasta")]
    min_len = min([len(d) for d in db])
    max_len = max([len(d) for d in db])

    cand_pheo_ITS = set(raw_data.subject)

    min_pos = {c : min( min(raw_data.loc[raw_data.subject == c].sstart), min(raw_data.loc[raw_data.subject == c].send)) for c in cand_pheo_ITS}
    max_pos = {c : max( max(raw_data.loc[raw_data.subject == c].sstart), max(raw_data.loc[raw_data.subject == c].send)) for c in cand_pheo_ITS}

    min_pos_sub = {c : min( min(raw_data.loc[raw_data.subject == c].qstart), min(raw_data.loc[raw_data.subject == c].send)) for c in cand_pheo_ITS}
    max_pos_sub = {c : max( max(raw_data.loc[raw_data.subject == c].sstart), max(raw_data.loc[raw_data.subject == c].send)) for c in cand_pheo_ITS}

    cand_pos = {c : (min(min_pos[c], max_pos[c]), max(min_pos[c], max_pos[c])) for c in cand_pheo_ITS}
    cand_pos = {k :v  for k, v in cand_pos.items() if (v[1]-v[0]) > min_len}
