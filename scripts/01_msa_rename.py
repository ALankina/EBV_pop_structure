# Script to scrape genbank for geographical & clinical metadata of CMV sequences
# current form gets data from each sequence in a fasta alignment, could provide an id_list using Entrez db search etc.
import os
import re
from Bio import SeqIO, Entrez, Seq, SeqRecord
import pycountry_convert as pc

Entrez.email = "oscar.charles.18@ucl.ac.uk"
# fasta file of minimum set
fasta_file = 'data/ebv_150k+_subset.fa'
regex= '^([^.]+)' # everything before fullstop
#print("id\tcontinent\tcountry\tdate\tisolate\tstrain")
start = 1
i = 0 
for record in SeqIO.parse(fasta_file, "fasta"):
    i = i+1
    if(i < start):
        continue
    # get genbank id
    #id = "KP745672"
    id = record.id
    id = re.search(regex, id)
    id = id.group()

    # get loc
    handle = Entrez.efetch(db="nucleotide", id=id, rettype = "gb")
    seq_record = SeqIO.read(handle, "gb")
    keys = seq_record.features[0].qualifiers
    country = keys.get('country')
    if country is None:
        country = "NA"
        continent = "NA"
    else:
        country = ''.join(country)
        # country -> continent
        country_alpha2 = pc.country_name_to_country_alpha2(country)
        country_continent_code = pc.country_alpha2_to_continent_code(country_alpha2)
        continent = pc.convert_continent_code_to_continent_name(country_continent_code)



    # get date
    date = keys.get('collection_date')
    if date is None:
        date = "NA"
    else:
        date = ''.join(date)


    # isolate, if exists then is clinical not passage
    isolate = keys.get('isolation_source')
    if isolate is None:
        isolate = "NA"
    else:
        isolate = ''.join(isolate)

    #strain - lab strains have this value
    strain = keys.get('strain')
    if strain is None:
        strain = "NA"
    else:
        strain = ''.join(strain)

    #print(id+"\t"+continent+"\t"+country+"\t"+date+"\t"+isolate+"\t"+strain)




    # now rename the fasta headers
    newname = continent+"_"+country+"_"+record.name
    newname = newname.replace(" ","-") # some aligners fall over with spaces
    print(">"+newname)
    print(record.seq)
