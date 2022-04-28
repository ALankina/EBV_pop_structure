# Script to scrape genbank for geographical & clinical metadata of CMV sequences
# current form gets data from each sequence in a fasta alignment, could provide an id_list using Entrez db search etc.
import os
import re
from Bio import SeqIO, Entrez, Seq, SeqRecord
import pycountry_convert as pc

Entrez.email = "oscar.charles.18@ucl.ac.uk"
accession_list_file = "data/ebv_genbank_over_150kb.tab"


accessions = []
with open(accession_list_file) as f:
    for line in f:
        accessions.append(line.strip("\n"))

print("accessions\tcontinent\tcountry\tdate\tisolate\tstrain")
for accession in accessions:
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype = "gb")
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

    print(accession+"\t"+continent+"\t"+country+"\t"+date+"\t"+isolate+"\t"+strain)
