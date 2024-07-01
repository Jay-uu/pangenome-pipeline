import sys, os
from Bio import Entrez
from os.path import join as pjoin
from xml.etree import ElementTree as ET
from tqdm import tqdm
import json
import pandas
from random import shuffle


def parse_exp_summary(sra_id):
    sra_summary = []
    while len(sra_summary) == 0:
        sra_summary = Entrez.read(Entrez.esummary(db="SRA", id = id))
    xxml = sra_summary[0]['ExpXml']
    md_dict = dict()
    try :
        md_dict['SRA_ID'] =ET.fromstring(sra_summary[0]['Runs']).attrib['acc']
    except :
        md_dict['SRA_ID'] =ET.fromstring("<Data>" + sra_summary[0]['Runs']+ "</Data>")[0].attrib['acc']

    dat = ET.fromstring("<Data>" + xxml + "</Data>")
    for child in dat:
        if child.tag == "Summary":
            for ll in child:
                if ll.tag == "Statistics":
                    md_dict['nb_reads'] = ll.attrib['total_spots']
                    md_dict['nb_bases'] = ll.attrib['total_bases']
        elif child.tag == "Submitter" :
            md_dict['contact_name'] = child.attrib['contact_name']
        elif child.tag == "Experiment" :
            md_dict['sample_name'] = child.attrib['name']
        elif child.tag == "Study" :
            md_dict['study'] = child.attrib['name']
        elif child.tag == "Organism" :
            md_dict['taxon'] = child.attrib.get('ScientificName')
            md_dict['taxon_id'] = child.attrib['taxid']
        elif child.tag == "Instrument" :
            md_dict['sequencer'] = ";".join(child.attrib.values())
        elif child.tag == "Library_descriptor":
            md_dict.update({cc.tag : cc.text for c in dat for cc in c if cc.tag not in [ 'LIBRARY_NAME', 'LIBRARY_LAYOUT' ] })
        elif child.tag == "Bioproject":
            md_dict['study_id'] = child.text
        elif child.tag == "Biosample":
            md_dict['sample_id'] = child.text
    md_dict['sample_attributes'] = fetch_sample_data(md_dict['sample_id'])
    return md_dict

def fetch_sample_data(sample_id):
    intern_id = Entrez.read(Entrez.esearch(db="Biosample", term = sample_id ))['IdList'][0]

    summary_dict = None
    while not summary_dict:
        try :
            summary = Entrez.read(Entrez.esummary(db="Biosample", id = intern_id))
            summary_dict = summary['DocumentSummarySet']['DocumentSummary'][0]
        except:
            continue
    samp_dict = dict()
    xxml = summary_dict['SampleData']
    dat = ET.fromstring("<Data>" + xxml + "</Data>")[0]
    atts = [l for l in dat if l.tag == "Attributes"][0]
    samp_dict.update({l.attrib['attribute_name'] : l.text  for l in atts})
    return samp_dict

#script , organism, username  = sys.argv

# all the metagenome categories in NCBI
# https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Undef&id=408169&lvl=3&keep=1&srchmode=1&unlock
organism = 'water metagenome'
username = 'fpusan@gmail.com'
all_orgs = ["aquaculture metagenome", "aquatic metagenome", "bog metagenome", "cold spring metagenome", "estuary metagenome", "freshwater metagenome", "glacier lake metagenome", "lagoon metagenome", "lake water metagenome", "oasis metagenome", "pond metagenome", "riverine metagenome"]


Entrez.email = username

query = '({orgs}) AND "strategy wgs"[Properties] AND "platform illumina"[Properties]'.format(orgs = " OR ".join(['"{organism}"[Organism]'.format(organism = a) for a in all_orgs]))
#query = '"strategy wgs"[Properties] AND "platform illumina"[Properties]'
sra_ids = Entrez.read(Entrez.esearch(db="SRA", term = query, retmax=10000000))['IdList']

#shuffle(sra_ids) #########
#sra_ids = sra_ids[:10000] #########################

sras_metadat = dict()
for id in tqdm(sra_ids):
    if id not in sras_metadat:
        try:
            sras_metadat[id] = parse_exp_summary(id)
        except:
            pass ##### ignore samples with encoding errors for now

with open("output/sra_data.json", "w") as handle:
    json.dump(sras_metadat, handle, indent=4, sort_keys=True)

for v in sras_metadat.values():
    v.update(v['sample_attributes'])
    del v['sample_attributes']
pandas.DataFrame.from_dict(sras_metadat, orient = "index").to_csv("output/sra_data.csv")
