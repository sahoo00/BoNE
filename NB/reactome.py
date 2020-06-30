import requests
import sys
import argparse
import json

import os
import os.path
import re

ap = argparse.ArgumentParser()
ap.add_argument('-i','--input', help = 'Input file')
ap.add_argument('-o','--out', help = 'Output file')
args = vars(ap.parse_args())

input = ""
out = ""

idhash = {}
cfile = "zhang-2014-rnaseq-nb-idx.txt";
if not os.path.isfile(cfile):
    print("Can't open file {0} <br>".format(cfile))
    exit()
fp = open(cfile, "r")
for line in fp:
    line = line.strip();
    ll = re.split("\t", line);
    idhash[ll[0]] = ll[2]
fp.close();

if (args['input']):
    input = args['input'].strip()
if (args['out']):
    out = args['out']

if not os.path.isfile(input):
    print "Can't open file {0} <br>".format(input);
    exit()
fp = open(input, "r")
lines = fp.readlines()
fp.close();

d = "\t".join(lines)
d = re.sub("[\r\n]", "", d)
d = " ".join([idhash[k] for k in re.split("\s", d) if k in idhash])

reactomeURI = 'http://www.reactome.org/AnalysisService/identifiers/projection?pageSize=100&page=1';
response = requests.post(reactomeURI, data = d, \
    headers = { "Content-Type": "text/plain",  "dataType" : "json" })
obj = json.loads(response.text)
print "\t".join(["name", "pValue", "fdr"])
for pathway in obj["pathways"]:
    print "\t".join([str(i) for i in [pathway["name"], pathway["entities"]["pValue"], pathway["entities"]["fdr"]]])
