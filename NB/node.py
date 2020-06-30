import requests
import sys
import argparse
import json

import os
import os.path
import re

ap = argparse.ArgumentParser()
ap.add_argument('-i','--input', help = 'Input file')
ap.add_argument('-n','--num', help = 'num')
ap.add_argument('-c','--cls', help = 'cluster file')
ap.add_argument('-o','--out', help = 'Output file')
args = vars(ap.parse_args())

input = ""
out = ""
num = -1
cls = "nb-net-cls.txt"

if (args['input']):
    input = args['input'].strip()
if (args['out']):
    out = args['out']
if (args['num']):
    num = int(args['num'])
if (args['cls']):
    cls = args['cls']

if not os.path.isfile(input):
    print "Can't open file {0} <br>".format(input);
    exit()
if not os.path.isfile(cls):
    print "Can't open file {0} <br>".format(cls);
    exit()

chash = {}
fp = open(cls, "r")
for line in fp:
    line = line.strip();
    l = line.split("\t")
    chash[l[0]] = l[2:]
fp.close();

obj = []
with open(input) as f:
    obj = json.load(f)

def uniq(mylist):
    used = set()
    unique = [x for x in mylist if x not in used and (used.add(x) or True)]
    return unique

nodes = []
for i in obj:
    if (i[0] is num and len(i) > 2):
        #print "\t".join([k[0] for k in i[2]])
        for k in i[2]:
            if (k[0] in chash):
                nodes += chash[k[0]]
    if (i[0] is num and len(i) > 3):
        for k in i[3]:
            l = k.split(" ");
            if (l[0] in chash):
                nodes += chash[l[0]]
            else:
                nodes += [l[0]]

if out == "":
    print " ".join(uniq(nodes))
else:
    fp = open(out, "w")
    fp.write(" ".join(uniq(nodes)))
    fp.close()
