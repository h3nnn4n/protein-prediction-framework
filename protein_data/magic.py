#!/usr/bin/python3.5

import subprocess
import requests
import gzip
import sys
import os

args = sys.argv[1:]
if len(args) < 1:
    print("Missing arguments")


def fetch(target):
    if os.path.exists(target):
        print("Target exists, skipping...")
        return

    print("Summoning %s" % target)

    os.makedirs(target)
    os.chdir(target)

    os.makedirs('output')
    os.makedirs('psspred')

    pdbid = target.upper()
    url = "http://files.rcsb.org/download/%s.pdb.gz" % pdbid
    r = requests.get(url)

    with open(target + '.pdb.gz', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    with gzip.open(target + '.pdb.gz', 'rb') as f_in, open(target + '.pdb', 'wb') as f_out:
        f_out.write(f_in.read())

    url = "http://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=%s&compressionType=uncompressed" % pdbid
    r = requests.get(url)

    with open(target + '.fasta', 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)

    pss_path = "/home/h3nnn4n/PSSpred/PSSpred.pl"

    if not os.path.isfile(pss_path):
        pss_path = "/home/h3nnn4n/rosetta/PSSpred/PSSpred.pl"

    os.chdir('psspred')

    subprocess.call([pss_path, '../' + target + '.fasta'])

    os.chdir("../")

    with open('psspred/seq.dat.ss', 'rt') as f_in, open(target + '.psipred.ss2', 'wt') as f_out:
        c = 0
        for line in f_in.readlines():
            if c == 0:
                c += 1
                f_out.write('# PSIPRED VFORMAT (PSIPRED V2.6 by David Jones)\n')
            else:
                f_out.write(line)

    frag_path = '../../tools/frag_picker/frag_picker.sh'

    subprocess.call([frag_path, target])

    os.chdir("../")


for target in args:
    fetch(target)
