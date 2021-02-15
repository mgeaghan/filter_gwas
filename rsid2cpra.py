import csv
import argparse
import sys
from copy import copy

def check_args(args=None):
    parser = argparse.ArgumentParser(description="Convert dbSNP rsIDs in a GWAS summary statistics file to CHROM:POS:REF:ALT.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-g', '--gwas', help="GWAS summary statistics file. Expects a header line.", required=True)
    parser.add_argument('-s', '--sep', help="Delimiter for GWAS file. Default: tab ('\t').", default="\t")
    parser.add_argument('-c', '--chrom', help="Column number of chromosome code.", required=True)
    parser.add_argument('-b', '--bp', help="Column number of base position.", required=True)
    parser.add_argument('-a', '--a1', help="Column number of A1 allele.", required=True)
    parser.add_argument('-A', '--a2', help="Column number of A2 allele.", required=True)
    parser.add_argument('-r', '--rsid', help="Column number of rsIDs.", required=True)
    parser.add_argument('-l', '--lookup', help="Lookup file. Contains 2 tab-delimited columns: rsID and CHR:POS:REF:ALT-formatted IDs.", required=True)
    parser.add_argument('-H', '--noheader', help="Specifies that input GWAS file has no header.", action='store_true')
    parser.add_argument('-o', '--output', help="Output file.", required=True)
    return(parser.parse_args(args))

def load_lookup(lookupFile):
    lookupDict = {}
    with open(lookupFile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for line in reader:
            rsid = line[0]
            cpra = line[1]
            if rsid in lookupDict:
                lookupDict[rsid].append(cpra)
            else:
                lookupDict[rsid] = [cpra]
    return(lookupDict)

def revcomp(seq):
    rc_dict = {'A': 'T',
               'T': 'A',
               'G': 'C',
               'C': 'G',
               'a': 'T',
               't': 'A',
               'g': 'C',
               'c': 'G'}
    rc_arr = [rc_dict[c] if c in rc_dict else None for c in seq]
    if None in rc_arr:
        return None
    rc_arr.reverse()
    return ''.join(rc_arr)

def convert_gwas(gwasFile, sep, chromNum, bpNum, a1Num, a2Num, rsidNum, lookupDict, header):
    gwasList = []
    sep=sep.replace('\\t', '\t')
    try:
        chromNum = int(chromNum) - 1
        bpNum = int(bpNum) - 1
        a1Num = int(a1Num) - 1
        a2Num = int(a2Num) - 1
        rsidNum = int(rsidNum) - 1
    except:
        print("ERROR: --chrom, --bp, --a1, --a2, and --rsid expect integer arguments. Exiting.")
        sys.exit()
    with open(gwasFile, 'r') as f:
        reader = csv.reader(f, delimiter=sep)
        if header:
            gwasList.append(next(reader))
        for line in reader:
            chrom = line[chromNum]
            bp = line[bpNum]
            a1 = line[a1Num]
            a2 = line[a2Num]
            rsid = line[rsidNum]
            var_cpra = ":".join([chrom, bp, a1, a2])
            var_cpar = ":".join([chrom, bp, a2, a1])
            try:
                var_cpra_rc = ":".join([chrom, bp, revcomp(a1), revcomp(a2)])
                var_cpar_rc = ":".join([chrom, bp, revcomp(a2), revcomp(a1)])
            except TypeError:
                continue
            try:
                cpra_lookup = lookupDict[rsid]
            except KeyError:
                continue
            # if var_cpra in cpra_lookup and var_cpar in cpra_lookup:
            #     continue
            if var_cpra in cpra_lookup:
                line[rsidNum] = var_cpra
            elif var_cpar in cpra_lookup:
                line[rsidNum] = var_cpar
            elif var_cpra_rc in cpra_lookup:
                line[rsidNum] = var_cpra_rc
            elif var_cpar_rc in cpra_lookup:
                line[rsidNum] = var_cpar_rc
            else:
                continue
            gwasList.append(line)
    return(gwasList)

if __name__ == "__main__":
    args = check_args(sys.argv[1:])
    lookup = load_lookup(args.lookup)
    gwas = convert_gwas(args.gwas, args.sep, args.chrom, args.bp, args.a1, args.a2, args.rsid, lookup, (not args.noheader))
    gwas = ['\t'.join([str(i) for i in j]) for j in gwas]
    gwas = '\n'.join(gwas) + '\n'
    with open(args.output, 'w') as f:
        f.write(gwas)