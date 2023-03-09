#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
from itertools import takewhile

def merge_strainsharing_data(ttabIn, ttabOut):
    profiles_list = []
    merged_tables = None

    for f in ttabIn:
        headers = [x.strip() for x in takewhile(lambda x: x.startswith('X'), open(f))]
        names = headers[-1].split(',')
        iIn = pd.read_csv(f, sep=',', skiprows=len(headers), names=names, usecols=[0,1,3], index_col=(0,1))
        profiles_list.append(iIn)

    merged_tables = pd.concat([merged_tables, pd.concat(profiles_list, axis=1)], axis=1)
    merged_tables.to_csv(ttabOut, sep='\t')

def merge_threshold_data(suptabIn, suptabOut):
    profiles_list = []
    merged_tables = None

    for f in suptabIn:
        headers = [x.strip() for x in takewhile(lambda x: x.startswith('SGB_ID'), open(f))]
        names = headers[-1].split(',')
        iIn = pd.read_csv(f, sep=',', skiprows=1, names=names, index_col=(0))
        profiles_list.append(iIn)

    merged_tables = pd.concat([merged_tables, pd.concat(profiles_list, axis=0)], axis=0)
    merged_tables.to_csv(suptabOut, sep='\t')

argp = argparse.ArgumentParser(prog="concat_strainsharing_thresholds.py",
                               description="Performs a table join on one or more output files.")
argp.add_argument("-s", metavar="input.csv", nargs="*", help="One or more csv tables to join")
argp.add_argument("-t", metavar="input.csv", nargs="*", help="One or more csv tables consisting of nGD data to join")
argp.add_argument('-o', metavar="output.txt", help="Name of output file in which joined nGD tables are saved")
argp.add_argument('-p', metavar="output.txt", help="Name of output file in which joined supplemental tables are saved")

argp.usage = (argp.format_usage() + "\nPlease make sure to supply file paths to the files to combine.\n\n" +
              "If combining 2 files (Table_strainsharing_t__SGB4951.csv and Table_strainsharing_t__SGB13979.csv) the call should be:\n" +
              "   ./concat_files.py Table_strainsharing_t__SGB4951.csv Table_strainsharing_t__SGB13979.csv > output.txt\n\n" +
              "A wildcard to indicate all .csv files that start with Table can be used as follows:\n" +
              "    ./concat_files.py Table_strainsharing*.csv > output.txt")

def main():
    args = argp.parse_args()

    if not args.s:
        if args.t: 
            print('Concat_files: No strainsharing table files to merge!')
            print('Concat_files: Threshold data is present, merging this')
            merge_threshold_data(args.t, open(args.p,'w') if args.p else sys.stdout)
        else:
            print('No threshold tables or strainsharing tables to merge!')
        return
    
    if not args.t:
        if args.s:
            print('Concat_files: No threshold table files to merge!')
            print('Concat_files: strainsharing data is present, merging this')
            merge_strainsharing_data(args.s, open(args.o, 'w') if args.o else sys.stdout)
        return
    
    merge_strainsharing_data(args.s, open(args.o, 'w') if args.o else sys.stdout)
    merge_threshold_data(args.t, open(args.p,'w') if args.p else sys.stdout)

if __name__ == '__main__':
    main()