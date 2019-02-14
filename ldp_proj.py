#!/usr/bin/env python3
#
# ldp projections
# 

# import system modules
import argparse

# import local modules
import prj_parser as pparse

# parse command line
parser = argparse.ArgumentParser()
parser.add_argument("file", help=".prj file to load")
args = parser.parse_args()

srs_dict = pparse.prj_reader(args.file)
if srs_dict == 0:
	print('Error reading the prj file, see error above!')
	exit(1)

pparse.prj_writer('test.prj', srs_dict)
