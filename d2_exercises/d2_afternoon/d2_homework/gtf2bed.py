#!/usr/bin/env python3

import sys

my_file = open( sys.argv[1])

for my_line in my_file:
    if "##" in my_line:
        continue
    tabs = my_line.split( "\t" )

    s = tabs[8].split(";")[2].lstrip('gene_name "').rstrip( '"' )

    print(tabs[0], tabs[3], tabs[4], s)

my_file.close()