#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np
from biom import load_table, Table
import itertools

###
#   Script to estimate core, accessory and unique OTUs of treatments in a given biom file.
###
#
#   COPYRIGHT DISCALIMER:
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#
#   Author: Tim Kahlke, tim.kahlke@audiotax.is
#   Date:   May 2017
#



def main(args):
    input_table= load_table(args.biom)
    
    treats = []
    if args.treatments:
        treats = _getT(args.treatments)
    else:
        treats = np.ndarray.tolist(input_table.ids())

    # get treatment indices
    inds = [np.ndarray.tolist(input_table.ids()).index(x) for x in treats]

    core = "test"
    for i in inds:
        tc = [y for (x,y) in enumerate(input_table.ids(axis="observation")) if input_table[x,i] >= args.abundance_minimum]
        if core != "test":
            core = [x for x in tc if x in core]
        else:
            core = tc
    output = open(args.output,"w")
    output.write("# Core OTUs of file %s: %d \n" % (args.biom,len(core)))
    
    for c in core:
        output.write("%s\t%s\n" %(c,";".join(input_table.metadata(axis="observation")[int(np.where(input_table.ids(axis="observation")==c)[0])]['taxonomy'])))

    output.close()








def _curve_val(table):
    co = 0
    for i in range(1,len(table.ids(axis="observation"))-1):
        print("\nCreating combinations of %d\n" % (i))
        combs = itertools.combinations(range(0,len(table.ids(axis="observation"))-1),i)
        for c in combs:
            co+=1
    print(str(co))


#    print(input_table.metadata(axis="observation")[-6]['taxonomy'])
#    print([method for method in dir(input_table) if callable(getattr(input_table, method))])
  


def _getT(fp):
    td =[] 
    with open(fp,"r") as f:
        for line in f:
            td.append(line.replace("\n",""))
    return td



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Estimation of core, accessory and unique OTU numbers")
    parser.add_argument("biom", help="OTU table in biom format")
    parser.add_argument("output", help="Directory for output files")
    parser.add_argument("-t","--treatments", help="File of treatments that should be considered. Otherwise all samples/treatments are used for the analysis")
    parser.add_argument("-m", "--abundance_minimum", help="Abundance minimum. If set only OTUs with given relative abundance are considered", type=float, default=0.0)

    args = parser.parse_args()
    main(args)

