#!/usr/bin/env python

import argparse
import os
import sys
import numpy as np
from biom import load_table, Table
import itertools
import logging

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
    logging.basicConfig(format='',level=logging.INFO)
    logger = logging.getLogger()


    # If treatmenst are given get treatment names
    # Additionally, if reps in treatment file get names of reps
    # else reps are empty and treats = whole data set
    if args.treatments:
        (treats,reps) = _getT(args.treatments)
    else:
        treats = np.ndarray.tolist(input_table.ids())
        reps = ()
    # get treatment indices/data set table indices
    inds = _get_inds(input_table,treats)

    # get abundance threshold
    sums = _get_sample_sums(input_table,inds)

    if args.replicate_threshold:
        core = _get_rep_core_otus(input_table,sums,reps,args.replicate_threshold,inds,args.abundance_minimum,args.exclude)
    else:
        core = _get_all_core_otus(input_table,sums,inds,args.abundance_minimum)

    logger.info("\n[STATUS] Done! No of core OTUs: %d\n\n" % (len(core)))

    output = open(args.output,"w")
    output.write("# Core OTUs of file %s: %d \n" % (args.biom,len(core)))
    
    for c in core:
        if args.print_taxonomy == "True": 
            output.write("%s\t%s\n" %(c,";".join(input_table.metadata(axis="observation")[int(np.where(input_table.ids(axis="observation")==c)[0])]['taxonomy'])))
        else:
            output.write("%s\n" % (c))

    output.close()


# get core wrt replicates/groups threshold
def _get_rep_core_otus(biom,sums,reps,rt,all_inds,ab_min,rta):
    core = -1 
    for r in reps:
        rep_core = []
        inds = _get_inds(biom,reps[r])
        min_samples = 0
        if rta == 'r':
            min_samples = rt
        else:
            min_samples = len(inds) - rt
        for (x,y) in enumerate(biom.ids(axis="observation")):
            tmp = [i for i in inds if biom[x,i] > sums[all_inds.index(i)] * ab_min]
            if len(tmp) >= min_samples:
                rep_core.append(y)

        if core == -1:
            core = rep_core
        else:
            core = [x for x in rep_core if x in core]
    return core 



# Get core for all treatments (regardless of replicates/groups)
def _get_all_core_otus(biom,sums,inds,ab_min):
    core = -1
    for i in inds:
        tc = [y for (x,y) in enumerate(biom.ids(axis="observation")) if biom[x,i] and biom[x,i] >= sums[inds.index(i)] * ab_min]
        if core == -1:
            core = tc
        else:
            core = [x for x in tc if x in core]

    return core


# return indices of given sample names
def _get_inds(biom,names):
    inds = [np.ndarray.tolist(biom.ids()).index(x) for x in names]
    return inds


# calculate thresholds based on either complete
# data set or treatments
def _get_sample_sums(biom,tind):

    # get array of sum per sample/treatments

    all_cs = biom.sum(axis="sample")
    return all_cs[tind[:None],]


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
  

# Get names of samples and, if given, replicates
def _getT(fp):
    td =[] 
    reps = {}
    with open(fp,"r") as f:
        for line in f:
            tabs = line.replace("\n","").replace("\r","").split("\t")
            if len(tabs) == 2:
                if tabs[1] in reps:
                    reps[tabs[1]].append(tabs[0])
                else:
                    reps[tabs[1]] = [tabs[0]]
                td.append(tabs[0])
            elif len(tabs) == 1:
                td.append(tabs[0])
    return (td,reps)



if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Estimation of core, accessory and unique OTU numbers")
    parser.add_argument("biom", help="OTU table in biom format")
    parser.add_argument("output", help="Directory for output files")
    parser.add_argument("-t","--treatments", help="File of treatments that should be considered. Otherwise all samples/treatments are used for the analysis")
    parser.add_argument("-m", "--abundance_minimum", help="Abundance minimum. If set only OTUs with given relative abundance are considered", type=float, default=0.0)
    parser.add_argument("-r", "--replicate_threshold", help="If set at least the given number of replicates/samples of the same group have to make the given abundance threshold",type=int)
    parser.add_argument("-p", "--print_taxonomy", help="If set to False no taxonomy will be printed to result file (default=True)", default="True")
    parser.add_argument("-x", "--exclude", help="If set to x the replicate threshold will be interpreted as the maximum outliers that can be excluded per treatment.", default="r", choices=['r','x'])

    args = parser.parse_args()
    main(args)

