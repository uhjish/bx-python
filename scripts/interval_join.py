#!/usr/bin/env python2.4

"""
Match up intersecting intervals from two files. This performs a "full join", 
any pair of intervals with any basewise overlap will be printed side-by-side.

usage: %prog bed1 bed2
"""

from __future__ import division

import psyco_full

import string
import sys

from optparse import OptionParser

import bx.intervals.io
import bx.intervals.intersection
import bx.intervals.operations.join

#def join(leftSet, rightSet, mincols=1, leftfill=True, rightfill=True, asfraction=False, matchStrand=STRAND_NEUTRAL, outColumns=[-1,-1]):

def main():
    usage = "usage: %prog [options] intervalFile1 intervalFile2"

    parser = OptionParser()
    parser.add_option("-f", "--fractional",
                        action="store_true", dest="fraction", default=False,
                        help="fractional intersections -- fraction of the smaller locus")
    parser.add_option("-o", "--overlap", type="float", dest="overlap", default=1.0,
                        help="minimum overlapping bases for intersection (or minimum fraction)")
    parser.add_option("-c", "--columns", dest="columns", default="0,0",
                        help="columns to output, 0 for all columns (ex: 4,0 gives c4 from file1 and all cols from file2)")
    parser.add_option("-s", "--strand", type="int", dest="strand", default=0,
                        help="strand matching: 0 - neutral, 2 - strict match, 1 - permissive match, -2 - strict complement, -1 - permissive complement")
    
    (options, args) = parser.parse_args()
    columns = map(int, options.columns.split(","))
    columns[0] -= 1
    columns[1] -= 1
    if (len(columns) != 2):
        sys.exit("Error! Columns must be defined as two integers -- ex: 4,4")
    if (options.fraction):
        if (options.overlap <= 0 or options.overlap > 1):
            sys.exit("Error! For fractional intersections, overlap must be in range (0,1]")
    if (options.strand >2 or options.strand < -2):
        sys.exit( "Error! Strand must be an integer in range [-2,2]" )
    
    
    intersecters = {}

    set1 = bx.intervals.io.GenomicIntervalReader( open(args[0]) )
    set2 = bx.intervals.io.GenomicIntervalReader( open(args[1]) )
    
    result = bx.intervals.operations.join.join(set1, set2, mincols=options.overlap, leftfill=False, rightfill=False, asfraction=options.fraction, matchStrand=options.strand, outColumns=columns)

    for pair in result:
        print "\t".join(pair)

if __name__ == "__main__":
    main()
