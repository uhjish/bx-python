#!/usr/bin/env python2.4

"""
usage: %prog score_file interval_file [out_file] [options] 
    -b, --binned: 'score_file' is actually a directory of binned array files
"""

from __future__ import division

import psyco_full
import sys
import os, os.path
from UserDict import DictMixin
import bx.wiggle
from bx.binned_array import BinnedArray, FileBinnedArray
from fpconst import isNaN
import cookbook.doc_optparse
import misc

class FileBinnedArrayDir( DictMixin ):
    """
    Adapter that makes a directory of FileBinnedArray files look like
    a regular dict of BinnedArray objects. 
    """
    def __init__( self, dir ):
        self.dir = dir
        self.cache = dict()
    def __getitem__( self, key ):
        value = None
        if key in self.cache:
            value = self.cache[key]
        else:
            fname = os.path.join( self.dir, "%s.ba" % key )
            if os.path.exists( fname ):
                value = FileBinnedArray( open( fname ) )
                self.cache[key] = value
        if value is None:
            raise KeyError( "File does not exist: " + fname )
        return value

def load_scores_wiggle( fname ):
    """
    Read a wiggle file and return a dict of BinnedArray objects keyed 
    by chromosome.
    """
    scores_by_chrom = dict()
    for chrom, pos, val in bx.wiggle.Reader( misc.open_compressed( fname ) ):
        if chrom not in scores_by_chrom:
            scores_by_chrom[chrom] = BinnedArray()
        scores_by_chrom[chrom][pos] = val
    return scores_by_chrom

def load_scores_ba_dir( dir ):
    """
    Return a dict-like object (keyed by chromosome) that returns 
    FileBinnedArray objects created from "key.ba" files in `dir`
    """
    return FileBinnedArrayDir( dir )
    
def main():

    # Parse command line
    options, args = cookbook.doc_optparse.parse( __doc__ )
    try:
        score_fname = args[0]
        interval_fname = args[1]
        if len( args ) > 2:
            out_file = open( args[2], 'w' )
        else:
            out_file = sys.stdout
        binned = bool( options.binned )
    except:
        cookbook.doc_optparse.exit()

    if binned:
        scores_by_chrom = load_scores_ba_dir( score_fname )
    else:
        scores_by_chrom = load_scores_wiggle( score_fname )

    for line in open( interval_fname ):
        fields = line.split()
        chrom, start, stop = fields[0], int( fields[1] ), int( fields[2] )
        total = 0
        count = 0
        min_score = 100000000
        max_score = -100000000
        for i in range( start, stop ):
            if chrom in scores_by_chrom and scores_by_chrom[chrom][i]:
                score = scores_by_chrom[chrom][i]
                if not isNaN( score ):
                    total += score
                    count += 1
                    max_score = max( score, max_score )
                    min_score = min( score, min_score )
        if count > 0:
            avg = total/count
        else:
            avg = "nan"
            min_score = "nan"
            max_score = "nan"
            
        print >> out_file, chrom, start, stop, avg, min_score, max_score

    out_file.close()

if __name__ == "__main__": main()
