#!/usr/bin/env python2.3

"""
NOT YET CORRECT!!!

usage: %prog start end maf_files...
"""

import psyco_full

import cookbook.doc_optparse

import bx.align.maf
import bx.align as align
from bx import misc
import bx.seq.nib
import os
import string
import sys

tree_tx = string.maketrans( "(),", "   " )

def main():

    sources = sys.argv[1].translate( tree_tx ).split()
    seq_db = load_seq_db( sys.argv[2] )
    index = bx.align.maf.MultiIndexed( sys.argv[3:] )

    out = bx.align.maf.Writer( sys.stdout )

    for line in sys.stdin:
        ref_src, start, end = line.split()[0:3]
        do_interval( sources, index, out, ref_src, int( start ), int( end ), seq_db )

    out.close()

def load_seq_db( fname ):
    db = {}
    for line in open( fname ):
        fields = line.split(',')
        src = fields[1] + "." + fields[2]
        seq = fields[4]
        db[src]=seq.strip()
    return db

def do_interval( sources, index, out, ref_src, start, end, seq_db ):

    assert sources[0].split('.')[0] == ref_src.split('.')[0], "%s != %s" % ( sources[0].split('.')[0], ref_src.split('.')[0] )

    base_len = end - start
        
    blocks = index.get( ref_src, start, end )
    # From low to high score
    blocks.sort( lambda a, b: cmp( a.score, b.score ) )

    mask = [ -1 ] * base_len

    # print len( blocks )
    # print blocks[0]
    
    ref_src_size = None
    for i, block in enumerate( blocks ):
        ref = block.get_component_by_src_start( ref_src )
        ref_src_size = ref.src_size
        assert ref.strand == "+"
        slice_start = max( start, ref.start )
        slice_end = min( end, ref.end )
        for j in range( slice_start, slice_end ):
            mask[j-start] = i

    #print mask

    tiled = []
    for i in range( len( sources ) ): tiled.append( [] )

    for ss, ee, index in intervals_from_mask( mask ):
        if index < 0:
            tiled[0].append( bx.seq.nib.NibFile( open( seq_db[ ref_src ] ) ).get( start+ss, ee-ss ) )
            for row in tiled[1:]:
                row.append( "-" * ( ee - ss ) )
        else:
            slice_start = start + ss
            slice_end = start + ee
            block = blocks[index]
            ref = block.get_component_by_src_start( ref_src )
            sliced = block.slice_by_component( ref, slice_start, slice_end ) 
            sliced = sliced.limit_to_species( sources )
            sliced.remove_all_gap_columns()
            for i, src in enumerate( sources ):
                comp = sliced.get_component_by_src_start( src )
                if comp:
                    tiled[i].append( comp.text )
                else:
                    tiled[i].append( "-" * sliced.text_size )
        
    a = align.Alignment()
    for i, name in enumerate( sources ):
        text = "".join( tiled[i] )
        size = len( text ) - text.count( "-" )
        if i == 0:
            if ref_src_size is None: ref_src_size = bx.seq.nib.NibFile( open( seq_db[ ref_src ] ) ).length
            c = align.Component( ref_src, start, end-start, "+", ref_src_size, text )
        else:
            c = align.Component( name + ".fake", 0, size, "?", size, text )
        a.add_component( c )

    out.write( a )


def intervals_from_mask( mask ):
    start = 0
    last = mask[0]
    for i in range( 1, len( mask ) ):
        if mask[i] != last:
            yield start, i, last
            start = i
            last = mask[i]
    yield start, len(mask), last

main()
