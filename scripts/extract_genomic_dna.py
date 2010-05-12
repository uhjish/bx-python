#!/usr/bin/env python
"""
usage: %prog $input $out_file1
    -1, --cols=N,N,N,N: Columns for start, end, strand in input file
    -o, --output_format=N: the data type of the output file
    -s, --seq_path=N: the directory containing the chromosome fasta files
    -l, --left_flank=N: extra bases on the  left
    -r, --right_flank=N: extra bases on the  right
"""
import sys, string, os, re
from bx.cookbook import doc_optparse
import bx.seq.nib
import bx.seq.twobit

assert sys.version_info[:2] >= ( 2, 4 )
    
def stop_err( msg ):
    sys.stderr.write( msg )
    sys.exit()

# Default chrom, start, end, strand cols for a bed file
BED_DEFAULT_COLS = 0, 1, 2, 5

def parse_cols_arg( cols ):
    """Parse a columns command line argument into a four-tuple"""
    if cols:
        # Handle case where no strand column included - in this case, cols
        # looks something like 1,2,3,
        if cols.endswith( ',' ):
            cols += '0'
        col_list = map( lambda x: int( x ) - 1, cols.split(",") )
        return col_list
    else:
        return BED_DEFAULT_COLS

def reverse_complement( s ):
    complement_dna = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "N":"N", "n":"n" }
    reversed_s = []
    for i in s:
        reversed_s.append( complement_dna[i] )
    reversed_s.reverse()
    return "".join( reversed_s )

def __main__():

    lflank = 0
    rflank = 0

    options, args = doc_optparse.parse( __doc__ )
    try:
        chrom_col, start_col, end_col, strand_col = parse_cols_arg( options.cols )
        output_format = options.output_format
        seq_path = options.seq_path
        if ( options.left_flank): lflank = int(options.left_flank)
        if ( options.right_flank): rflank = int( options.right_flank)
        input_filename, output_filename = args
    except:
        doc_optparse.exception()
    includes_strand_col = strand_col >= 0
    strand = None
    nibs = {}
    twobits = {}
    if not os.path.exists( seq_path ):
        # If this occurs, we need to fix the metadata validator.
        print "No sequences are available for '%s', request them by reporting this error."

    skipped_lines = 0
    first_invalid_line = 0
    invalid_line = ''
    fout = open( output_filename, "w" )
    warnings = []
    warning = ''
    twobitfile = None
    dbkey=seq_path
     
    for i, line in enumerate( open( input_filename ) ):
        line = line.rstrip( '\r\n' )
        if line and not line.startswith( "#" ):
            fields = line.split( '\t' )
            try:
                chrom = fields[chrom_col]
                ostart = int( fields[start_col] )
                oend = int( fields[end_col] )
                start = ostart - lflank
                end = oend + rflank
                if includes_strand_col:
                    strand = fields[strand_col]
            except:
                warning = "Invalid chrom, start or end column values. "
                warnings.append( warning )
                skipped_lines += 1
                if not invalid_line:
                    first_invalid_line = i + 1
                    invalid_line = line
                continue
            if start > end:
                warning = "Invalid interval, start '%d' > end '%d'.  " % ( start, end )
                warnings.append( warning )
                skipped_lines += 1
                if not invalid_line:
                    first_invalid_line = i + 1
                    invalid_line = line
                continue

            if strand not in ['+', '-']:
                strand = '+'
            sequence = ''

            if seq_path and os.path.exists( "%s/%s.nib" % ( seq_path, chrom ) ):
                if chrom in nibs:
                    nib = nibs[chrom]
                else:
                    nibs[chrom] = nib = bx.seq.nib.NibFile( file( "%s/%s.nib" % ( seq_path, chrom ) ) )
                try:
                    sequence = nib.get( start, end-start )
                except:
                    warning = "Unable to fetch the sequence from '%d' to '%d' for build '%s'. " %( start, end-start, dbkey )
                    warnings.append( warning )
                    skipped_lines += 1
                    if not invalid_line:
                        first_invalid_line = i + 1
                        invalid_line = line
                    continue
            elif seq_path and os.path.isfile( seq_path ):
                if not(twobitfile):
                    twobitfile = bx.seq.twobit.TwoBitFile( file( seq_path ) )
                try:
                    sequence = twobitfile[chrom][start:end]
                except:
                    warning = "Unable to fetch the sequence from '%d' to '%d' for build '%s'. " %( start, end-start, dbkey )
                    warnings.append( warning )
                    skipped_lines += 1
                    if not invalid_line:
                        first_invalid_line = i + 1
                        invalid_line = line
                    continue
            else:
                warning = "Chromosome by name '%s' was not found for build '%s'. " % ( chrom, dbkey )
                warnings.append( warning )
                skipped_lines += 1
                if not invalid_line:
                    first_invalid_line = i + 1
                    invalid_line = line
                continue
            if sequence == '':
                warning = "Chrom: '%s', start: '%s', end: '%s' is either invalid or not present in build '%s'. " %( chrom, start, end, dbkey )
                warnings.append( warning )
                skipped_lines += 1
                if not invalid_line:
                    first_invalid_line = i + 1
                    invalid_line = line
                continue
            if includes_strand_col and strand == "-":
                sequence = reverse_complement( sequence )
            sequence = sequence[0:lflank].lower() + sequence[lflank:len(sequence)-rflank+1].upper() + sequence[len(sequence)-rflank+1:len(sequence)].lower()

            if output_format == "fasta" :
                l = len( sequence )        
                c = 0
                fields = [dbkey, str( chrom ), str( ostart ), str( oend ), strand]
                meta_data = "_".join( fields )
                fout.write( ">%s\n" % meta_data )
                while c < l:
                    b = min( c + 50, l )
                    fout.write( "%s\n" % str( sequence[c:b] ) )
                    c = b
            else: # output_format == "interval"
                meta_data = "\t".join( fields )
                fout.write( "%s\t%s\n" % ( meta_data, str( sequence ) ) )

    fout.close()

    if warnings:
        warn_msg = "%d warnings, 1st is: " % len( warnings )
        warn_msg += warnings[0]
        print warn_msg
    if skipped_lines:
        print 'Skipped %d invalid lines, 1st is #%d, "%s"' % ( skipped_lines, first_invalid_line, invalid_line )

if __name__ == "__main__": __main__()
