"""
Support for reading and writing genomic intervals from delimited text files.
"""

import sys
from itertools import *
from bx.tabular.io import *
from bx.bitset import *

class MissingFieldError( ParseError ):
    pass

class FieldFormatError( ParseError ):
    def __init__( self, *args, **kwargs):
        ParseError.__init__( self, *args, **kwargs )
        self.expected = kwargs.get("expected",None)
    def __str__( self ):
        if self.expected:
            return ParseError.__str__( self ) + ", " + self.expected + " expected"
        else:
            return ParseError.__str__( self )

class StrandFormatError( ParseError ):
    pass

class GenomicInterval( TableRow ):
    """
    A genomic interval stored in a set of fields (a row of a table)
    """
    def __init__( self, reader, fields, chrom_col, start_col, end_col, strand_col, default_strand, fix_strand=False ):
        TableRow.__init__( self, reader, fields )
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.nfields = nfields = len( fields )
        # Parse chrom/source column
        if chrom_col >= nfields:
            raise MissingFieldError( "No field for chrom_col (%d)" % chrom_col )
        self.chrom = fields[chrom_col]
        # Parse start column and ensure it is an integer
        if start_col >= nfields:
            raise MissingFieldError( "No field for start_col (%d)" % start_col )
        try:
            self.start = int( fields[start_col] )
        except ValueError, e:
            raise FieldFormatError( "Could not parse start_col: " + str( e ), expected="integer" )
        # Parse end column and ensure it is an integer
        if end_col >= nfields:
            raise MissingFieldError( "No field for end_col (%d)" % end_col )
        try:
            self.end = int( fields[end_col] )
        except ValueError, e:
            raise FieldFormatError( "Could not parse end_col: " + str( e ), expected="integer" )
        # Ensure start <= end
        if self.end < self.start:
            raise ParseError( "Start is greater than End. Interval length is < 1." )
        # Parse strand and ensure it is valid
        if strand_col >= nfields or strand_col < 0:
            # This should probable be immutable since the fields are 
            # not updated when it is set
            self.strand = default_strand
        else:
            strand = fields[strand_col]
            if strand not in ( "+", "-", "."):
                if fix_strand:
                    strand = "+"
                else: raise StrandFormatError( "Strand must be either '+' or '-'" )
            self.strand = strand
    def __setattr__( self, name, value ):
        if name == "chrom":
            self.fields[self.chrom_col] = str( value )
        elif name == "start":
            self.fields[self.start_col] = str( value )
    'chr1\\tfoo\\t30\\t100\\txxx'
    >>> assert type( elements[2] ) is GenomicInterval
    >>> assert type( elements[3] ) is Comment
    >>> assert type( elements[4] ) is GenomicInterval
    """
    def __init__( self, input, chrom_col=0, start_col=1, end_col=2, strand_col=5, 
                  default_strand=".", return_header=True, return_comments=True, force_header=None, fix_strand=False, comment_lines_startswith = ["#", "track "] ):
        TableReader.__init__( self, input, return_header, return_comments, force_header, comment_lines_startswith )
        self.chrom_col = chrom_col
        self.start_col = start_col
        self.end_col = end_col
        self.strand_col = strand_col
        self.default_strand = default_strand
        self.fix_strand = fix_strand
    def parse_row( self, line ):
        return GenomicInterval( self, line.split( "\t" ), self.chrom_col, 
                                self.start_col, self.end_col,
                                self.strand_col, self.default_strand, fix_strand=self.fix_strand )

    def binned_bitsets( self , upstream_pad=0, downstream_pad=0, lens={} ):
        # The incoming lens dictionary is a dictionary of chromosome lengths
        # which are used to initialize the bitsets.
        last_chrom = None
        last_bitset = None
        bitsets = dict()
        for interval in self:
            if type( interval ) == GenomicInterval:
                chrom = interval[self.chrom_col]
                if chrom != last_chrom:
                    if chrom not in bitsets:
                        size = lens.get( chrom, MAX )
                        try:
                            bbs = BinnedBitSet( size )
        # It is assumed that the reader is an interval reader, i.e. it has chr_col, start_col, end_col and strand_col attributes. 
        NiceReaderWrapper.__init__( self, reader.input, chrom_col=reader.chrom_col, start_col=reader.start_col, end_col=reader.end_col, strand_col=reader.strand_col)
        self.lens = lens
    def next( self ):
        while True:
            rval = NiceReaderWrapper.next( self )
            if type( rval ) == GenomicInterval and rval.end > self.lens.get( rval.chrom, MAX ): # MAX_INT is defined in bx.bitset
                try:
                    # This will only work if reader is a NiceReaderWrapper
                    self.skipped += 1
                    # no reason to stuff an entire bad file into memmory
                    if self.skipped < 10:
                        self.skipped_lines.append( ( self.linenum, self.current_line, str( e ) ) )
