
"""
Match up intersecting intervals from two files. This performs a "full join", 
any pair of intervals with any basewise overlap will be printed side-by-side.

usage: %prog [options] queryFile wordFile1 wordFile2 wordFile3...
"""

from __future__ import division

import psyco_full

import string
import sys
import threading

from optparse import OptionParser

import bx.intervals.io
import bx.intervals.intersection
import bx.intervals.operations.rightwords

#def getpairs(leftSet, rightSet, rightCol, mincols=1, asfraction=False, matchStrand=STRAND_NEUTRAL, skipChrNames=True, skipStrandNames=True):

def main():
    usage = "usage: %prog [options] queryFile wordFile1 wordFile2 wordFile3..."

    parser = OptionParser()
    parser.add_option("-b", "--background",
                        dest="background", default=None,
                        help="intervals to use for estimating backgrund frequency of associations")
    parser.add_option("-p", "--pairs",
                        action="store_true", dest="pairs", default=False,
                        help="output individual pairs instead of summary counts for each query word")
    parser.add_option("-f", "--fractional",
                        action="store_true", dest="fraction", default=False,
                        help="fractional intersections -- fraction of the smaller locus")
    parser.add_option("-q", "--querycolumn", type="int", dest="querycol", default=4,
                        help="column in the query interval set to associate the mined words with (default=4 sensible for BED files)")
    parser.add_option("-o", "--overlap", type="float", dest="overlap", default=1.0,
                        help="minimum overlapping bases for intersection (or minimum fraction)")
    parser.add_option("-s", "--strand", type="int", dest="strand", default=0,
                        help="strand matching: 0 - neutral, 2 - strict match, 1 - permissive match, -2 - strict complement, -1 - permissive complement")
    
    (options, args) = parser.parse_args()

    if (options.fraction):
        if (options.overlap <= 0 or options.overlap > 1):
            sys.exit("Error! For fractional intersections, overlap must be in range (0,1]")

    if (options.strand >2 or options.strand < -2):
        sys.exit( "Error! Strand must be an integer in range [-2,2]" )
        
    if (options.pairs and options.background !=None):
        sys.exit( "Error! Cannot have a background file when outputting pairs")
    
    options.querycol = options.querycol-1

    queryFile = args[0]
    wordFileArgs = {}
    for argdex in range(1, len(args)):
        argStrs = args[argdex].split(",")
        curArgs = [options.fraction, options.overlap, options.strand]
        if len(argStrs) ==1:
            wordFileArgs[ argStrs[0] ] = curArgs
            continue
        if (string.upper(argStrs[1][0]) == 'F'):
            curArgs[0] = True
        else:
            curArgs[0] = False
        try:
            curArgs[1] = float(argStrs[2])
        except:
            print "in override parameters for file [" +str(argStrs)+ "] -- the third argument must be numeric!"
        try:
            curArgs[2] = int(argStrs[3])
        except:
            print "in override parameters for file [" +str(argStrs)+ "] -- the fourth argument must be an integer in [-2,2]!"
        wordFileArgs[ argStrs[0] ] = curArgs

    result = []
    qcounters = []
    bcounters = []
    for wfile, wargs in wordFileArgs.items():
        qset = makeFormatSafeIntervalReader( queryFile )
        wset = makeFormatSafeIntervalReader( wfile )
        qcounts = bx.intervals.operations.rightwords.HarvestRightWords ( qset, wset, isbackground=False, qcol=options.querycol, 
                        mincols=wargs[1], asfraction=wargs[0], matchstrand=wargs[2], skipchrnames=True, skipstrandnames=True, aspairs=options.pairs)
        qcounts.start()
        qcounters.append( qcounts )
        if options.background != None: 
            wset = makeFormatSafeIntervalReader( wfile )
            bset = makeFormatSafeIntervalReader( options.background )
            bcounts = bx.intervals.operations.rightwords.HarvestRightWords ( bset, wset, isbackground=True, qcol=options.querycol, 
                            mincols=wargs[1], asfraction=wargs[0], matchstrand=wargs[2], skipchrnames=True, skipstrandnames=True, aspairs=options.pairs)
            bcounts.start()
            bcounters.append(bcounts)
    
    qsize = countvalidlines(queryFile)
    if (options.background != None):
        bsize = countvalidlines(options.background)
    

    if ( options.background == None ):
        for idx in range(0,len(qcounters)):
            qcounters[idx].join()
            wfilename = wordFileArgs.keys()[idx]
            for items in qcounters[idx].getWords():
                print wfilename + "\t" + "\t".join(items)
    else:
        sigcounters = []
        for idx in range(0,len(bcounters)):
            qcounters[idx].join()
            bcounters[idx].join()
            sigcounters.append(threading.Thread(target=printsignificancescores, args=tuple([qcounters[idx], bcounters[idx], qsize, bsize, wordFileArgs.keys()[idx]])))
            sigcounters[idx].start()

def printsignificancescores( qharvester, bharvester, qsize, bsize, wfilename ):
    qw = list(qharvester.getWords())
    bw = list(bharvester.getWords())
    scores = bx.intervals.operations.rightwords.getsignificance( qw, bw, qsize, bsize )
    for line in scores:
        print wfilename + "\t" + "\t".join(line)


def countvalidlines( filename ):
    count = 0
    for line in open(filename):
        if line[0] != '#':
            count += 1
    
    return count

    

class FileFormatError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def makeFormatSafeIntervalReader( filename ):
    #checks if the file (based on name) is bed, gff or gtf
    # and returns an appropriate GenomicIntervalReader
    if ( string.upper(filename).endswith( tuple(['GTF','GFF']))):
        # typical GTF or GFF see: http://genome.ucsc.edu/FAQ/FAQformat#format1
        return bx.intervals.io.GenomicIntervalReader( open(filename), chrom_col=0, start_col=3, end_col=4, strand_col=6 )
    if ( string.upper(filename).endswith( tuple(['bed','BED']))):
        # typical BED file see: http://genome.ucsc.edu/FAQ/FAQformat#format3
        return bx.intervals.io.GenomicIntervalReader( open(filename) )
    if ( string.upper(filename).endswith( tuple(['BEDTAG', 'BTG']) ) ):
        # first six columns are same as bed -- then groups of tags, values, or tag value pairs
        # tag and value in pair separated by spaces
        # pairs or singletons in group separated by commas or semicolons
        # groups separated by tabs
        return bx.intervals.io.GenomicIntervalReader( open(filename) )
    raise FileFormatError('Files must be BEDTAG, BED, GFF, or GTF and end with the appropriate three-letter extension!')

if __name__ == "__main__":
    main()

