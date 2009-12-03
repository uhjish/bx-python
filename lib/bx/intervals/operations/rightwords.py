"""
Join two sets of intervals using their overlap as the key.  The
intervals MUST be sorted by chrom(lexicographically),
start(arithmetically) and end(arithmetically).  This works by simply
walking through the inputs in O(n) time.
"""

import psyco_full

import math
import re
import shlex
import scipy.stats
import threading
import traceback
import fileinput
from warnings import warn
from cStringIO import StringIO

from bx.intervals.io import *
from bx.intervals.operations import *
from quicksect import IntervalTree


def getpairs(leftSet, rightSet, leftCol, mincols=1, asfraction=False, matchStrand=STRAND_NEUTRAL, skipChrNames=True, skipStrandNames=True):
    # Read leftSet into memory:
    leftlen = 0
    rightlen = 0
    leftStrandCol = -1
    minoverlap = mincols
    leftTree = IntervalTree()
    rightCols = list()
    for item in leftSet:
        if type( item ) is GenomicInterval:
            leftTree.insert( item, leftSet.linenum, item.fields )
            if leftlen == 0: leftlen = item.nfields
            if leftStrandCol == -1: leftStrandCol = item.strand_col

    for interval in rightSet:
        if rightlen == 0 and type( interval ) is GenomicInterval:
            rightlen = interval.nfields
            rightCols = range(rightlen)
            #remove the useless columns
            rightCols.remove( interval.start_col )
            rightCols.remove( interval.end_col )
            if skipChrNames:
                rightCols.remove( interval.chrom_col )
            if skipStrandNames:
                rightCols.remove( interval.strand_col )
        if not (type( interval ) is GenomicInterval):
            yield interval
        else:
            result = []
            leftTree.intersect( interval, lambda node: result.append( node ) )
            overlap_not_met = 0
            rightbases = interval.end - interval.start
            for item in result:
                leftbases = item.end - item.start
                if (asfraction==True):
                    if leftbases < rightbases:
                        mincols = leftbases
                    else:
                        mincols = rightbases
                    mincols = math.floor(mincols * minoverlap)
                if item.start in range(interval.start,interval.end+1) and item.end not in range(interval.start,interval.end+1):
                    overlap = interval.end-item.start
                elif item.end in range(interval.start,interval.end+1) and item.start not in range(interval.start,interval.end+1):
                    overlap = item.end-interval.start
                elif item.start in range(interval.start,interval.end+1) and item.end in range(interval.start,interval.end+1):
                    overlap = item.end-item.start
                else:   #the intersecting item's start and end are outside the interval range
                    overlap = interval.end-interval.start
                if overlap < mincols:
                    overlap_not_met += 1
                    continue
                else:
                    #check strand
                    strandMatched = STRAND_INTEGER_VALUES[interval.strand] * STRAND_INTEGER_VALUES[item.other[leftStrandCol]]
                    if (strandMatched == -1 and matchStrand > 0):
                        #needed match but found a complement
                        overlap_not_met += 1
                        continue
                    if (strandMatched == 1 and matchStrand < 0):
                        #needed complement but found a match
                        overlap_not_met += 1
                        continue
                    if (strandMatched == 0 and (matchStrand < -1 or matchStrand > 1)):
                        #strict criteria but only permissive match found
                        overlap_not_met += 1
                        continue
                #strand criteria met
                setattr( item, "visited", True )
                leftTerm = item.other[leftCol]
                for col in rightCols:
                    #take each field that's not a number
                    #split it on semicolons, commas, and spaces
                    #output the word and the leftTerm as being associated
                    #curcol = re.sub("\;|\,","\t",interval.fields[col])
                    curcol= interval.fields[col]
                    lexer = shlex.shlex(curcol)
                    lexer.whitespace='\t\r\n\,\;'
                    lexer.wordchars += ":'"
                    lexer.whitespace_split=True
                    lexer.quotes='"'
                        
                    for item in lexer:
                        item = item.strip()
                        if (item == "."): continue
                        try:
                            float(item) 
                        except ValueError:
                            yield [item, leftTerm]
                            
def tabsplit(strToSplit):
    return strToSplit.split("\t")

def getcounts( pairs, showMembers=False ):
    pairs = map(tabsplit,set(map("\t".join, pairs)))
    counts = {}
    words = {}
    for pair in pairs:
        try:
            counts[ pair[0] ] += 1
            if showMembers:
                words[ pair[0] ].write(pair[1])
                words[ pair[0] ].write(",")
        except KeyError:
            counts[ pair[0] ]= 1
            if showMembers:
                words[ pair[0] ] = StringIO()
                words[ pair[0] ].write(pair[1])
                words[ pair[0] ].write(",")
    for word, count in counts.iteritems():
        if showMembers:
            yield [word, str(count), words[ word ].getvalue()]
        else:
            yield [word, str(count)]

def getsignificance( qcounts, bgcounts, qsize, bgsize ):
    bgdict = {}
    for bct in bgcounts:
        bgdict[ bct[0] ] = bct[1]
        
    for qct in qcounts:
        try:
            bgct = bgdict[ qct[0] ]
            pval = scipy.stats.binom_test( qct[1], qsize, float(bgct) / float(bgsize) )
        except:
            pval= -1
        yield [ str(qct[0]), str(qct[1]), str(pval), str(bgct), str(qct[2]) ]

class Significator(threading.Thread):
    def __init__( self, qharvester, bharvester, qsize, bsize ):
        self.qharvester = qharvester
        self.bharvester = bharvester
        self.qsize = qsize
        self.bsize = bsize
        self.rows = []
        threading.Thread.__init__(self)
        self.cond=threading.Condition()
        self.done=False
        
    def run(self):
        self.cond.acquire()
        self.rows = getsignificance( self.qharvester.getWords(), self.bharvester.getWords(), self.qsize, self.bsize )
        self.done=True
        self.cond.notify()
        self.cond.release()
    
    def getWords(self):
        self.cond.acquire()
        while not self.done:
            self.cond.wait()
        self.cond.release()
        return self.rows

         
class HarvestRightWords(threading.Thread):
    def __init__( self, qset, wset, isbackground, qcol, mincols, asfraction, matchstrand, skipchrnames, skipstrandnames, aspairs ):
        self.qset = qset
        self.wset = wset
        self.isbackground = isbackground
        self.qcol = qcol
        self.mincols = mincols
        self.asfraction = asfraction
        self.matchstrand = matchstrand
        self.skipchrnames = skipchrnames
        self.skipstrandnames = skipstrandnames
        self.aspairs = aspairs
        self.words = []
        threading.Thread.__init__(self)
        self.cond=threading.Condition()
        self.done=False
        
    def run(self):
        self.cond.acquire()
        self.words = getpairs( self.qset, self.wset, self.qcol, self.mincols, self.asfraction, self.matchstrand, self.skipchrnames, self.skipstrandnames ) 
        if self.aspairs == False: self.words = getcounts( self.words, showMembers=(not self.isbackground) )
        self.done=True
        self.cond.notify()
        self.cond.release()
    
    def getWords(self):
        self.cond.acquire()
        while not self.done:
            self.cond.wait()
        self.cond.release()
        return self.words
    


