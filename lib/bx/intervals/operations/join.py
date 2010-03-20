"""
Join two sets of intervals using their overlap as the key.  The
intervals MUST be sorted by chrom(lexicographically),
start(arithmetically) and end(arithmetically).  This works by simply
walking through the inputs in O(n) time.
"""

import psyco_full

import math
import traceback
import fileinput
from warnings import warn


from types import ListType

from bx.intervals.io import *
from bx.intervals.operations import *
from quicksect import IntervalTree


def join(leftSet, rightSet, mincols=1, leftfill=True, rightfill=True, asfraction=False, matchStrand=STRAND_NEUTRAL, outColumns=[-1,-1]):
    # Read rightSet into memory:
    rightlen = 0
    leftlen = 0
    rightStrandCol = -1
    minoverlap = mincols
    rightTree = IntervalTree()
    
    for item in rightSet:
        if type( item ) is GenomicInterval:
            rightTree.insert( item, rightSet.linenum, item.fields )
            if rightlen == 0: rightlen = item.nfields
            if rightStrandCol == -1: rightStrandCol = item.strand_col

    for interval in leftSet:
        if leftlen == 0 and type( interval ) is GenomicInterval:
            leftlen = interval.nfields
        if not (type( interval ) is GenomicInterval):
            yield interval
        else:
            result = []
            rightTree.intersect( interval, lambda node: result.append( node ) )
            overlap_not_met = 0
            leftbases = interval.end - interval.start
            for item in result:
                rightbases = item.end - item.start
                if (asfraction==True):
                    if rightbases < leftbases:
                        mincols = rightbases
                    else:
                        mincols = leftbases
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
                    strandMatched = STRAND_INTEGER_VALUES[interval.strand] * STRAND_INTEGER_VALUES[item.other[rightStrandCol]]
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
                yield(getSelectedColumns( interval.fields, item.other, outColumns ))
            if (len(result) == 0 or overlap_not_met == len(result)) and rightfill:
                yield(getSelectedColumns( interval.fields, rightlen, outColumns ))
    if leftfill:
        def report_unvisited( node, results ):
            if not hasattr(node, "visited"):
                results.append( node )
        results = []
        rightTree.traverse( lambda x: report_unvisited( x, results ) )
        for item in results:
            yield(getSelectedColumns( leftlen, item.other, outColumns))


def getSelectedColumns( left, right, fields=[-1,-1] ):
    #left is a list or the length of the missing list
    #right is a list or the length of the missing list
    #fields tell me which columns to show you -- -1,-1 means all from both
    outfields = list()    
    try:
        leftLen=len(left)
        if (fields[1]==-1):
            map(outfields.append,list(left))
        else:
            outfields.append(left[fields[0]])
    except TypeError:
        leftLen=left
        if (fields[0]==-1):
            for x in range(leftLen): outfields.append(".")
        else:
            outfields.append(".")
    try:
        rightLen=len(right)
        if (fields[1]==-1):
            map(outfields.append,list(right))
        else:
            outfields.append(right[fields[1]])
    except TypeError:
        rightLen=right
        if (fields[1]==-1):
            for x in range(rightLen): outfields.append(".")
        else:
            outfields.append(".")
    
    return(outfields)

def interval_cmp(a, b):
    interval1 = a[0]
    interval2 = b[0]
    if not (type( interval1 ) is GenomicInterval and type( interval2 ) is GenomicInterval):
        return 0
    # Both are intervals
    if interval1.chrom == interval2.chrom:
        center1 = interval1.start + ((interval1.end - interval1.start) / 2)
        center2 = interval2.start + ((interval2.end - interval2.start) / 2)
        return center1 - center2
    else:
        if interval1.chrom > interval2.chrom:
            return 1
        else:
            return -1

    return 0

def findintersect(interval, sortedlist, mincols):
    # find range of intervals that intersect via a binary search
    # find lower bound
    x = len(sortedlist) / 2
    n = int(math.pow(2,math.ceil(math.log(len(sortedlist),2))))

    not_found = True
    not_done = True
    while not_found and not_done:
        n = n / 2
        if n == 0:
            n = 1
            not_done = False
        if x >= len(sortedlist):
            x -= n
        elif x < 0:
            x += n
        else:
            if findoverlap(sortedlist[x][0], interval) >= mincols:
                not_found = False
            else:
                comp = interval_cmp(sortedlist[x], [interval, 0])
                if comp > 0:
                    x -= n
                else:
                    x += n

    print "\t".join(sortedlist[x][0].fields)
    print "not_found = " + str(not_found)
    if not_found:
        return 0,-1

    lowerbound = x
    middlebound = x
    upperbound = x
    while (lowerbound > -1) and (findoverlap(sortedlist[lowerbound-1][0],interval) >= mincols):
        lowerbound -= 1
    while (upperbound+1 < len(sortedlist)) and (findoverlap(sortedlist[upperbound+1][0],interval) >= mincols):
        upperbound += 1

    return lowerbound, upperbound

def findoverlap(a, b):
    # overlapping
    if a.chrom == b.chrom:
        return min(a.end, b.end) - max(a.start, b.start)
    else:
        return 0
