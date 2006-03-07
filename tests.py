import unittest
import bx.bitset_tests
import bx.align.score_tests
import bx.intervals.io
import bx.phylo.newick_tests

suite = unittest.TestSuite( [ bx.phylo.newick_tests.suite,
                              bx.bitset_tests.suite, 
                              bx.align.score_tests.suite,
                              bx.intervals.io.suite ] )