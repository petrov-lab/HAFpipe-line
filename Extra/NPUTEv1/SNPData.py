# Copyright (c) 2007, 2010, 2012, Adam Roberts, Leonard McMillan
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

import time
from numpy import *

#Infinite Number (for acc array)
INF = 2**16-1 # set for unsigned 16-bit integers used in acc array

#Allele Types
MAJ = '0'
MIN = '1'
I_MAJ = 'Z'
I_MIN = 'X'

class SNPData:
    '''
    This class provides a datastructure for SNP datasets being imputed with NPUTE. Input
    files should be either CSV or TXT files with lines as SNPs (separated by crlf) and
    columns as samples.  The alleles in the SNPs may be separated by commas or nothing.
    The SNPs must be in order by chromosome position.  The SNPs must be ternary with
    a majority allele, a minority allele (not required), and a unknowns (not required).
    Majoriy and minority alleles can be represented by any characters, but an unknown
    must be a '?' character.
    '''

    def __init__(self, inFile):
        '''
        Constructor reads in data file, generates mismatch vectors and upper
        triangular matrix extract indices, and reports stats of data.
        '''
        start = time.time()
        self.readInData(inFile)
        self.genMismatchVectors()
        self.genExtractIndices()
        t = int(time.time()-start+0.5)
        print 'Number of Samples: %d' % self.numSamps
        print 'Number of SNPs: %d' % len(self.snps)
        print 'Number of SDPs: %d' % len(self.sdps)
        print 'Time to Process: %d m %d s\n' % (t/60,t%60)

    def readInData(self, inFile):
        '''
        Reads in data file and normalizes alleles.
        '''
        print "Reading in SNP data from '%s'..." % inFile,
        snps = []
        nucs = []
        self.numSamps = -1        
        
        lines = file(inFile, 'r').readlines()
        for i in xrange(len(lines)):
            line = self.removeExtraChars(lines[i])
            if self.numSamps == -1:
                self.numSamps = len(line)
            elif self.numSamps != len(line):
                print '\nSNP %d has an inconsistent number of samples.' % (i+1)
                sys.exit(1)
            major, minor = self.getAlleles(line, i)
            if minor == '':
                snp = line.replace(major, MAJ)
            else:    
                snp = (line.replace(major, MAJ)).replace(minor, MIN)
            snps += [snp]
            nucs += [(major,minor)]
            
        self.snps = array(snps)
        self.sdps = set(snps)
        self.nucs = array(nucs)
        print "Done"

    def getAlleles(self, s, i):
        '''
        Returns the majority and minority allele for the SNP.  Ties are broken by making the first
        known allele the majority.
        EDITS: SIG 10/23/2019 --> only considers ACGT as possible major or minor alleles, accommodates ambiguous base codes (ie for heterozygous sites)
        '''
        s = s.upper().strip()
        q = s.count('?')
        major = ''
        minor = ''
        bases=("A" "C" "G" "T")
        cts={}
        for b in bases:
            cts[b]=s.count(b)
        major = sorted(cts, key=cts.get, reverse=True)[0]
        minor = sorted(cts, key=cts.get, reverse=True)[1]
        
        return major, minor
                
        
    def removeExtraChars(self,s):
        '''
        Helper function for reading in data.  Removes formatting characters.
        '''
        s = s.replace(',','')
        s = s.replace(' ','')
        s = s.replace('\n','')
        s = s.replace('\r','')
        s = s.replace('\t','')
        return s
        
    def genMismatchVectors(self):
        '''
        Generates a pair-wise mismatch vector for each SDP and stores it in a dictionary.
        Match = 0
        Mismatch = 2
        Unknown = 1
        '''
        print "Generating pair-wise mismatch vectors...",

        n = self.numSamps
        o = [1 for i in range(n)]
        self.vectors = dict()

        for sdp in self.sdps:
            dM = []   
            p = []
            c = []
            for i in xrange(n):
                if sdp[i] == '1':
                    p += [2]
                    c += [0]
                elif sdp[i] == '0':
                    p += [0]
                    c += [2]
                else:
                    c += [1]
                    p += [1]

            q = 0             
            for i in xrange(n):
                if p[i] == 0:
                    dM += p[i+1:]
                elif c[i] == 0:
                    dM += c[i+1:]
                else:
                    dM += o[i+1:]
            
            self.vectors[sdp] = array(dM,uint16)
        print "Done"
        
    def incorporateChanges(self):
        '''
        Uses the changes dictionary to replace unknowns with imputed values.
        '''
        print "Incorporating imputed values into the SNP data...",
        snps = self.snps
        for loc,val in self.changes.iteritems():
            locI, samp = loc
            snp = snps[locI]
            
            if val == MIN:
                iVal = I_MIN
            elif val == MAJ:
                iVal = I_MAJ
            else:
                print "unknown imputed allele "+val+" at "+loc
            
            snps[locI] = snp[0:samp] + iVal + snp[samp+1:]
        print "Done"         
        
    def outputData(self, outFile): # revised to remove bottleneck: Serge Batalov
        '''
        Outputs the SNP data to the specified file in csv format.  Lower-case values are imputed.
        '''
        print "Writing imputed data to '%s'..." % outFile,
               
        fi = file(outFile,'w')
        for i in xrange(len(self.snps)):
            snp = self.snps[i]
            out = ''
            for allele in snp:
                if allele == MAJ:
                    nuc = self.nucs[i,0]
                elif allele == MIN:
                    nuc = self.nucs[i,1]
                elif allele == I_MAJ:
                    nuc = self.nucs[i,0].lower()
                elif allele == I_MIN:
                    nuc = self.nucs[i,1].lower()
                else:
                    nuc = allele  ### EDIT SIG 10/23/2019 - allows ambiguous symbols to remain in SNP table
                out += ',' + nuc
            fi.write(out[1:] + '\n')
        print "Done"
        
    def genExtractIndices(self):
        '''
        Generate lookup table of indices for the rows of a matrix stored in upper-triangular form.
        '''
        print 'Generating indices for triangular matrix row extraction...',
        numSamps = self.numSamps
        
        extractIndices = zeros((numSamps,numSamps),int)
        for samp in xrange(numSamps):
            indices = []  

            j = samp-1
            for i in xrange(samp):
                indices += [j]
                j += numSamps -(i+2)

            indices += range(j,j+numSamps-samp)
            extractIndices[samp] = indices
        
        self.extractIndices = extractIndices
        print 'Done'

    
    def extractRow(self,tM,rowNum):
        '''
        Extracts a row from an upper-triangular mismatch matrix for the dataset.  Sets value
        for imputed sample to inf so that it will not be chosen to impute itself during
        window tests.
        '''
        row = take(tM,self.extractIndices[rowNum])
        row[rowNum] = INF
        return row        
