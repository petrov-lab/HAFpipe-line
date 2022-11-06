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

from numpy import *

class CircularQueue:
    '''
    This class implements a data structure for multiple, nested circular queues
    using a numpy array.  The class has several extra useful functions used by NPUTE
    to implement sliding window(s).
    '''
    def __init__(self,Ls,width):
        '''
        Constructor initializes datastructure.
        '''
        height = max(Ls)*2+2
        self.height = height
        self.Ls = Ls
        self.queue = zeros((height,width),uint16)
        self.half = max(Ls)
        self.mid = -1


    def getEnds(self,i):
        '''
        Returns the bounding rows of the nested queue specified by i.
        Used to calculate mismatch vector over window i.
        '''
        top = self.queue[(self.mid+self.Ls[i])%self.height]
        bottom = self.queue[(self.mid-self.Ls[i]-1)%self.height]
        return top,bottom

    def enqueue(self,e):
        '''
        Adds element e to the end of the largest queue.
        '''
        self.incrementMid()
        nextIn = (self.mid+self.half)%self.height
        self.queue[nextIn] = e

    def incrementMid(self):
        '''
        Increments the index to the center element of all queues.
        '''
        self.mid = (self.mid+1) % self.height

    def getMid(self):
        '''
        Returns the center element of all queues. Used to remove bias
        when imputing called values during window testing.
        '''
        return self.queue[self.mid]
