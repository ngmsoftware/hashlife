#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 20:00:06 2020

@author: ninguem
"""

import sys
from numpy import log2
import re

def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size



def createMatrix(size):
    
    return [[0 for i in range(size)] for j in range(size)]


def centerMatrix(M, m):
    """ center matri m into matrix M
        middle elements of M will be replaced
        matrices are assumed to be square
        len(m) < len(M)
        result is in M itself (M is changed)
    """
    
    size1 = len(M)
    size2 = len(m)
    
    pos = (size1-size2+1)//2
    
    for i in range(len(m)):
        for j in range(len(m)):
            M[pos+i][pos+j] = m[i][j]
        


def nextPowerOf2(size):
    return int(1+log2(size)-1e-6)



def openLifFile(self, fileName, padding=10):
    """
    """
    
    fid = open(fileName, "rt")

    coords = []
    patterns = []

    keepGoing = True
    while keepGoing:
        
        line = fid.readline()
        if re.search('#P',line) != None:
            coordsList = [int(x) for x in re.findall('-*[0-9]+',line)]
            if coords != None:
                coords.append((coordsList[0], coordsList[1]))
                patternTemp = []
                patternState = True
                while patternState:
                    x = fid.tell()
                    line = fid.readline()
                    if len(line) == 0:
                        patternState = False
                        keepGoing = False
                    else:
                        if line[0] != '#':
                            patternTemp.append(line)    
                        else:
                            fid.seek(x)
                            patternState = False
                patterns.append(patternTemp)
                
                
        if line == '':
            keepGoing = False
        
        
    
    fid.close()


    minX = coords[0][0]
    maxX = coords[0][0]
    minY = coords[0][1]
    maxY = coords[0][1]
    for coord, pattern in zip(coords, patterns):

        patternLengthX = len(pattern[0])
        for patternLine in pattern:
            if len(patternLine) > patternLengthX:
                patternLengthX = len(patternLine)
        patternLengthY = len(pattern)
        
        x1 = coord[0]
        x2 = coord[0]+patternLengthX
        
        y1 = coord[1]
        y2 = coord[1]+patternLengthY
        
        
        if x2>maxX:
            maxX = x2
        if y2>maxY:
            maxY = y2

        if x1<minX:
            minX = x1
        if y1<minY:
            minY = y1



    MSize = max( [ maxX-minX, maxY-minY  ] )+2*padding

    M = [[0 for i in range(MSize)] for j in range(MSize)]

    newCoords = []
    for coord in coords:
        newCoords.append( (padding+coord[0] - minX, padding+coord[1] - minY) )


    for newCoord, pattern in zip(newCoords, patterns):
        
        for patternLine, j in zip(pattern, range(len(pattern))):
            for pixel, i in zip(patternLine, range(len(patternLine))):
                if pixel == '*':
                    M[newCoord[1]+j][newCoord[0]+i] = 1

    return M
