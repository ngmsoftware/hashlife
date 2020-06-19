#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 21:16:11 2020

@author: ninguem
"""
from math import log2, sqrt
import matplotlib.pyplot as plt
from copy import deepcopy
import re
import time
from utils import *
from classes import *


CEL_SIZE = 12
MEMO_SCALE = 0.3

COUNT = [0]

memoRes = {}


def isEqual(node1, node2):
    if node1.depth != node2.depth:
        return False
    
    if node1.depth == 1:
        return node1.sw == node2.sw and node1.se == node2.se and node1.nw == node2.nw and node1.ne == node2.ne
    
    return isEqual(node1.sw, node2.sw) and isEqual(node1.se, node2.se) and isEqual(node1.nw, node2.nw) and isEqual(node1.ne, node2.ne)


def stringify(node):
    
    if node.depth == 1:
        return '%d%d%d%d'%(node.sw, node.se, node.nw, node.ne)
    
    return stringify(node.sw) + stringify(node.se) + stringify(node.nw) + stringify(node.ne)


class Node:
    """  Class that holds a node """
    
    def __init__(self, sw = None, se = None, nw = None, ne = None, x = 0, y = 0):
        self.x = x
        self.y = y
        
        self.sw = sw
        self.se = se
        self.nw = nw
        self.ne = ne
        if type(sw) == type(self):
            self.depth = sw.depth+1
            self.area = sw.area + se.area + nw.area + ne.area
        else:
            self.depth = 1
            self.area = sw + se + nw + ne
            
            
            
            
            
    def __eq__(self, node):
        return isEqual(self, node)

    def __ne__(self, node):
        return not isEqual(self, node)

    def __hash__(self):
        return hash(stringify(self))
        
    
        
    
    
    
class HashLife:
    def __init__(self):


        
        # auxiliary matrices to run canonimal rules
        #
        # These are pre-alocated to try to speed up the computation

        self.CM = [[0 for i in range(6)] for j in range(6)]
        self.CMr = [[0 for i in range(6)] for j in range(6)]
        
        










    def createNode(self, sw, se, nw, ne):
        """
        Creates a node. 
        Although we could create a node with Node(...) we rather use this
        function because uppon creation of a node, we already cache it in the 
        memoNodes.
        Also, if we try to create a level 1 node, we instead return the
        canonical one (which is already just a pointer)
        Even more, if we try to create a node that is already in the the memo
        we just return the pointer.
        This function alone is responsible for an AMAZING memory saving and 
        works so well that we never have a repeated node.
        """
        
        node = Node(sw, se, nw, ne)
    
        if node.depth>1:
            node.x = sw.x
            node.y = sw.y
        
        return node













    def generateCanonical0(self, depth):
        """
        Generate 0 filled node.
        It starts from a canonical 0 and builds the way up to depth
        """
        
        if depth == 0:
            return 0
        
        n1 = Node(0,0,0,0)
        for i in range(2,depth+1):
            n1 = self.createNode(deepcopy(n1), deepcopy(n1), deepcopy(n1), deepcopy(n1))
        
        return n1





    def updateNodesCoordinates(self, node, x0, y0):
        
        node.x = x0
        node.y = y0

        subNodeLen = 2**(node.depth-1)

        if node.depth >= 2: 
            self.updateNodesCoordinates(node.sw, x0, y0)
            self.updateNodesCoordinates(node.se, x0+subNodeLen, y0)
            self.updateNodesCoordinates(node.nw, x0, y0+subNodeLen)
            self.updateNodesCoordinates(node.ne, x0+subNodeLen, y0+subNodeLen)





    def addBorder(self, node):
        
        depth = node.depth
        
        nodeBorder = self.generateCanonical0(depth-1)
        
        
        resSW = self.createNode(deepcopy(nodeBorder), deepcopy(nodeBorder), deepcopy(nodeBorder), deepcopy(node.sw))
        resSE = self.createNode(deepcopy(nodeBorder), deepcopy(nodeBorder), deepcopy(node.se), deepcopy(nodeBorder))
        resNW = self.createNode(deepcopy(nodeBorder), deepcopy(node.nw), deepcopy(nodeBorder), deepcopy(nodeBorder))
        resNE = self.createNode(deepcopy(node.ne), deepcopy(nodeBorder), deepcopy(nodeBorder), deepcopy(nodeBorder))
    
        
    
        return self.createNode(resSW, resSE, resNW, resNE)
        






    def openMcFile(self, fileName):
        """
        Reads a macrocell from a Golly "Macrocell file format" 
        """
        
        
        
        
        def lineToMatrix(line):
            #
            # Converts a line of chars into a GoL matrix according to the 
            # specification of the compressed line format
            #
            
            M = [[0 for x in range(8)] for y in range(8)]
            pixelLines = line.split('$')
            for pixelLine,j in zip(pixelLines, range(8)):
                i = 0
                for pixel in pixelLine:
                    M[j][i] = 1 if pixel == '*' else 0
                    i += 1
            
            return M
        
        
        
        def matrixToNode(M):
            #
            # Returns a node from a GoL pixel matrix
            #          
            
            subNodes = []
            
            for i in range(4):
                for j in range(4):
                    sw = M[len(M)-1-2*i-0][2*j+0]
                    se = M[len(M)-1-2*i-0][2*j+1]
                    nw = M[len(M)-1-2*i-1][2*j+0]
                    ne = M[len(M)-1-2*i-1][2*j+1]
                    
                    subNodes.append( Node(sw, se, nw, ne) )
    
            SW = self.createNode(subNodes[0], subNodes[1], subNodes[4], subNodes[5])
            SE = self.createNode(subNodes[2], subNodes[3], subNodes[6], subNodes[7])
            NW = self.createNode(subNodes[8], subNodes[9], subNodes[12],subNodes[13])
            NE = self.createNode(subNodes[10], subNodes[11], subNodes[14], subNodes[15])
    
            node = self.createNode(SW, SE, NW, NE)
   
            return node
    
    
    
    
    
    
    
        fid = open(fileName, "rt")
    
    
        NODES = {}
        idx = 1
        
        NODES[0] = self.generateCanonical0(3)
    
    
        keepGoing = True
        while keepGoing:
    
            #
            # While there are lines to read...
            #
            
            line = fid.readline()
            
            if len(line)>0:
            
                # if its not a comment or a rule specifier...
                if line[0] != '#' and line[0] != '[':
                    
                    if line[0] == '.' or line[0] == '*' or line[0] == '$':
                        # if its a raw compressed matrix line
                        # create a node directly
                        
                        M = lineToMatrix(line)
                        
                        node = matrixToNode(M)


                        
                    else:
                        # else, read the indices and create a super node using
                        # the previous ones created                                                  
                        
                        v = [int(x) for x in line.split(' ')]
                        
                        if v[0] == 1:
                            node = Node(v[3], v[4], v[1], v[2])
                                
                        else:
                        
                            SW = NODES[v[3]] if v[3] != 0 else self.generateCanonical0(v[0]-1)
                            SE = NODES[v[4]] if v[4] != 0 else self.generateCanonical0(v[0]-1)
                            NW = NODES[v[1]] if v[1] != 0 else self.generateCanonical0(v[0]-1)
                            NE = NODES[v[2]] if v[2] != 0 else self.generateCanonical0(v[0]-1)
                            
                            
                            node = self.createNode( SW, SE, NW, NE)
                        
                        
                    # add the recently created node to the list of nodes
                    # the order in the list matters completelly because they
                    # are referenced in the file
                    
                    NODES[idx] = node
                       
                    idx += 1
                        
            else:
                keepGoing = False
            
        fid.close()
    
    
    
        # add two borders just to have enough space to step for some generations
        res = self.addBorder( NODES[idx-1] ) 
    
    
        # update the worldDepth (top level depth)
        self.worldDepth = res.depth
    
    
        return res
    





    



    def getCenterNode(self, node):
        """
        Returns the depth-1 center node of a node
        """
        return self.createNode(node.sw.ne, node.se.nw, node.nw.se, node.ne.sw)




    def processAuxiliaryMatrix(self):
        """
        Process the auxiliary matrix with regular GoL algorithm.
        Looping and processing each pixel.
        """
        
        M = self.CM
        Mr = self.CMr
        
        for i in range(1,5):
            for j in range(1,5):
                n1 = M[i-1][j-1]
                n2 = M[i][j-1]
                n3 = M[i+1][j-1]
                n4 = M[i+1][j]
                n5 = M[i+1][j+1]
                n6 = M[i][j+1]
                n7 = M[i-1][j+1]
                n8 = M[i-1][j]
            
                res = M[i][j]
        
                nAlive = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8
                if nAlive < 2 or nAlive > 3:
                    res = 0
                else:
                    if nAlive == 3:
                        res = 1
    
                
                Mr[i][j] = res
    






    def nodeToAuxiliaryMatrix(self, node):
        """
        Assemble a node with the auxiliary matrix pixels 
        """
        
        self.CM[1][1] = node.sw.sw
        self.CM[1][2] = node.sw.se
        self.CM[1][3] = node.se.sw
        self.CM[1][4] = node.se.se
        
        self.CM[2][1] = node.sw.nw
        self.CM[2][2] = node.sw.ne
        self.CM[2][3] = node.se.nw
        self.CM[2][4] = node.se.ne

        self.CM[3][1] = node.nw.sw
        self.CM[3][2] = node.nw.se
        self.CM[3][3] = node.ne.sw
        self.CM[3][4] = node.ne.se
        
        self.CM[4][1] = node.nw.nw
        self.CM[4][2] = node.nw.ne
        self.CM[4][3] = node.ne.nw
        self.CM[4][4] = node.ne.ne
        


    def auxiliaryMatrixToNode(self):
        """
        Assemble a node with the auxiliary matrix pixels 
        """
        
        M = self.CMr
        
        sw = self.createNode(M[1][1], M[1][2], M[2][1], M[2][2])
        se = self.createNode(M[1][3], M[1][4], M[2][3], M[2][4])
        nw = self.createNode(M[3][1], M[3][2], M[4][1], M[4][2])
        ne = self.createNode(M[3][3], M[3][4], M[4][3], M[4][4])

        return self.createNode(sw, se, nw, ne)


    



    def stepNode(self, node, plot=False):
        """ 
        Main HashLife function!
        Takes a node and returns the result as the center (depth-1) 
        The return node is depth - 1.
        The function ignores the depth-2 border and processes only the center
        node.
        
        According to the hashLife algorithm, if node is the minimun depth (2
        in out case), we process the cneter an return a canonical node with
        the result. This corresponds to one step in the GoL.
        If the node depth is 3 or more, we recursivelly call setpNode.
        Initially we call stepNode for the 9 depth-2 nodes to gather enought
        information to assemble teh center node. To assemble the center node
        we could use the sub nodes  from the 9 results and asemble a depth-1 
        center result. Which would be one generation ahead.
        Instead, since we computed all the 9 subnodes, we actually have enough
        information to assemble 4 depth-1 subnodes whose centers would form the 
        resulting depth-1 node. So, we call setepNode AGAIN on those 4 nodes
        and the center will be ahead in time. Hence when we assemble the
        resultinh depth-1 node form those denth-1 new center results we advance 
        2**(depth-2) in time!!!         
        """
        
        if plot:
            AX.cla()
            hashPlot.drawMemo(memoRes, [0.3, 0.7,0.3], [0.7, 0.3,0.3], [0,0,0], [0,0,1], 5, 2**(rootNode.depth)+1, 0, 2**(rootNode.depth), MEMO_SCALE, currentNode)
            hashPlot.drawNode(currentNode, [0.2,0.2,0.2], CEL_SIZE)
            hashPlot.drawRootGrid(5, [0,0,0])
            hashPlot.drawNodeGrid(node, [0,0.6,1])
            plt.pause(0.01)
            F.savefig('output/anim/%d.png'%(COUNT[0]))
            COUNT[0] += 1
    
        
        if node.area == 0:
            result = self.getCenterNode(node)
        
        else:
    
            # depth 2 means we do things in the usual way
            # make the small matrix, process the pixels and return the result
            if node.depth == 2:
        
                self.nodeToAuxiliaryMatrix(node)
                
                self.processAuxiliaryMatrix()
                
                result = self.getCenterNode( self.auxiliaryMatrixToNode() )
                
                resultSize = 2*(node.depth-1)
                result.x = node.x + resultSize/2
                result.y = node.y + resultSize/2
    
                if plot:
                    AX.cla()
                    if node.area>0:
                        alpha = 0.2
                    else:
                        alpha = 1.0
                    hashPlot.drawMemo(memoRes, [0.3, 0.7,0.3], [0.7, 0.3,0.3], [0,0,0], [0,0,1], 5, 2**(rootNode.depth)+1, 0, 2**(rootNode.depth), MEMO_SCALE, node)
                    hashPlot.drawNode(currentNode, [0.2,0.2,0.2], CEL_SIZE, alpha)
                    hashPlot.drawRootGrid(5, [0,0,0])
                    hashPlot.drawNodeGrid(node, [0,0,1])
    
                    hashPlot.drawNode( node , [0,1,0], CEL_SIZE)
        
                    tmpNode = self.addBorder(result)
                    
                    self.updateNodesCoordinates(tmpNode, node.x, node.y)
        
                    hashPlot.drawNode( tmpNode , [1,0,0], CEL_SIZE)
                    
                    plt.pause(0.01)
                    F.savefig('output/anim/%d.png'%(COUNT[0]))
                    COUNT[0] += 1
                    
                    if node.area > 0:
                        plt.pause(0.01)
                        for i in range(3):
                            F.savefig('output/anim/%d.png'%(COUNT[0]))
                            COUNT[0] += 1
                    
            else:
                
                
                node11 = self.createNode(node.sw.sw, node.sw.se, node.sw.nw, node.sw.ne)
                node21 = self.createNode(node.sw.nw, node.sw.ne, node.nw.sw, node.nw.se)
                node31 = self.createNode(node.nw.sw, node.nw.se, node.nw.nw, node.nw.ne)
            
                node12 = self.createNode(node.sw.se, node.se.sw, node.sw.ne, node.se.nw)
                node22 = self.createNode(node.sw.ne, node.se.nw, node.nw.se, node.ne.sw)
                node32 = self.createNode(node.nw.se, node.ne.sw, node.nw.ne, node.ne.nw)
                
                node13 = self.createNode(node.se.sw, node.se.se, node.se.nw, node.se.ne)
                node23 = self.createNode(node.se.nw, node.se.ne, node.ne.sw, node.ne.se)
                node33 = self.createNode(node.ne.sw, node.ne.se, node.ne.nw, node.ne.ne)
                
                
                
                # step the auxiliary nodes!
                
                res11 = self.stepNode(node11, plot)
                res12 = self.stepNode(node12, plot)
                res13 = self.stepNode(node13, plot)
                res21 = self.stepNode(node21, plot)
                res22 = self.stepNode(node22, plot)
                res23 = self.stepNode(node23, plot)
                res31 = self.stepNode(node31, plot)
                res32 = self.stepNode(node32, plot)
                res33 = self.stepNode(node33, plot)
    
                    
                sw = self.getCenterNode( self.createNode( res11, res12, res21, res22 ) )
                se = self.getCenterNode( self.createNode( res12, res13, res22, res23 ) )
                nw = self.getCenterNode( self.createNode( res21, res22, res31, res32 ) )
                ne = self.getCenterNode( self.createNode( res22, res23, res32, res33 ) )
    
    
                result = self.createNode(sw, se, nw, ne)
    
    
                if plot:
                    AX.cla()
                    if node.area>0:
                        alpha = 0.2
                    else:
                        alpha = 1.0
                    hashPlot.drawMemo(memoRes, [0.3, 0.7,0.3], [0.7, 0.3,0.3], [0,0,0], [0,0,1], 5, 2**(rootNode.depth)+1, 0, 2**(rootNode.depth), MEMO_SCALE, node)
                    hashPlot.drawNode(currentNode, [0.2,0.2,0.2], CEL_SIZE, alpha)
                    hashPlot.drawRootGrid(5, [0,0,0])
                    hashPlot.drawNodeGrid(node, [0,0,1])
    
                    hashPlot.drawNode( node , [0,1,0], CEL_SIZE)
        
                    tmpNode = self.addBorder(result)
                    
                    self.updateNodesCoordinates(tmpNode, node.x, node.y)
        
                    hashPlot.drawNode( tmpNode , [1,0,0], CEL_SIZE)
                    
                    plt.pause(0.01)
                    F.savefig('output/anim/%d.png'%(COUNT[0]))
                    COUNT[0] += 1
                    
                    if node.area > 0:
                        plt.pause(0.01)
                        for i in range(3):
                            F.savefig('output/anim/%d.png'%(COUNT[0]))
                            COUNT[0] += 1


        key = deepcopy(node)
        value = deepcopy(result)
        self.updateNodesCoordinates(key, 0, 0)
        self.updateNodesCoordinates(value, 0, 0)
        memoRes[key] = value

            
        return result
    
    









class HashPlot:
    """
    Utility class used to plot and draw nodes
    """

    def __init__(self):
        self.P = None


    def node2image(self, node, maxSize = 1024):
        """
        Takes a node and returns the node image. 
        This is not the matrix with cells. It takes the maxSize and computes
        the depth that it must go down to generate an image close to that size
        If it does not go to 0 depth (pixel level) it uses the area to turn a
        non-empty subnode into a pixel in the image
        """
        
        imgDepth = int(log2(maxSize)+1e-6)
        minDepth = node.depth - imgDepth
        if minDepth < 0:
            minDepth = 0
            imgDepth = node.depth
        
    
        def recurse(node, corner, U, minDepth):
            depth = node.depth-minDepth
    
            m1 = 1
            m2 = 0
    
            if minDepth != 0:
            
            
                if node.depth == minDepth and minDepth > 0:
                    U[corner[1]][corner[0]] = 1 if node.area > 0 else 0
    
                    
                else:
                    recurse(node.sw, (corner[0]               , corner[1]               ), U, minDepth)
                    recurse(node.se, (corner[0] + 2**(depth-1), corner[1]               ), U, minDepth)
                    recurse(node.nw, (corner[0]               , corner[1] + 2**(depth-1)), U, minDepth)
                    recurse(node.ne, (corner[0] + 2**(depth-1), corner[1] + 2**(depth-1)), U, minDepth)
                
            else:
                
                if node.depth == 1:
                    if node.sw == 1:
                        U[corner[1]][corner[0]] = m1
                    else:
                        U[corner[1]][corner[0]] = m2
                            
                        
                    if node.se == 1:
                        U[corner[1]][corner[0]+1] = m1
                    else:
                        U[corner[1]][corner[0]+1] = m2
                        
                    if node.nw == 1:
                        U[corner[1]+1][corner[0]] = m1
                    else:
                        U[corner[1]+1][corner[0]] = m2
                            
                        
                    if node.ne == 1:
                        U[corner[1]+1][corner[0]+1] = m1
                    else:
                        U[corner[1]+1][corner[0]+1] = m2
                        
                else:
                    recurse(node.sw, (corner[0]               , corner[1]               ), U, minDepth)
                    recurse(node.se, (corner[0] + 2**(depth-1), corner[1]               ), U, minDepth)
                    recurse(node.nw, (corner[0]               , corner[1] + 2**(depth-1)), U, minDepth)
                    recurse(node.ne, (corner[0] + 2**(depth-1), corner[1] + 2**(depth-1)), U, minDepth)
    
    
        s = 2**imgDepth
        U = [[0]*s for i in range(s)]
        recurse(node, (0,0), U, minDepth)
    
    
        return U



    def drawNodeGrid(self, node, color, x = 0, y = 0, lineWidth = 3, scale = 1):
        
        nodeSize = 2**node.depth
        
        AX.plot([x+(node.x-0.5)*scale, 
                 x+(node.x+nodeSize-0.5)*scale, 
                 x+(node.x+nodeSize-0.5)*scale, 
                 x+(node.x-0.5)*scale, 
                 x+(node.x-0.5)*scale ], 
                [y+(node.y-0.5)*scale, 
                 y+(node.y-0.5)*scale, 
                 y+(node.y+nodeSize-0.5)*scale, 
                 y+(node.y+nodeSize-0.5)*scale, 
                 y+(node.y-0.5)*scale ], color=color, lineWidth=lineWidth)
        AX.plot([x+(node.x+nodeSize/4-0.5)*scale, 
                 x+(node.x+3*nodeSize/4-0.5)*scale, 
                 x+(node.x+3*nodeSize/4-0.5)*scale, 
                 x+(node.x+nodeSize/4-0.5)*scale, 
                 x+(node.x+nodeSize/4-0.5)*scale],
                [y+(node.y+nodeSize/4-0.5)*scale, 
                 y+(node.y+nodeSize/4-0.5)*scale, 
                 y+(node.y+3*nodeSize/4-0.5)*scale, 
                 y+(node.y+3*nodeSize/4-0.5)*scale, 
                 y+(node.y+nodeSize/4-0.5)*scale ], color=color, lineWidth=lineWidth)






    def drawNode(self, node, color, size, alpha=1, x = 0, y = 0, scale = 1):
        """
        Draw a node 
        """


        if node.depth == 1:
            if node.sw:
                AX.plot([x+(node.x)*scale],[y+(node.y)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.se:
                AX.plot([x+(node.x+1)*scale],[y+(node.y)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.nw:
                AX.plot([x+(node.x)*scale],[y+(node.y+1)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.ne:
                AX.plot([x+(node.x+1)*scale],[y+(node.y+1)*scale],'o', color=color, markerSize=size, alpha=alpha)
            
            return
    
        if node.depth == 2:
            if node.sw.sw:
                AX.plot([x+(node.x)*scale],[y+(node.y)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.sw.se:
                AX.plot([x+(node.x+1)*scale],[y+(node.y)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.sw.nw:
                AX.plot([x+(node.x)*scale],[y+(node.y+1)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.sw.ne:
                AX.plot([x+(node.x+1)*scale],[y+(node.y+1)*scale],'o', color=color, markerSize=size, alpha=alpha)
                
            if node.se.sw:
                AX.plot([x+(node.x+2)*scale],[y+(node.y)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.se.se:
                AX.plot([x+(node.x+1+2)*scale],[y+(node.y)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.se.nw:
                AX.plot([x+(node.x+2)*scale],[y+(node.y+1)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.se.ne:
                AX.plot([x+(node.x+1+2)*scale],[y+(node.y+1)*scale],'o', color=color, markerSize=size, alpha=alpha)

            if node.nw.sw:
                AX.plot([x+(node.x)*scale],[y+(node.y+2)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.nw.se:
                AX.plot([x+(node.x+1)*scale],[y+(node.y+2)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.nw.nw:
                AX.plot([x+(node.x)*scale],[y+(node.y+1+2)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.nw.ne:
                AX.plot([x+(node.x+1)*scale],[y+(node.y+1+2)*scale],'o', color=color, markerSize=size, alpha=alpha)

            if node.ne.sw:
                AX.plot([x+(node.x+2)*scale],[y+(node.y+2)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.ne.se:
                AX.plot([x+(node.x+1+2)*scale],[y+(node.y+2)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.ne.nw:
                AX.plot([x+(node.x+2)*scale],[y+(node.y+1+2)*scale],'o', color=color, markerSize=size, alpha=alpha)
            if node.ne.ne:
                AX.plot([x+(node.x+1+2)*scale],[y+(node.y+1+2)*scale],'o', color=color, markerSize=size, alpha=alpha)
        else:
            
            self.drawNode(node.sw, color, size, alpha, x, y, scale)
            self.drawNode(node.se, color, size, alpha, x, y, scale)
            self.drawNode(node.nw, color, size, alpha, x, y, scale)
            self.drawNode(node.ne, color, size, alpha, x, y, scale)




    def drawRootGrid(self, depth, color, linewidth=1):
        L = 2**depth
        for i in range(L):
            AX.plot([i-0.5, i-0.5],[0-0.5, L-0.5],color=color, linewidth = linewidth)
            AX.plot([0-0.5, L-0.5],[i-0.5, i-0.5],color=color, linewidth = linewidth)

            if i==L/8 or i==3*L/8 or i==5*L/8 or i==7*L/8:
                AX.plot([i-0.5, i-0.5],[0-0.5, L-0.5],color=color, linewidth = 1+linewidth)
                AX.plot([0-0.5, L-0.5],[i-0.5, i-0.5],color=color, linewidth = 1+linewidth)

            if i==L/4 or i==3*L/4:
                AX.plot([i-0.5, i-0.5],[0-0.5, L-0.5],color=color, linewidth = 2+linewidth)
                AX.plot([0-0.5, L-0.5],[i-0.5, i-0.5],color=color, linewidth = 2+linewidth)

            if i==L/2 or i==0:
                AX.plot([i-0.5, i-0.5],[0-0.5, L-0.5],color=color, linewidth = 3+linewidth)
                AX.plot([0-0.5, L-0.5],[i-0.5, i-0.5],color=color, linewidth = 3+linewidth)
                
    
        AX.plot([L-0.5, L-0.5],[0-0.5, L-0.5],color=color, linewidth = 3+linewidth)
        AX.plot([0-0.5, L-0.5],[L-0.5, L-0.5],color=color, linewidth = 3+linewidth)



    def drawMemo(self, memo, colorNode1, colorNode2, colorGrid1, colorGrid2, cellSize, x, y, maxY, scale, selectedKey):
        maxX = 0
        Y = y
        
        for key, value in zip(memo.keys(), memo.values()):
            width = 2**key.depth

            colorGrid = colorGrid1
            lineWidth = 1
            if isEqual(key, selectedKey):
                colorGrid = colorGrid2
                lineWidth = 3

        
            self.drawNode(key, colorNode1, cellSize, 1, x, Y, scale)
            self.drawNodeGrid(key, colorGrid, x, Y, lineWidth, scale)
            self.drawNode(value, colorNode2, cellSize, 1, x+scale*(width/4), Y+scale*(width/4), scale)
            
            
            if scale*width>maxX:
                maxX = scale*width
                
            Y += scale*(width+0.5)
            if Y>maxY:
                Y = 0
                x += maxX+0.5
                maxX = 0


F = plt.figure(figsize=(20.0,10.0))
AX = F.add_axes([0, 0, 1, 1])
AX.axis('equal')

hashLife = HashLife()
hashPlot = HashPlot()

L = 64
AX.plot([0,L,L,0,0],[0,0,L,L,0])


rootNode = hashLife.openMcFile('untitled.mc')

rootNode = hashLife.getCenterNode(hashLife.getCenterNode(rootNode ))

hashLife.updateNodesCoordinates(rootNode, 0, 0)

currentNode = deepcopy(rootNode)

for i in range(2):
    
    hashLife.updateNodesCoordinates(currentNode, 0, 0)

    currentNode = hashLife.stepNode(currentNode, True)

    currentNode = hashLife.addBorder(currentNode)

    memoRes = {}

    for i in range(4):
        F.savefig('output/anim/%d.png'%(COUNT[0]))
        COUNT[0] += 1
