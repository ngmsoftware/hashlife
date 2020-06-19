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






class Node:
    """  Class that holds a node """
    
    def __init__(self, sw = None, se = None, nw = None, ne = None):
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
        return id(self.sw) == id(node.sw)  and  id(self.se) == id(node.se) and id(self.nw) == id(node.nw) and id(self.ne) == id(node.ne)

    def __ne__(self, node):
        return not self.__eq__(node)

    def __hash__(self):
        return hash((id(self.sw), id(self.se), id(self.nw), id(self.ne)))
        
    
        
    
    
    
class HashLife:
    def __init__(self):

        self.memoNodes = {}
        self.memoRes = {}


        self.generation = 0
        self.worldDepth = 0
        self.fastForward = False





        # Canonical Nodes

        # CN1:       CN2:        CN3:       CN4: 
        #      . .        . x        x .        x x
        #      . .        . .        . .        . .
        # CN5:       CN6:        CN7:       CN8: 
        #      . .        . x        x .        x x
        #      . x        . x        . x        . x
        # CN9:       CN10:       CN11:      CN12: 
        #      . .        . x        x .        x x
        #      x .        x .        x .        x .
        # CN13:      CN14:       CN15:      CN16: 
        #      . .        . x        x .        x x
        #      x x        x x        x x        x x
        
        
        self.CN1 = Node(0,0,0,0)
        self.CN2 = Node(0,0,0,1)
        self.CN3 = Node(0,0,1,0)
        self.CN4 = Node(0,0,1,1)
        self.CN5 = Node(0,1,0,0)
        self.CN6 = Node(0,1,0,1)
        self.CN7 = Node(0,1,1,0)
        self.CN8 = Node(0,1,1,1)
        self.CN9 = Node(1,0,0,0)
        self.CN10 = Node(1,0,0,1)
        self.CN11 = Node(1,0,1,0)
        self.CN12 = Node(1,0,1,1)
        self.CN13 = Node(1,1,0,0)
        self.CN14 = Node(1,1,0,1)
        self.CN15 = Node(1,1,1,0)
        self.CN16 = Node(1,1,1,1)
        
   
        # List of canonical nodes to be used for canonization of pixels
        
        self.CNList = [self.CN1, self.CN2, self.CN3, self.CN4, self.CN5, self.CN6, self.CN7, self.CN8, self.CN9, self.CN10, self.CN11, self.CN12, self.CN13, self.CN14, self.CN15, self.CN16]
        
        
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
        
        if type(sw) == type(0):
            node = self.pickCanonical(sw, se, nw, ne)
        else:
            node = Node(sw, se, nw, ne)
        
        if node in self.memoNodes:
            return self.memoNodes[node]
        
        self.memoNodes[node] = node
        
        return node






    def pickCanonical(self, sw, se, nw, ne):
        """
        This function gets the for corners bits and returns a pointer to one of
        the canonical nodes (this saves a lot of memory too because we never
        repeat a depth 1 node)
        """
        sw = 1 if sw == 1 else 0
        se = 1 if se == 1 else 0
        nw = 1 if nw == 1 else 0
        ne = 1 if ne == 1 else 0
        
        for CNCandidate in self.CNList:
            if sw == CNCandidate.sw and se == CNCandidate.se and nw == CNCandidate.nw and ne == CNCandidate.ne:
                CN = CNCandidate
        
        return CN








    def generateCanonical0(self, depth):
        """
        Generate 0 filled node.
        It starts from a canonical 0 and builds the way up to depth
        """
        
        if depth == 0:
            return 0
        
        n1 = self.CN1
        for i in range(2,depth+1):
            n1 = self.createNode(n1, n1, n1, n1)
        
        return n1











    def addBorder(self, node):
        
        depth = node.depth
        
        nodeBorder = self.generateCanonical0(depth-1)
        
        
        resSW = self.createNode(nodeBorder, nodeBorder, nodeBorder, node.sw)
        resSE = self.createNode(nodeBorder, nodeBorder, node.se, nodeBorder)
        resNW = self.createNode(nodeBorder, node.nw, nodeBorder, nodeBorder)
        resNE = self.createNode(node.ne, nodeBorder, nodeBorder, nodeBorder)
    
        
    
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
                    
                    subNodes.append( self.pickCanonical(sw, se, nw, ne) )
    
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
                            node = self.pickCanonical(v[3], v[4], v[1], v[2])
                                
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
        res = self.addBorder( self.addBorder( NODES[idx-1] ) )
    
    
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


    



    def stepNode(self, node):
        """ 
        Main HashLife function!
        Takes a node and returns the result as the center (depth-1) 
        The return node is depth - 1.
        The function ignores the depth-2 border and processes only the center
        node.
        
        According to the hashLife algorithm, if node is the minimun depth (2
        in our case), we process the cneter an return a canonical node with
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
        
        

        # increment generation if we are processing the top most 
        # macrocell (whole universe)
        if node.depth == self.worldDepth:
            self.generation += 2**(self.worldDepth-2)

        # if area is 0, there is nothing to process...    
        if node.area == 0:
            return self.getCenterNode(node)
    
    
        # if we did that already... get from the memo!
        if node in self.memoRes:
            return self.memoRes[node]
    
    
        # depth 2 means we do things in the usual way
        # make the small matrix, process the pixels and return the result
        if node.depth == 2:
    
            self.nodeToAuxiliaryMatrix(node)
            
            self.processAuxiliaryMatrix()
            
            result = self.getCenterNode( self.auxiliaryMatrixToNode() )
        
        else:
            
            #   +------+--+------+--+------+
            #   |      |  |      |  |      |            
            #   |  n31 |  | n32  |  |  n33 |
            #   |      |  |      |  |      |            
            #   +------+--+------+--+------+
            #   |      |  |      |  |      | 
            #   +------+--+------+--+------+
            #   |      |  |      |  |      |            
            #   |  n21 |  | n22  |  | n23  |
            #   |      |  |      |  |      |            
            #   +------+--+------+--+------+
            #   |      |  |      |  |      | 
            #   +------+--+------+--+------+
            #   |      |  |      |  |      |            
            #   |  n11 |  | n12  |  | n13  |
            #   |      |  |      |  |      |            
            #   +------+--+------+--+------+
            #
            #   9 Auxiliary nodes that are needed to be processed to have enough
            # information to assemble the depth-1 center result
            
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
            
            res11 = self.stepNode(node11)
            res12 = self.stepNode(node12)
            res13 = self.stepNode(node13)
            res21 = self.stepNode(node21)
            res22 = self.stepNode(node22)
            res23 = self.stepNode(node23)
            res31 = self.stepNode(node31)
            res32 = self.stepNode(node32)
            res33 = self.stepNode(node33)


            # To combine the 9 results, we could subdivide one more time to get
            # sub-sub nodes and use 16 of then to assemble the denter result.
            #
            # Instead we assemble 4 subnodes whose centers will form the 
            # sw, se, nw, ne of the center cesult.
            #
            # But wait... those four new subnodes (whose centers would form)
            # our result are depth-1 as the 9 auxiliary subnodes, so why not
            # call the setpNode again? They are ready and in the right format
            # since we are going to use theiy centers to assemble the result.
            # so lets use the center in the FUTURE to assemble the result.
            #
            # It turns out that this second call to setpNode makes the final 
            # result step 2**(depth-2) ahead in time!
            #
            # For instance: If we call setpNode in an 8x8 node (depth 3)
            # we will generate 9 nodes 4x4 and call setpNode so we have 9 nodes
            # one generation ahead.
            # then, from this already one generation ahead nodes we assemble 4
            # 4x4 nodes (which would be one generation ahead) and call stepNode
            # again, so we have another gneenration and the 4 results that will
            # be used to assemble the final result will be 2 generations ahead
            #
            # If we have a 16x16 then we have:
            #
            # 9 8x8 subnodes resulting in 2 generations ahead (because they are
            # 8x8 as before)
            #            
            # From those 9 we assemble 4 more 8x8 subnodes which are already 2
            # in the future and call again which will be 2 more (now 4) 
            # generations in the future
            #
            # If we have a 32x32 then we have:
            #
            # 9 16x16 subnodes resulting in 4 generations ahead (because they are
            # 16x16 as before)
            #            
            # From those 9 we assemble 4 more 16x16 subnodes which are already 
            # 4 in the future and call again which will be 4 more (now 8) 
            # generations in the future
            #
            #  and so on... So, each level up we DOUBLE the generations
            #



            # If fastForward, we do the trick of calling setpNode again
            # if not, we just assemble the center node using the 4 centers

            if self.fastForward:


                sw = self.stepNode( self.createNode( res11, res12, res21, res22 ) )
                se = self.stepNode( self.createNode( res12, res13, res22, res23 ) )
                nw = self.stepNode( self.createNode( res21, res22, res31, res32 ) )
                ne = self.stepNode( self.createNode( res22, res23, res32, res33 ) )
                

            else:
                
                sw = self.getCenterNode( self.createNode( res11, res12, res21, res22 ) )
                se = self.getCenterNode( self.createNode( res12, res13, res22, res23 ) )
                nw = self.getCenterNode( self.createNode( res21, res22, res31, res32 ) )
                ne = self.getCenterNode( self.createNode( res22, res23, res32, res33 ) )


            result = self.createNode(sw, se, nw, ne)


        #     memoize the result
        self.memoRes[node] = result
        
        return result
    
    


    def getNodeList(self, node, nodeList = None):
        """
        Takes a node and assemble a list (mora as a hash to avoid duplicates)
        of all subnodes in all levels
        """        
        
        if nodeList == None:
            nodeList = {}
        
        if node.area > 0:
        
            if node.depth == 1:
                nodeList[node] = node
            else:
                if node.sw.area > 0:
                    nodeList[node.sw] = node.sw
                if node.se.area > 0:
                    nodeList[node.se] = node.se
                if node.nw.area > 0:
                    nodeList[node.nw] = node.nw
                if node.ne.area > 0:
                    nodeList[node.ne] = node.ne
                
                self.getNodeList(node.sw, nodeList)
                self.getNodeList(node.se, nodeList)
                self.getNodeList(node.nw, nodeList)
                self.getNodeList(node.ne, nodeList)
        
        return nodeList






    def cleanBorders(self, node):
        """
        Clean the empty nodes in the border of a root node. This is more than
        than the getCenterNode because it keeps taking the depth-a border
        until this part have any live cell
        """
        
        resNode = deepcopy(node)
        
        keepGoing = True
        while keepGoing:
            centerNode = self.getCenterNode(resNode)
            if resNode.area == centerNode.area:
                resNode = centerNode
            else:
                keepGoing = False

        return resNode













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
#                    U[corner[1]][corner[0]] = node.area
    
                    
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




    def drawNodeList(self, memo):
        """
        Draw all the subnodes in all levels of a root node.
        Assemble a subplot matrix and draw all the nodes into the subaxis
        """
        
        nNodes = len(memo)
        rows = int(sqrt(nNodes))
        cols = 1+int(nNodes/rows)
        
        widthAxis = 1.0/cols
        heightAxis = 1.0/rows
        
        print(cols,rows)

        nodeList = [x for x in memo]
        
        nodeList = sorted(nodeList, key=lambda x: x.depth)

        count = 0
        for i in range(rows):
            for j in range(cols):

                if count < nNodes:
                    node = nodeList[count]
            
                    U = self.node2image(node, 128)
                    
                    axX = 0.05*widthAxis + j/cols
                    axY = 1.0 - i/rows - heightAxis - 0.05*heightAxis
                    
                    P = plt.axes([axX, axY, widthAxis*0.9, heightAxis*0.9])
                    plt.imshow(U, aspect='auto')
                    P.xaxis.set_visible(False)
                    P.yaxis.set_visible(False)
                    

                count += 1




    def printNode(self, node):
        """
        Prints an ascii representation of a node
        """
    
        U = self.node2image(node,  128)
    
        strRes = ''
        U.reverse()
        for u in U:
            for v in u:
                strRes += '@ ' if v == 1 else '. '
            strRes += '\n'
        
    
        print(strRes)




    def drawNode(self, node, maxSize = 1024):
        """
        Draw a node 
        """
    
        U = self.node2image(node, maxSize)
            
    
        U.reverse()
    
        if self.P == None:
            self.P = plt.imshow(U, aspect='auto')
            plt.set_cmap('Greys')
            # plt.plot([len(U)/2, len(U)/2],[5, len(U)-1], 'k')
            # plt.plot([0, len(U)-1], [len(U)/2, len(U)/2], 'k')
        else:
            self.P.set_data(U)



        #plt.show()
        plt.pause(0.001)







    def plotNode(self, node, minDepth = 1):
        """
        Draws a node as a scatterplot
        """
    
        
        def recurse(node, pos, xdata, ydata):
            
            depth = node.depth
            
            if depth == minDepth:
                if depth == 1:
                    sw = node.sw
                    se = node.se
                    nw = node.nw
                    ne = node.ne
                else:
                    sw = 1 if node.sw.area > 0 else 0
                    se = 1 if node.se.area > 0 else 0
                    nw = 1 if node.nw.area > 0 else 0
                    ne = 1 if node.ne.area > 0 else 0
                    
                xdata.append(pos[0]) if sw == 1 else None
                xdata.append(pos[0]) if se == 1 else None
                xdata.append(pos[0]+1) if nw == 1 else None
                xdata.append(pos[0]+1) if ne == 1 else None
                ydata.append(pos[1]) if sw == 1 else None
                ydata.append(pos[1]+1) if se == 1 else None
                ydata.append(pos[1]) if nw == 1 else None
                ydata.append(pos[1]+1) if ne == 1 else None
            else:
                recurse(node.sw, (pos[0], pos[1]), xdata, ydata)
                recurse(node.se, (pos[0]+2**depth, pos[1]), xdata, ydata)
                recurse(node.nw, (pos[0], pos[1]++2**depth), xdata, ydata)
                recurse(node.ne, (pos[0]+2**depth, pos[1]+2**depth), xdata, ydata)
    
    
        xdata = []
        ydata = []
    
    
        recurse(node, (0,0), xdata, ydata)
            
        if self.P == None:
            self.P = plt.plot(xdata,ydata, 's')[0]
        else:
            self.P.set_xdata(xdata)
            self.P.set_ydata(ydata)



        #plt.show()
        plt.pause(0.001)







F = plt.figure(figsize=(10.0,10.0))
plt.axes([0, 0, 1, 1])

hashLife = HashLife()
hashPlot = HashPlot()






a = time.time()
#rootNode = hashLife.openMcFile('c5-adjustable-rake.mc')
rootNode = hashLife.openMcFile('metapixel-galaxy.mc')
#rootNode = hashLife.openMcFile('oscilator1.mc')
#rootNode = hashLife.openMcFile('koks galaxy.mc')
#rootNode = hashLife.openMcFile('simple.mc')
#rootNode = hashLife.openMcFile('sample.mc')
#rootNode = hashLife.openMcFile('extensible-low-period.mc')
#rootNode = hashLife.openMcFile('line-puffer-superstable.mc')
#rootNode = hashLife.addBorder( hashLife.addBorder( hashLife.addBorder( hashLife.openMcFile('ark2.mc') ) ) )
#rootNode = hashLife.openMcFile('p103079214841.mc') 
#rootNode = hashLife.openMcFile('p690-PT-Cordership-gun.mc') 
#rootNode = hashLife.openMcFile('pseudo-p34-gun.mc') 
#rootNode = hashLife.openMcFile('4 osc.mc')
#rootNode = hashLife.openMcFile('acorn.mc')
print('%.4fs'%(time.time()-a))

hashLife.worldDepth = rootNode.depth


#hashPlot.drawNodeList( hashLife.getNodeList( hashLife.cleanBorders(rootNode )) )



print('%d bytes'%(get_size(rootNode)))

# a = time.time()
# print(countPixels(rootNode))   
# print('%.4fs'%(time.time()-a))

hashLife.fastForward = True

obsceneSpeed = True

counter = 0

for i in range(200):
    
 
    
    a = time.time()
    rootNode = hashLife.addBorder( hashLife.stepNode(rootNode) )
    print('step time: %.4f s'%(time.time()-a))





    #
    # HyperSpeed
    #
    if obsceneSpeed and hashLife.fastForward:
        rootNode = hashLife.addBorder( rootNode )
        hashLife.worldDepth += 1

 
    # nodeToDraw = deepcopy(rootNode)
    # for i in range(counter):
    #     nodeToDraw = hashLife.getCenterNode(nodeToDraw)
    # hashPlot.drawNode( hashLife.getCenterNode(nodeToDraw), False, 3)



    
    a = time.time()

    # nodeToDraw = deepcopy(rootNode)
    # for j in range(counter):
    #     nodeToDraw = hashLife.getCenterNode(nodeToDraw)
    # hashPlot.P = None
    # hashPlot.drawNode(nodeToDraw, False, [5, 5, 4, 3, 1, 1, 0, 0, 0, 0][counter])
    # plt.savefig('output/%d.png'%(i))


    hashPlot.drawNode(hashLife.cleanBorders( rootNode ), 1024)
    
    # for puffer
    #hashPlot.drawNode( hashLife.getCenterNode(rootNode), False, 10)
    
    # for the metapixel
    #hashPlot.drawNode(hashLife.getCenterNode(rootNode), False, 7)
    
    print('draw time: %.4f s'%(time.time()-a))

    print(len(hashLife.memoRes), len(hashLife.memoNodes))

