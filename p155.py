#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 17 09:14:11 2020

@author: ninguem
"""

import sympy
from functools import reduce
import time



def fracReduce(frac):
    fac1 = sympy.factorint(frac[0]) if frac[0] != 1 else {1:1}
    fac2 = sympy.factorint(frac[1]) if frac[1] != 1 else {1:1}
    
    for primeFactor in fac1.keys():
        if primeFactor in fac2:
            if fac1[primeFactor] >= fac2[primeFactor]:
                fac1[primeFactor] -= fac2[primeFactor]
                fac2[primeFactor] = 0
            else:
                fac2[primeFactor] -= fac1[primeFactor]
                fac1[primeFactor] = 0
         
           
         
    res1 = reduce( lambda x,y : x*y,  [x[0]**x[1] for x in zip(fac1.keys(), fac1.values()) ] )
    res2 = reduce( lambda x,y : x*y,  [x[0]**x[1] for x in zip(fac2.keys(), fac2.values()) ] )

    return (res1, res2)


def parallel(op1, op2, memo):
    if (op1, op2) in memo or (op2, op1) in memo:
        return memo[(op1, op2)]
    
    res = ( op1[0]*op2[0], op1[0]*op2[1] + op1[1]*op2[0] )
    res = fracReduce(res)
    
    memo[(op1, op2)] = res
    memo[(op2, op1)] = res
    
    return res

def serie(op1, op2, memo):
    if (op1, op2) in memo or (op2, op1) in memo:
        return memo[(op1, op2)]

    res = ( op1[0]*op2[1] + op1[1]*op2[0], op1[1]*op2[1] )

    res = fracReduce(res)

    memo[(op1, op2)] = res
    memo[(op2, op1)] = res

    return res


def evalCircuit(ckt, memo):
    
    if ckt in memo:
        return memo[ckt]
    
    
    res = (0,1)

    cktPointer = 0
    stack = []

    memoResS = {}
    memoResP = {}

    while cktPointer < len(ckt):

        op = ckt[cktPointer]
        
        if op == 'P':
            op1 = stack.pop()
            op2 = stack.pop()
            res = parallel(op1, op2, memoResP)
            stack.append(res)
            
        if op == 'S':
            op1 = stack.pop()
            op2 = stack.pop()
            res = serie(op1, op2, memoResS)
            stack.append(res)
            
        if op == 'r':
            stack.append((1,1))
            
        cktPointer += 1

    res = fracReduce(stack[0])

    memo[ckt] = res

    return res






    






def generateAllPOssibleCircuits(maxNRs, currentCircuit, currentPendingRs, evalStackSize, memoCkt, memoValues, memoEval, currentAmount):


    # if currentCircuit in memoCkt:
    #     return memoCkt[currentCircuit]
        

    if len(currentCircuit) == 2*maxNRs-1:
        res = evalCircuit(currentCircuit, memoEval)
        if res in memoValues:
            result = currentAmount
        else:
            memoValues[res] = 1
            result = currentAmount + 1
        
        memoCkt[currentCircuit] = 1
        
    else:
    
        if maxNRs == currentPendingRs:
            
            result = generateAllPOssibleCircuits(maxNRs, currentCircuit + 'S', currentPendingRs, evalStackSize - 1, memoCkt, memoValues, memoEval, currentAmount) + generateAllPOssibleCircuits(maxNRs, currentCircuit + 'P', currentPendingRs, evalStackSize - 1, memoCkt, memoValues, memoEval, currentAmount)
            
        else:
    
            if evalStackSize == 1:
                
                res = evalCircuit(currentCircuit, memoEval)
                if res in memoValues:
                    result = generateAllPOssibleCircuits(maxNRs, currentCircuit + 'r', currentPendingRs + 1, evalStackSize + 1, memoCkt, memoValues, memoEval, currentAmount)
                else:
                    memoValues[res] = 1
                    result = 1 + generateAllPOssibleCircuits(maxNRs, currentCircuit + 'r', currentPendingRs + 1, evalStackSize + 1, memoCkt, memoValues, memoEval, currentAmount)

                memoCkt[currentCircuit] = 1

                    
            else:
                result = generateAllPOssibleCircuits(maxNRs, currentCircuit + 'r', currentPendingRs + 1, evalStackSize + 1, memoCkt, memoValues, memoEval, currentAmount) + generateAllPOssibleCircuits(maxNRs, currentCircuit + 'S', currentPendingRs, evalStackSize - 1, memoCkt, memoValues, memoEval, currentAmount) +  generateAllPOssibleCircuits(maxNRs, currentCircuit + 'P', currentPendingRs, evalStackSize - 1, memoCkt, memoValues, memoEval, currentAmount)


    #memoCkt[currentCircuit] = result

    return result 














# memoCkt = {}
# memoValues = {}
# memoEval = {}

# for i in range(1,8):
#     N = generateAllPOssibleCircuits(i, "r", 1, 1, memoCkt, memoValues, memoEval, 0)
#     print(len(memoValues))






def processLists2(maxNR):

    memoS = {}    
    memoP = {}    

    bigList = {(1,1):1}

    smallList = {(1,1):1}
    
    keepGoing = True
    while keepGoing:
        keepGoing = False
        newList = {}
        for bigL in bigList.keys():
            for smallL in smallList.keys():
                nR = bigList[bigL] + smallList[smallL]
                if nR <= maxNR:
                    rSerie = serie(bigL, smallL, memoS)
                    rParallel = parallel(bigL, smallL, memoP)

 
                    if rSerie not in bigList:
                        keepGoing = True
                        newList[rSerie] = nR
                    else:
                        bigList[rSerie] = nR if bigList[rSerie] > nR else bigList[rSerie]

                    if rParallel not in bigList:
                        keepGoing = True
                        newList[rParallel] = nR 
                    else:
                        bigList[rParallel] = nR if bigList[rParallel] > nR else bigList[rParallel]

        
        newNewList = {}
        for bigL in newList.keys():
            for smallL in newList.keys():
                nR = newList[bigL] + newList[smallL]
                if nR <= maxNR:
                    rSerie = serie(bigL, smallL, memoS)
                    rParallel = parallel(bigL, smallL, memoP)
                    
                    if rSerie not in newList:
                        keepGoing = True
                        newNewList[rSerie] = nR
                    else:
                        newList[rSerie] = nR if newList[rSerie] > nR else newList[rSerie]

                    if rParallel not in newList:
                        keepGoing = True
                        newNewList[rParallel] = nR 
                    else:
                        newList[rParallel] = nR if newList[rParallel] > nR else newList[rParallel]

        
        
        bigList = {**smallList, **bigList}
        smallList = {**newNewList, **newList}
        
    return bigList
        
        
        
        
        
        
        
        
        
        
        
        
def processLists(maxNR):

    memoS = {}    
    memoP = {}  

    bigList = {(1,1):1}

    keepGoing = True
    while keepGoing:
        keepGoing = False
        newList = {}
        for l1 in bigList.keys():
            for l2 in bigList.keys():
    
                nR = bigList[l2] + bigList[l1]
                if nR <= maxNR:
                    rSerie = serie(l1, l2, memoS)
                    rParallel = parallel(l1,l2, memoP)

                    if rSerie not in bigList:
                        keepGoing = True
                        newList[rSerie] = nR
                    else:
                        bigList[rSerie] = nR if bigList[rSerie] > nR else bigList[rSerie]

                    if rParallel not in bigList:
                        keepGoing = True
                        newList[rParallel] = nR 
                    else:
                        bigList[rParallel] = nR if bigList[rParallel] > nR else bigList[rParallel]


        bigList = {**newList, **bigList}
        
        
    return bigList
        
        
a = time.time()
L = processLists(4)
a = time.time()-a
        
        
        
        
        
        
        
        
        
        
        
    
    



