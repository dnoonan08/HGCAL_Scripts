#!/usr/bin/env python

import numpy as np
import pandas as pd
import awkward
import time

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--threshold','--thr', default=1.5,choices=[0.5,1.0,1.5,2.0,2.5],type=float)
parser.add_argument('-N', default=50000, type=int)
parser.add_argument('--layer', default=5, type=int)
parser.add_argument('--job', default=0, type=int)
parser.add_argument('--fill', default=-1, type=int)
parser.add_argument('--fillAgain', default=-1, type=int)
parser.add_argument('--add', default=0, type=int)
args = parser.parse_args()

fillMod = args.fill>0
fillAgainMod = args.fillAgain>0

fill = args.fill
if not fillMod:
    fill = args.N + 10

fillAgain = args.fillAgain
if not fillAgainMod:
    fillAgain = args.N + 10

addTC = args.add


threshold = args.threshold
layer = args.layer
N = args.N

occDF = pd.read_csv("HGCAL_Occupancies_JB_v9imp2.csv")

occDF = occDF[(occDF.threshold==threshold) & (occDF.layer==layer)]
#occDF = occDF[occDF.motherboard==1]
occDF.set_index(occDF.wafer,inplace=True)

wafers = occDF.wafer.unique()


from numba import vectorize, int64, float64

@vectorize([int64(int64),float64(float64)])
def getBits(nTC):
    header   = 4 + 5 + 1 ### header(4) + bxID(5) + addressMode(1)
    towerSum = 8
    if nTC < 8:
        total = header + towerSum + 3 + 13*nTC
    else:
        total = header + towerSum + 48 + 7*nTC

    return total

# @vectorize([int64(int64),float64(float64)])
# def getBits(nTC):
#     header   = 20
#     towerSum = 8
#     if nTC < 8:
#         total = header + towerSum + 13*nTC
#     elif nTC < 42:
#         total = header + towerSum + 48 + 7*nTC
#     else:
#         total = header + towerSum + 7*48

#     return total



@vectorize([int64(int64,int64),int64(int64,float64)])
def moveBuffer(x,nLink):
    global Buffer
    Buffer[nLink].content[x:(x+98)] = Buffer[nLink].content[x+1:(x+99)]
    return 0

@vectorize([int64(int64,int64),int64(int64,float64)])
def moveBuffer_Latency(x,nLink):
    global Buffer_Latency
    Buffer_Latency[nLink].content[x:x+98] = Buffer_Latency[nLink].content[x+1:x+99]
    return 0

@vectorize([int64(int64,int64),int64(int64,float64)])
def clearBuffer(x,nLink):
    global Buffer
    Buffer[nLink].content[x:(x+100)] = np.array([0]*100)
    return 0

@vectorize([int64(int64,int64),int64(int64,float64)])
def clearBuffer_Latency(x,nLink):
    global Buffer_Latency
    Buffer_Latency[nLink].content[x:(x+100)] = np.array([0]*100)
    return 0


totstart = time.time()
tstart = time.time()

#links = range(1,13)
links = [x*.5 for x in range(2,25)]

Buffer = {}
Buffer_Latency = {}
OverflowCount = {}
ResetCount = {}
latencyTracker = {}

countBX = {}

overflowThreshold = 12

bitCount = []
tcCount = []

buffLen = 100

for nLink in links:
    Buffer[nLink] = awkward.fromiter( np.array([[0]*buffLen]*len(wafers),np.int16) ) 
    Buffer_Latency[nLink] = awkward.fromiter( np.array([[-1]*buffLen]*len(wafers),np.int16) ) 
    OverflowCount[nLink] = np.array([0]*len(wafers),np.int32)
    ResetCount[nLink] = np.array([0]*len(wafers),np.int32)

    Buffer[nLink].stops = Buffer[nLink].starts.copy()
    Buffer_Latency[nLink].stops = Buffer_Latency[nLink].starts.copy()

    latencyTracker[nLink] =  awkward.fromiter( np.array([[-1]*N]*len(wafers),np.int16) )
    countBX[nLink] = latencyTracker[nLink].starts.copy()

elapsed = time.time()-tstart
print(elapsed)

tstart = time.time()
for BX in range(1,N+1):
    occDF['nTC'] = np.random.poisson(occDF.mean_occupancy)
    occDF.nTC += addTC

    if BX%fill==0:
        if fillMod:
            occDF.nTC = 48

    if BX%fill==fillAgain:
        if fillAgainMod:
            occDF.nTC = 48

    occDF['nBits'] = getBits(occDF.nTC)

    tcCount.append(occDF.nTC.values.copy())
    bitCount.append(occDF.nBits.values.copy())
    for nLink in links:
        bufferNotEmpty = Buffer[nLink].starts<Buffer[nLink].stops

        Buffer[nLink].content[ Buffer[nLink].starts[bufferNotEmpty] ] -= int(nLink*32)

        #check if the leading buffer is empty
        isEmpty = Buffer[nLink].content[ Buffer[nLink].starts ]< 0

        while isEmpty.sum()>0:

            #get indices of empty buffers
            emptyIdx = Buffer[nLink].starts[isEmpty]

            #Check how long the buffer is
            lenBuffer = Buffer[nLink].stops - Buffer[nLink].starts
            #if buffer is empty, and there are at least two entries, start emptying the second entry as well
            pushForward = Buffer[nLink].starts[(isEmpty) & (lenBuffer>1)]        
            Buffer[nLink].content[ pushForward + 1 ] += Buffer[nLink].content[ pushForward ]
            moveBuffer(Buffer[nLink].starts[isEmpty],nLink)

            Buffer[nLink].stops[isEmpty] -= 1
            Buffer_Latency[nLink].stops[isEmpty] -= 1

            latencyTracker[nLink].content[countBX[nLink][isEmpty]] = Buffer_Latency[nLink].content[Buffer_Latency[nLink].starts[isEmpty]]
            countBX[nLink][isEmpty] += 1
            
            moveBuffer_Latency(Buffer_Latency[nLink].starts[isEmpty],nLink)
            isEmpty = Buffer[nLink].content[ Buffer[nLink].starts ]<0


        #add bits to buffer                                                                                                                 


        bits = occDF.nBits.values

#        bits = (Buffer_Latency[nLink].content[Buffer_Latency[nLink].starts]>=overflowThreshold)*28 + occDF.nBits.values*(Buffer_Latency[nLink].content[Buffer_Latency[nLink].starts]<overflowThreshold)

        OverflowCount[nLink] += Buffer_Latency[nLink].content[Buffer_Latency[nLink].starts]>=overflowThreshold
        
        Buffer[nLink].content[ Buffer[nLink].stops ] += bits
        Buffer[nLink].stops += 1

        Buffer_Latency[nLink].content[Buffer_Latency[nLink].stops] = 0
        Buffer_Latency[nLink].stops += 1
        Buffer_Latency[nLink].content[Buffer_Latency[nLink].content>-1] += 1

        #find buffers already over 75
        bufferFull = Buffer_Latency[nLink].content[Buffer_Latency[nLink].starts]>75
        if bufferFull.sum()>0:
            clearBuffer(Buffer[nLink].starts[bufferFull],nLink)
            clearBuffer_Latency(Buffer_Latency[nLink].starts[bufferFull],nLink)
            ResetCount[nLink] += bufferFull


elapsed = time.time()-tstart
print(elapsed)
print(elapsed/N*1.)
print(time.time()-totstart)


tcCount = np.array(tcCount)
bitCount = np.array(bitCount)

occDF.nBits = getBits(occDF.mean_occupancy)


import pickle
BufferOutput = "Output_Buffer_thr%s_layer%i_job%i.pkl"%(("%.1f"%threshold).replace('.','p'),layer,args.job)
if addTC>0:
    BufferOutput = BufferOutput.replace('_job','_add%iTC_job'%args.add)
if fillMod:
    BufferOutput = BufferOutput.replace('_job','_fillEvery%i_job'%args.fill)
    if fillAgainMod:
        BufferOutput = BufferOutput.replace('_job','_againEvery%i_job'%args.fillAgain)


with open(BufferOutput,'wb') as f:
    pickle.dump(tcCount,f)
    pickle.dump(bitCount,f)
    pickle.dump(latencyTracker,f)
    pickle.dump(OverflowCount,f)
    pickle.dump(ResetCount,f)
    pickle.dump(wafers,f)
    pickle.dump(occDF.nBits.values,f)
