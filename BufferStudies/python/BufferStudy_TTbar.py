#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import awkward
import time
import pickle
import os
from subprocess import Popen,PIPE

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--threshold','--thr', default=1.5,choices=[0.5,1.0,1.5,2.0,2.5],type=float)
parser.add_argument('-N', default=50000, type=int)
parser.add_argument('--layer', default=5, type=int)
parser.add_argument('--job', default=0, type=int)
parser.add_argument('--fill', default=-1, type=int)
parser.add_argument('--fillAgain', default=-1, type=int)
parser.add_argument('--fillSize', default=48, type=int)
parser.add_argument('--add', default=0, type=int)
parser.add_argument('--wordSize', '--word',default=32, type=int)
parser.add_argument('--newFormat', default = False, action='store_true')
parser.add_argument('--newFormatNoHeader', default = False, action='store_true')
parser.add_argument('--newFormatPaul', default = False, action='store_true')
parser.add_argument('--links', dest='linkList',action='append')

args = parser.parse_args()


threshold = args.threshold
layer = args.layer
N = args.N
fill = args.fill
fillAgain = args.fillAgain
fillSize = args.fillSize
addTC = args.add
job=args.job

wordSize=args.wordSize

# threshold=1.5
# layer=5
# N = 500
# fill = 0
# fillAgain = 0
# addTC = 0
# job = 0

thrStr = ("thr%.1f"%threshold).replace('.','p')

fillMod = fill>0
fillAgainMod = fillAgain>0

if not fillMod:
    fill = N + 10
if not fillAgainMod:
    fillAgain = N + 10



waferOccupancyFile = 'waferOccupancy_layer%i_%s.pkl'%(layer,thrStr)

if not os.path.exists(waferOccupancyFile):
    print("Copying file locally:")
    cmd = 'xrdcp root://cmseos.fnal.gov//store/user/dnoonan/HGCAL_Concentrator/%s .'%waferOccupancyFile
    print(cmd)
    stdout,stderr = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()

with open(waferOccupancyFile,'rb') as _file:
    waferList = pickle.load(_file,encoding='bytes')
    waferOccup = pickle.load(_file,encoding='bytes')

nWafer = len(waferList)
for i in range(nWafer):
    np.random.shuffle(waferOccup[i])
    
buffLen = 50


from numba import vectorize, int64, float64

@vectorize([int64(int64,int64),float64(float64,int64)])
def getBitsOrig(nTC,wordSize=-1):
    header   = 4 + 5 + 1 ### header(4) + bxID(5) + addressMode(1)                                                                                           
    towerSum = 8

    if nTC < 7:
        # head + tower + number of cells + nTC * (7 for data + 6 for address)
        total = header + towerSum + 3 + 13*nTC
    else:
        # head + tower + trigger bit map + nTC * (7 for data)
        total = header + towerSum + 48 + 7*nTC

    # round up to word size
    total = int(np.ceil(total/wordSize)*wordSize)
    #force to 16 bits if empty in 16 bit word size
    if nTC==0 and wordSize==16:
        total=16

    return total

#address mode (2 bits)
#  0 - empty, no TC information
#  1 - only one TC (send 6-bit address + 7-bit energy)
#  2 - less than or equal 7 TC, send by address mode (6-bit address + 7-bit energy per TC), include 3 bit for number of TC's
#  3 - full TC, send 48-bit bitmap + 7 

@vectorize([int64(int64,int64),float64(float64,int64)])
def getBitsNew(nTC,wordSize):
    header   = 4 + 5 + 2 ### header(4) + bxID(5) + addressMode(2)                                                                                           
    towerSum = 8
    if nTC == 0:
        total = header + towerSum
    elif nTC == 1:
        total = header + towerSum + 13
    elif nTC < 7:
        total = header + towerSum + 3 + 13*nTC
    else:
        total = header + towerSum + 48 + 7*nTC
    total = int(np.ceil(total/wordSize)*wordSize)

    return total

#address mode (2 bits)
#  0 - empty, no TC information
#  1 - only one TC (send 6-bit address + 7-bit energy)
#  2 - less than or equal 7 TC, send by address mode (6-bit address + 7-bit energy per TC), include 3 bit for number of TC's
#  3 - full TC, send 48-bit bitmap + 7 

@vectorize([int64(int64,int64),float64(float64,int64)])
def getBitsNoHeader(nTC,wordSize):
    header   = 5 + 2 ### header(4) + bxID(5) + addressMode(2)                                                                                           
    towerSum = 8
    if nTC == 0:
        total = header + towerSum
    elif nTC == 1:
        total = header + towerSum + 13
    elif nTC < 7:
        total = header + towerSum + 3 + 13*nTC
    else:
        total = header + towerSum + 48 + 7*nTC
    total = int(np.ceil(total/wordSize)*wordSize)

    return total

@vectorize([int64(int64,int64),float64(float64,int64)])
def getBitsPaul(nTC, wordSize):
    nWords16 = { 0  : 0,
                 1  : 2,
                 2  : 4, 3  : 4,
                 4  : 6,
                 5  : 8, 6  : 8, 
                 7  : 12, 8  : 12, 9  : 12,
                 10 : 14, 11 : 14, 12 : 14,
                 13 : 16, 14 : 16, 15 : 16,
                 16 : 18, 17 : 18, 18 : 18,
                 19 : 20, 20 : 20, 21 : 20,
                 22 : 22, 23 : 22, 24 : 22,
                 25 : 24, 26 : 24, 27 : 24,
                 28 : 26, 29 : 26, 30 : 26,
                 31 : 28, 32 : 28, 33 : 28,
                 34 : 30, 35 : 30, 36 : 30,
                 37 : 32, 38 : 32, 39 : 32,
                 40 : 34, 41 : 34, 42 : 34,
                 43 : 36, 44 : 36, 45 : 36,
                 46 : 38, 47 : 38, 48 : 38,
             }
    return nWords16[int(nTC)]*16



getBits = getBitsOrig
if args.newFormat:
    getBits = getBitsNew
if args.newFormatNoHeader:
    getBits = getBitsNoHeader
if args.newFormatPaul:
    getBits = getBitsPaul


@vectorize([int64(int64,int64),int64(int64,float64)])
def moveBuffer(x,nLink):
    global Buffer
    global buffLen
    Buffer[nLink].content[x:(x+buffLen-2)] = Buffer[nLink].content[x+1:(x+buffLen-1)]
    # Buffer[nLink].content[x:(x+98)] = Buffer[nLink].content[x+1:(x+99)]
    return 0

@vectorize([int64(int64,int64),int64(int64,float64)])
def moveBuffer_Latency(x,nLink):
    global Buffer_Latency
    global buffLen
    Buffer_Latency[nLink].content[x:x+buffLen-2] = Buffer_Latency[nLink].content[x+1:x+buffLen-1]
#    Buffer_Latency[nLink].content[x:x+98] = Buffer_Latency[nLink].content[x+1:x+99]
    return 0

@vectorize([int64(int64,int64),int64(int64,float64)])
def clearBuffer(x,nLink):
    global Buffer
    global buffLen
    Buffer[nLink].content[x:(x+buffLen)] = np.array([0]*buffLen)
    return 0

@vectorize([int64(int64,int64),int64(int64,float64)])
def clearBuffer_Latency(x,nLink):
    global Buffer_Latency
    global buffLen
    Buffer_Latency[nLink].content[x:(x+buffLen)] = np.array([0]*buffLen)
    return 0


totstart = time.time()
tstart = time.time()
if args.linkList is None:
    links = [x*.5 for x in range(2,25)]
else:
    links = [float(x) for x in args.linkList]

Buffer = {}
Buffer_Latency = {}
OverflowCount = {}
ResetCount = {}
latencyTracker = {}

countBX = {}

overflowThreshold = 12

bitCount = []
tcCount = []

for nLink in links:
    Buffer[nLink] = awkward.fromiter( np.array([[0]*buffLen]*nWafer,np.int16) )
    Buffer_Latency[nLink] = awkward.fromiter( np.array([[-1]*buffLen]*nWafer,np.int16) )
    OverflowCount[nLink] = np.array([0]*nWafer,np.int32)
    ResetCount[nLink] = np.array([0]*nWafer,np.int32)

    Buffer[nLink].stops = Buffer[nLink].starts.copy()
    Buffer_Latency[nLink].stops = Buffer_Latency[nLink].starts.copy()

    latencyTracker[nLink] =  awkward.fromiter( np.array([[-1]*N]*nWafer,np.int16) )
    countBX[nLink] = latencyTracker[nLink].starts.copy()

elapsed = time.time()-tstart
print(elapsed)

shuffleTime = 1000

tstart = time.time()
iBX = 0
for BX in range(1,N+1):
    if BX%shuffleTime == (shuffleTime-1):
        iBX = 0
        for i in range(nWafer):
            np.random.shuffle(waferOccup[i])
            
    iBX += 1
            
    _nTC = waferOccup[:,iBX]
    _nTC += addTC

    if BX%fill==0:
        if fillMod:
            _nTC = np.array([fillSize]*nWafer,np.int16)
    if BX%fill==fillAgain:
        if fillAgainMod:
            _nTC = np.array([fillSize]*nWafer,np.int16)

    _nBits = getBits(_nTC,wordSize)
    tcCount.append(_nTC.copy())
    bitCount.append(_nBits.copy())
    
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

        OverflowCount[nLink] += Buffer_Latency[nLink].content[Buffer_Latency[nLink].starts]>=overflowThreshold

        Buffer[nLink].content[ Buffer[nLink].stops ] += _nBits
        Buffer[nLink].stops += 1

        Buffer_Latency[nLink].content[Buffer_Latency[nLink].stops] = 0
        Buffer_Latency[nLink].stops += 1
        Buffer_Latency[nLink].content[Buffer_Latency[nLink].content>-1] += 1
        
        #find buffers already over buffLen-2
        bufferFull = Buffer_Latency[nLink].content[Buffer_Latency[nLink].starts]>(buffLen-2)
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

meanOcc = np.array([waferOccup[i].mean() for i in range(nWafer)])
meanBits = getBits(meanOcc,-1)


BufferOutput = "Output_Buffer_TTbar200PU_thr%s_layer%i_job%i.pkl"%(("%.1f"%threshold).replace('.','p'),layer,job)
if addTC>0:
    BufferOutput = BufferOutput.replace('_job','_add%iTC_job'%addTC)
if fillMod:
    BufferOutput = BufferOutput.replace('_job','_fill%iEvery%i_job'%(fillSize,fill))
    if fillAgainMod:
        BufferOutput = BufferOutput.replace('_job','_againEvery%i_job'%fillAgain)
if wordSize>0:
    BufferOutput = BufferOutput.replace('_job','_words%i_job'%wordSize)
if args.newFormat:
    BufferOutput = BufferOutput.replace('_job','_newDataFormat_job')
if args.newFormatNoHeader:
    BufferOutput = BufferOutput.replace('_job','_newDataFormatNoHeader_job')
if args.newFormatPaul:
    BufferOutput = BufferOutput.replace('_job','_newDataFormatFromPaul_job')

with open(BufferOutput,'wb') as f:
    pickle.dump(tcCount,f)
    pickle.dump(bitCount,f)
    pickle.dump(latencyTracker,f)
    pickle.dump(OverflowCount,f)
    pickle.dump(ResetCount,f)
    pickle.dump(waferList,f)
    pickle.dump(meanOcc,f)
    pickle.dump(meanBits,f)

