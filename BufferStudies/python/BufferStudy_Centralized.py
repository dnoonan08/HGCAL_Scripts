#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import awkward
import time
import pickle


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--threshold','--thr', default=1.5,choices=[0.5,1.0,1.5,2.0,2.5],type=float)
parser.add_argument('-N', default=50000, type=int)
parser.add_argument('--layer', default=5, type=int)
parser.add_argument('--job', default=0, type=int)
parser.add_argument('--fill', default=-1, type=int)
parser.add_argument('--fillAgain', default=-1, type=int)
parser.add_argument('--fillFull', default=False, action='store_true')
parser.add_argument('--add', default=0, type=int)
args = parser.parse_args()


threshold = args.threshold
layer = args.layer
N = args.N
fill = args.fill
fillAgain = args.fillAgain
addTC = args.add
job=args.job


fillMod = fill>0
fillAgainMod = fillAgain>0

if not fillMod:
    fill = N + 10
if not fillAgainMod:
    fillAgain = N + 10

thrStr = ("%.1f"%threshold).replace('.','p')

if not os.path.exists('mbOccupancy_thr%s.pkl'%thrStr):
    cmd = 'xrdcp root://cmseos.fnal.gov//store/user/dnoonan/HGCAL_Concentrator/mbOccupancy_thr%s.pkl .'%thrStr
    stdout,stderr = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()
if not os.path.exists('mbOccupancyDropOne_thr%s.pkl'%thrStr):
    cmd = 'xrdcp root://cmseos.fnal.gov//store/user/dnoonan/HGCAL_Concentrator/mbOccupancyDropOne_thr%s.pkl .'%thrStr
    stdout,stderr = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE).communicate()

with open('mbOccupancy_thr%s.pkl'%thrStr,'rb') as _file:
    mbOccupancy = pickle.load(_file,encoding='bytes')
    
with open('mbOccupancyDropOne_thr%s.pkl'%thrStr,'rb') as _file:
    mbOccupancyDropOne = pickle.load(_file,encoding='bytes')


nMB = 18 if layer==5 else 24
motherboards = np.array(range(1,nMB+1))
nTC = np.array([mbOccupancy[layer][i] for i in motherboards],np.int16)
nTCDropOne = np.array([mbOccupancy[layer][i] for i in motherboards],np.int16)



for i in range(nMB):
    np.random.shuffle(nTC[i])
    np.random.shuffle(nTCDropOne[i])
    


from numba import vectorize, int64, float64

@vectorize([int64(int64),float64(float64)])
def getBits(nTC):
    header   = 4 + 5 + 1 ### header(4) + bxID(5) + addressMode(1)                                                                                           
    towerSum = 24
    if nTC < 18:
        total = header + towerSum + 5 + 15*nTC
    else:
        total = header + towerSum + 144 + 7*nTC

    return total


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

#links = [x*.125 for x in range(8,25)]
links = range(1,22)


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
    Buffer[nLink] = awkward.fromiter( np.array([[0]*buffLen]*nMB,np.int16) )
    Buffer_Latency[nLink] = awkward.fromiter( np.array([[-1]*buffLen]*nMB,np.int16) )
    OverflowCount[nLink] = np.array([0]*nMB,np.int32)
    ResetCount[nLink] = np.array([0]*nMB,np.int32)

    Buffer[nLink].stops = Buffer[nLink].starts.copy()
    Buffer_Latency[nLink].stops = Buffer_Latency[nLink].starts.copy()

    latencyTracker[nLink] =  awkward.fromiter( np.array([[-1]*N]*nMB,np.int16) )
    countBX[nLink] = latencyTracker[nLink].starts.copy()

elapsed = time.time()-tstart
print(elapsed)


tstart = time.time()
iBX = 0
for BX in range(1,N+1):
    if BX%180000 == 179999:
        iBX = 0
        for i in range(nMB):
            np.random.shuffle(nTC[i])
            np.random.shuffle(nTCDropOne[i])
            
    iBX += 1
            
    _nTC = nTC[:,iBX]
    _nTC += addTC

    if BX%fill==0:
        if fillMod:
            if args.fillFull:
                _nTC = 48*3
            else:
                _nTC = nTCDropOne[:,iBX] + 48
    if BX%fill==fillAgain:
        if fillAgainMod:
            if args.fillFull:
                _nTC = 48*3
            else:
                _nTC = nTCDropOne[:,iBX] + 48

    _nBits = getBits(_nTC)
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

meanOcc = np.array([nTC[i].mean() for i in range(nMB)])
meanBits = getBits(meanOcc)


BufferOutput = "Output_Buffer_MotherBoard_thr%s_layer%i_job%i.pkl"%(("%.1f"%threshold).replace('.','p'),layer,job)
if addTC>0:
    BufferOutput = BufferOutput.replace('_job','_add%iTC_job'%addTC)
if fillMod:
    BufferOutput = BufferOutput.replace('_job','_fillEvery%i_job'%fill)
    if fillAgainMod:
        BufferOutput = BufferOutput.replace('_job','_againEvery%i_job'%fillAgain)

with open(BufferOutput,'wb') as f:
    pickle.dump(tcCount,f)
    pickle.dump(bitCount,f)
    pickle.dump(latencyTracker,f)
    pickle.dump(OverflowCount,f)
    pickle.dump(ResetCount,f)
    pickle.dump(motherboards,f)
    pickle.dump(meanOcc,f)
    pickle.dump(meanBits,f)


OverflowCount


ResetCount




