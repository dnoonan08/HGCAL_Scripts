import argparse
import pandas as pd
import numpy as np
import os
import re

def toBin(x,nDigits=5):
    return format(x, '#0%ib'%(nDigits+2))[2:]
toBin = np.vectorize(toBin)

def replaceHeader(val, header, nDigits=5):
    return header + val[nDigits:]
replaceHeader = np.vectorize(replaceHeader)

def shiftBits(x,nDigits=11):
    return x<<nDigits
shiftBits= np.vectorize(shiftBits)


def writeOutputData(inputData, outputDir, SamplingOrder, headerValues, TX_SYNC_WORD, smallOutput):

    if not smallOutput:
        #read in data
        InputEPORTRX = pd.read_csv(f'{inputData}/EPORTRX_output.csv',index_col='entry')
        for c in InputEPORTRX.columns:
            InputEPORTRX[c] = InputEPORTRX[c].apply(int,args=(2,))
    
        InputEPORTRX.astype(int).loc[SamplingOrder].to_csv(f'{outputDir}/MuxFixCalib_Input_ePortRX.csv',index=False)
        del InputEPORTRX
    
        InputF2F = pd.read_csv(f'{inputData}/F2F_output.csv',index_col='entry')
        InputF2F.loc[SamplingOrder].to_csv(f'{outputDir}/MuxFixCalib_PreCalibration_F2F.csv',index=False)
        del InputF2F
    
        InputCALQ = pd.read_csv(f'{inputData}/CALQ_output.csv',index_col='entry')
        InputCALQ.loc[SamplingOrder].to_csv(f'{outputDir}/Algorithm_Input_CalQ.csv',index=False)
        del InputCALQ
    
        Threshold_AddMap = pd.read_csv(f'{inputData}/threshold_address.csv',index_col='entry')
        Threshold_AddMap.loc[SamplingOrder].to_csv(f'{outputDir}/Algorithm_Output_AddrMap.csv',index=False)
        del Threshold_AddMap
    
        Threshold_Charge = pd.read_csv(f'{inputData}/threshold_charge.csv',index_col='entry')
        Threshold_Charge.loc[SamplingOrder].to_csv(f'{outputDir}/Algorithm_Output_ChargeQ.csv',index=False)
        del Threshold_Charge
    
        Threshold_Wafer = pd.read_csv(f'{inputData}/threshold_wafer.csv',index_col='entry')
        Threshold_Wafer.loc[SamplingOrder,["NTCQ"]].to_csv(f'{outputDir}/Algorithm_Output_NTCQ.csv',index=False)
        Threshold_Wafer.loc[SamplingOrder,["SUM_NOT_TRANSMITTED"]].to_csv(f'{outputDir}/Algorithm_Output_Sum.csv',index=False,header='SUM')
        del Threshold_Wafer
    
        BC_Addresses = pd.read_csv(f'{inputData}/bc_address.csv',index_col='entry')
        BC_Addresses.columns = [f'BC_TC_MAP_{i}' for i in range(48)]
        BC_Addresses.loc[SamplingOrder].to_csv(f'{outputDir}/Algorithm_Output_BC_TC_map.csv',index=False)
        del BC_Addresses
    
        BC_Charge = pd.read_csv(f'{inputData}/bc_charge.csv',index_col='entry')
        BC_Charge.loc[SamplingOrder].to_csv(f'{outputDir}/Algorithm_Output_BC_Charge.csv',index=False,header=[f'BC_CHARGE_{i}' for i in range(48)])
        del BC_Charge
    
        Cols_XTC4_9 = [f'XTC4_9_{i}' for i in range(12)]
        Cols_XTC4_7 = [f'XTC4_7_{i}' for i in range(12)]
        Cols_MAX4_ADDR = [f'MAX4_ADDR_{i}' for i in range(12)]
        Cols_XTC16_9 = [f'XTC16_9_{i}' for i in range(3)]
        Cols_MAX16_ADDR = [f'MAX16_ADDR_{i}' for i in range(3)]
    
        STC_Sum = pd.read_csv(f'{inputData}/stc_sum.csv',index_col='entry')
        STC_Sum.columns = Cols_XTC4_9 + Cols_XTC16_9 + Cols_XTC4_7
        STC_Sum.loc[SamplingOrder,Cols_XTC4_7].to_csv(f'{outputDir}/Algorithm_Output_XTC4_7.csv',index=False)
        STC_Sum.loc[SamplingOrder,Cols_XTC4_9].to_csv(f'{outputDir}/Algorithm_Output_XTC4_9.csv',index=False)
        STC_Sum['XTC16_9'] = shiftBits(STC_Sum.XTC16_9_0,0) + shiftBits(STC_Sum.XTC16_9_1,9) + shiftBits(STC_Sum.XTC16_9_2,16)
        STC_Sum.loc[SamplingOrder,['XTC16_9']].to_csv(f'{outputDir}/Algorithm_Output_XTC16_9.csv',index=False)
        del STC_Sum
    
        STC_Idx = pd.read_csv(f'{inputData}/stc_idx.csv',index_col='entry')
        STC_Idx.columns = Cols_MAX4_ADDR + Cols_MAX16_ADDR
        STC_Idx.loc[SamplingOrder,Cols_MAX4_ADDR].to_csv(f'{outputDir}/Algorithm_Output_MAX4_ADDR.csv',index=False)
        STC_Idx.loc[SamplingOrder,Cols_MAX16_ADDR].to_csv(f'{outputDir}/Algorithm_Output_MAX16_ADDR.csv',index=False)
        del STC_Idx
    
        RPT_Chg = pd.read_csv(f'{inputData}/repeater_charge.csv',index_col='entry')
        RPT_Chg.loc[SamplingOrder].to_csv(f'{outputDir}/Algorithm_Output_RepeaterQ.csv',index=False)
        del RPT_Chg
    
        print('   Starting formatting')
        head = toBin(headerValues)
    
        cols =[f'WORD_{i}' for i in range(28)]
    
        Threshold_Format = pd.read_csv(f'{inputData}/threshold_formatblock.csv',index_col='entry',dtype=np.object)
        sampledThresholdFormat = Threshold_Format.loc[SamplingOrder]
        sampledThresholdFormat['FRAMEQ'] = replaceHeader(sampledThresholdFormat.FRAMEQ,head)
        sampledThresholdFormat['FRAMEQ_TRUNC'] = replaceHeader(sampledThresholdFormat.FRAMEQ_TRUNC,head)
        sampledThresholdFormat['WORD_0'] = replaceHeader(sampledThresholdFormat.WORD_0,head)
        sampledThresholdFormat['NULLVAL'] = shiftBits(headerValues) + int(TX_SYNC_WORD)
        sampledThresholdFormat['wordCount'] = sampledThresholdFormat[cols].notna().sum(axis=1)
        sampledThresholdFormat['FRAMEQ_TRUNC'] = sampledThresholdFormat['FRAMEQ_TRUNC'].apply(int,args=(2,))
        for c in cols:
            sampledThresholdFormat.loc[sampledThresholdFormat[c].notna(),c] = sampledThresholdFormat.loc[sampledThresholdFormat[c].notna(),c].apply(int,args=(2,))
            sampledThresholdFormat.loc[sampledThresholdFormat[c].isna(),c] = sampledThresholdFormat.loc[sampledThresholdFormat[c].isna(),'NULLVAL']        
        sampledThresholdFormat[cols+['FRAMEQ_TRUNC','wordCount']].to_csv(f'{outputDir}/Formatter_Output_ThresholdSum.csv',index=False)
        del sampledThresholdFormat
        del Threshold_Format
        
        BC_Format = pd.read_csv(f'{inputData}/bc_formatblock.csv',index_col='entry',dtype=np.object)
        sampledBCFormat = BC_Format.loc[SamplingOrder]
        sampledBCFormat['FRAMEQ'] = replaceHeader(sampledBCFormat.FRAMEQ,head)
        sampledBCFormat['WORD_0'] = replaceHeader(sampledBCFormat.WORD_0,head)
        sampledBCFormat['NULLVAL'] = shiftBits(headerValues) + int(TX_SYNC_WORD)
        sampledBCFormat['wordCount'] = sampledBCFormat[cols].notna().sum(axis=1)
        for c in cols:
            sampledBCFormat.loc[sampledBCFormat[c].notna(),c] = sampledBCFormat.loc[sampledBCFormat[c].notna(),c].apply(int,args=(2,))
            sampledBCFormat.loc[sampledBCFormat[c].isna(),c] = sampledBCFormat.loc[sampledBCFormat[c].isna(),'NULLVAL']        
        sampledBCFormat[cols + ['wordCount']].to_csv(f'{outputDir}/Formatter_Output_BC.csv',index=False)
        del sampledBCFormat
        del BC_Format
    
        STC_Format = pd.read_csv(f'{inputData}/stc_formatblock.csv',index_col='entry',dtype=np.object)
        sampledSTCFormat = STC_Format.loc[SamplingOrder]
        sampledSTCFormat['FRAMEQ'] = replaceHeader(sampledSTCFormat.FRAMEQ,head)
        sampledSTCFormat['WORD_0'] = replaceHeader(sampledSTCFormat.WORD_0,head)
        sampledSTCFormat['NULLVAL'] = shiftBits(headerValues) + int(TX_SYNC_WORD)
        sampledSTCFormat['wordCount'] = sampledSTCFormat[cols].notna().sum(axis=1)
        for c in cols:
            sampledSTCFormat.loc[sampledSTCFormat[c].notna(),c] = sampledSTCFormat.loc[sampledSTCFormat[c].notna(),c].apply(int,args=(2,))
            sampledSTCFormat.loc[sampledSTCFormat[c].isna(),c] = sampledSTCFormat.loc[sampledSTCFormat[c].isna(),'NULLVAL']        
        sampledSTCFormat[cols + ['wordCount']].to_csv(f'{outputDir}/Formatter_Output_STC.csv',index=False)
        del sampledSTCFormat
        del STC_Format
    
        RPT_Format = pd.read_csv(f'{inputData}/repeater_formatblock.csv',index_col='entry',dtype=np.object)
        sampledRPTFormat = RPT_Format.loc[SamplingOrder]
        sampledRPTFormat['FRAMEQ'] = replaceHeader(sampledRPTFormat.FRAMEQ,head)
        sampledRPTFormat['WORD_0'] = replaceHeader(sampledRPTFormat.WORD_0,head)
        sampledRPTFormat['NULLVAL'] = shiftBits(headerValues) + int(TX_SYNC_WORD)
        sampledRPTFormat['wordCount'] = sampledRPTFormat[cols].notna().sum(axis=1)
        for c in cols:
            sampledRPTFormat.loc[sampledRPTFormat[c].notna(),c] = sampledRPTFormat.loc[sampledRPTFormat[c].notna(),c].apply(int,args=(2,))
            sampledRPTFormat.loc[sampledRPTFormat[c].isna(),c] = sampledRPTFormat.loc[sampledRPTFormat[c].isna(),'NULLVAL']        
        sampledRPTFormat[cols + ['wordCount']].to_csv(f'{outputDir}/Formatter_Output_RPT.csv',index=False)
        del sampledRPTFormat
        del RPT_Format
    else:
        Threshold_Wafer = pd.read_csv(f'{inputData}/threshold_wafer.csv',index_col='entry')
        Threshold_Wafer.loc[SamplingOrder,["NTCQ"]].to_csv(f'{outputDir}/Algorithm_Output_NTCQ.csv',index=False)
        del Threshold_Wafer

        print('   Starting formatting')
        head = toBin(headerValues)
        cols =[f'WORD_{i}' for i in range(28)]

        Threshold_Format = pd.read_csv(f'{inputData}/threshold_formatblock.csv',index_col='entry',dtype=np.object)
        sampledThresholdFormat = Threshold_Format.loc[SamplingOrder]
        sampledThresholdFormat['FRAMEQ'] = replaceHeader(sampledThresholdFormat.FRAMEQ,head)
        sampledThresholdFormat['FRAMEQ_TRUNC'] = replaceHeader(sampledThresholdFormat.FRAMEQ_TRUNC,head)
        sampledThresholdFormat['WORD_0'] = replaceHeader(sampledThresholdFormat.WORD_0,head)
        sampledThresholdFormat['NULLVAL'] = shiftBits(headerValues) + int(TX_SYNC_WORD)
        sampledThresholdFormat['wordCount'] = sampledThresholdFormat[cols].notna().sum(axis=1)
        sampledThresholdFormat['FRAMEQ_TRUNC'] = sampledThresholdFormat['FRAMEQ_TRUNC'].apply(int,args=(2,))
        for c in cols:
            sampledThresholdFormat.loc[sampledThresholdFormat[c].notna(),c] = sampledThresholdFormat.loc[sampledThresholdFormat[c].notna(),c].apply(int,args=(2,))
            sampledThresholdFormat.loc[sampledThresholdFormat[c].isna(),c] = sampledThresholdFormat.loc[sampledThresholdFormat[c].isna(),'NULLVAL']        
        sampledThresholdFormat[cols+['FRAMEQ_TRUNC','wordCount']].to_csv(f'{outputDir}/Formatter_Output_ThresholdSum.csv',index=False)
        del sampledThresholdFormat
        del Threshold_Format


def sampleData(directoryName="./", outputName="outputTest", nBX=5000, samplingList=None, extraName=None, smallOutput=False):

    waferInfo = re.search('\w+_D(\d+)L(\d+)W(\d+)',directoryName)
    if waferInfo is None:
        waferInfo = re.search('\w+_D(\d+)L(\d+)U(\d+)V(\d+)',directoryName)
        

    if not waferInfo is None:
        _subdetector = int(waferInfo.groups()[0])
        _layer = int(waferInfo.groups()[1])
        if len(waferInfo.groups())==3:
            _wafer = int(waferInfo.groups()[2])
        else:
            _waferU = int(waferInfo.groups()[2])
            _waferV = int(waferInfo.groups()[3])
            _wafer=100*_waferU+_waferV
    else:
        _subdetector = 0
        _layer=0
        _wafer=0


    wafer = f'D{_subdetector}L{_layer}W{_wafer}'
    inputData = directoryName

    TX_SYNC_WORD='00000000000'




    #make output directory
    outputDir = f'{outputName}_{wafer}'
    if extraName is not None:
        outputDir = f'{outputName}_{wafer}_{extraName}'
    os.makedirs(outputDir,exist_ok=True)
    os.makedirs(f'{outputDir}/Idle',exist_ok=True)

    InputCALQ = pd.read_csv(f'{inputData}/CALQ_output.csv',index_col='entry')

    #get sampling order
    if samplingList is None:
        SamplingOrder = []    
        SamplingOrder = np.random.choice(InputCALQ.index.values,nBX)
    else:
        SamplingOrder = pd.read_csv(samplingList)['entry'].values
#        SamplingOrder = np.array(samplingLinst)
        nBX = len(SamplingOrder)
    print(SamplingOrder)

    del InputCALQ

    np.savetxt(f'{outputDir}/EventSampleOrder.csv', SamplingOrder,fmt='%d', header='entry',comments="")
    headerValues = np.arange(len(SamplingOrder))%32

    # output basic wafer info
    with open(f'{outputDir}/WaferInfo.txt','w') as _file:
        _file.write(f'InputData {inputData}\n')
        _file.write(f'Subdet {_subdetector}\n')
        _file.write(f'Layer {_layer}\n')
        _file.write(f'Wafer {_wafer}\n')


    if not smallOutput:
        #### Track Registers
        print('getting registers')
    
        threshReg = pd.read_csv(f'{inputData}/thres_{wafer}.csv')
        threshReg = threshReg.loc[threshReg.index.repeat(nBX)]
        threshReg.columns = [f'THRESHV_{i}' for i in range(48)]
        threshReg.to_csv(f'{outputDir}/Algorithm_Input_Threshold.csv',index=False)
        del threshReg
    
        calibReg = pd.read_csv(f'{inputData}/calib_{wafer}.csv')
        calibReg = calibReg.loc[calibReg.index.repeat(nBX)]
        calibReg.columns = [f'CALV_{i}' for i in range(48)]
        calibReg = (calibReg*2**11).astype(int)
        calibReg.to_csv(f'{outputDir}/Calibration_Input_Calibration.csv',index=False)
        del calibReg
    
        registers = pd.read_csv(f'{inputData}/registers_{wafer}.csv')
    #    isHDM = registers.isHDM.loc[0]
        registers = registers.loc[registers.index.repeat(nBX)]
        registers.columns = ['HIGH_DENSITY','DROPPED_BITS']
        registers[['HIGH_DENSITY']].astype(int).to_csv(f'{outputDir}/Algorithm_Input_HighDensity.csv',index=False)
        registers[['DROPPED_BITS']].astype(int).to_csv(f'{outputDir}/Algorithm_Input_DroppedBits.csv',index=False)
        del registers
    
        type_DF = pd.DataFrame({'TYPE':[0]*nBX})
        type_DF.to_csv(f'{outputDir}/Algorithm_Input_Type_TS.csv',index=False)
        type_DF['TYPE'] = 1
        type_DF.to_csv(f'{outputDir}/Algorithm_Input_Type_STC.csv',index=False)
        type_DF['TYPE'] = 2
        type_DF.to_csv(f'{outputDir}/Algorithm_Input_Type_BC.csv',index=False)
        type_DF['TYPE'] = 3
        type_DF.to_csv(f'{outputDir}/Algorithm_Input_Type_RPT.csv',index=False)
        del type_DF
    
        #need to track is HDM
        np.savetxt(f'{outputDir}/Algorithm_Input_Header.csv', headerValues,fmt='%d', header='HEADER',comments="")
        np.savetxt(f'{outputDir}/Algorithm_Output_Header.csv', headerValues,fmt='%d', header='HEADER',comments="")
    
    
        pd.DataFrame({'TX_SYNC_WORD':[int(TX_SYNC_WORD)]*nBX}).to_csv(f'{outputDir}/Formatter_Buffer_Input_Tx_Sync_Word.csv',index=False)

    print('Starting to write data')
    writeOutputData(inputData, outputDir, SamplingOrder, headerValues, TX_SYNC_WORD, smallOutput)

    if not smallOutput:
        #write idle word data as well
        print('Starting to write idle data')
        writeOutputData(inputData=f'{inputData}/Idle', outputDir=f'{outputDir}/Idle', SamplingOrder=[0], headerValues=[0], TX_SYNC_WORD=TX_SYNC_WORD)

    return outputDir


import glob
def alterCSVFiles(outputDir):
    csvFiles = glob.glob(f'{outputDir}/*csv')

    for fName in csvFiles:
        with open(fName, "r") as sources:
            lines = sources.readlines()
        with open(fName, "w") as sources:
            for line in lines:
                sources.write(line.replace(',',', '))




if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--inputDir", default = None,dest="inputDir", required=True,help="input directory to sample (required)")
    parser.add_argument('-o',"--output", default = "outputData",dest="outputName", help="file name for output (default: outputData)")
    parser.add_argument('-N',"--nBX", type=int, default = 100,dest="nBX", help="number of BX to simulate (default: 100)")
    parser.add_argument('--samplingList', '--sampleList', default = None,dest="samplingList", help="file containing the list of BX to sample (default: None)")
    parser.add_argument('--extra', '--extraText', default = None,dest="extraText", help="extra text ot append to directory name (ex: version number, other notes to track) (default: None)")
    parser.add_argument('--CSVNoSpace',dest='replaceCSVDelim', action='store_false', default=True, help='keep the output csv files as comma delimited with no trailing space')
    parser.add_argument('--small',dest='smallOutput', action='store_true', default=False, help='keep only a limited set of the outputs to reduce size (threshold formmating)')
    args = parser.parse_args()

    outputDir = sampleData(directoryName = args.inputDir, outputName = args.outputName, nBX = args.nBX, samplingList = args.samplingList, extraName=args.extraText, smallOutput=args.smallOutput)

    if not args.smallOutput:
        alterCSVFiles(outputDir)

        linkFiles = [["Formatter_Buffer_Input_Type_BC.csv",              "Algorithm_Input_Type_BC.csv"],
                     ["Formatter_Buffer_Input_Type_RPT.csv",             "Algorithm_Input_Type_RPT.csv"],
                     ["Formatter_Buffer_Input_Type_STC.csv",             "Algorithm_Input_Type_STC.csv"],
                     ["Formatter_Buffer_Input_Type_TS.csv",              "Algorithm_Input_Type_TS.csv"],
                     ["Formatter_Buffer_Input_ChargeQ.csv",              "Algorithm_Output_ChargeQ.csv"],
                     ["Formatter_Buffer_Input_AddrMap.csv",              "Algorithm_Output_AddrMap.csv"],
                     ["Formatter_Buffer_Input_NTCQ.csv",                 "Algorithm_Output_NTCQ.csv"],
                     ["Formatter_Buffer_Input_Sum.csv",                  "Algorithm_Output_Sum.csv"],
                     ["Formatter_Buffer_Input_XTC4_9.csv",               "Algorithm_Output_XTC4_9.csv"],
                     ["Formatter_Buffer_Input_XTC4_7.csv",               "Algorithm_Output_XTC4_7.csv"],
                     ["Formatter_Buffer_Input_XTC16_9.csv",              "Algorithm_Output_XTC16_9.csv"],
                     ["Formatter_Buffer_Input_MAX4_ADDR.csv",            "Algorithm_Output_MAX4_ADDR.csv"],
                     ["Formatter_Buffer_Input_MAX16_ADDR.csv",           "Algorithm_Output_MAX16_ADDR.csv"],
                     ["Formatter_Buffer_Input_BC_Charge.csv",            "Algorithm_Output_BC_Charge.csv"],
                     ["Formatter_Buffer_Input_BC_TC_map.csv",            "Algorithm_Output_BC_TC_map.csv"],
                     ["Formatter_Buffer_Input_RepeaterQ.csv",            "Algorithm_Output_RepeaterQ.csv"]]
    
        for f1, f2 in linkFiles:
            if not os.path.exists(f'{outputDir}/{f1}'):
                os.symlink(f'{f2}',f'{outputDir}/{f1}')
            
