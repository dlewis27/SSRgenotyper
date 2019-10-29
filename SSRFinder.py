# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 10:39:12 2019

@author: Dan-L
"""

#find fasta repeats
#reverse complelment taken into account (SAM already changes it to reverse comp)
#single letter repeats are removed (i.e. GGGG)
#repeats required to have 2bp on both sides that match reference flank
#lowercase seqs will not be considered
#offset is set to 75


import sys
import re
import itertools
import pandas as pd
from collections import Counter
from Bio.Seq import Seq
import time
from Bio import SeqIO


fastaRef = sys.argv[1]
#mapping done with BWA (no options selected)
#fasta file that made sam (second file) file MUST BE MAPPED to the referece file (first file)
#there are no optical duplicates
#sam files must be filtered for quality and headers removed (samtools view -q 45 fileName)
#then only grab first 10 columns (cut -f1-10 fileName)
#ls *.f1 > files.txt

samsDoc = sys.argv[2]
samFiles = []

with open(samsDoc) as sd:
    input = sd.readlines()
    for i in input:
        samFiles.append(i)

majorMinorRatio = .2
lowerBound = 1 # for left most position
upperBound = 72 #for left most postion (recommmended to be the flanking length)
flankAdded = 75 # how much upstream and downstream added
flankOffset = flankAdded - 30 #for the repeat searching window
minRefFreq = 4 # the minimum frequency of an SSR found in the reference to be reported and have matching reads processed
minSSRfreq = 3 # the minimum number of SSR repeats on a strand needed to be reported (recommended to be 3)
minNumReads = 2 #the minimum number of reads with the SSR frequency needed to be reported
numFlankNucs = 10


stat0 = 0
stat1 = 0
stat1_5 = 0
stat2 = 0
stat3 = 0


refDict = SeqIO.to_dict(SeqIO.parse(fastaRef, "fasta"))


#if flag = 16 then read is reverse com

def prepSam(samFile):
    #samFile is string
    samData = pd.read_csv(samFile, delimiter='\t',encoding='utf-8', header = None)
    samData = samData.drop(samData.columns[[4,5,6,7,8]], axis=1)
    samData.columns = ["ReadName", "Flag", "RefReadName", "LeftMostPos", "ReadSeq"]
    return samData

def genRepeats(repeatSize):
    returnable= []
    for x in itertools.product('ATGC', repeat=repeatSize):
        returnable.append(''.join(x))
    return(list(returnable))

def allCharactersSame(s) : 
    n = len(s) 
    for i in range(1, n) : 
        if s[i] != s[0] : 
            return False
    return True

def getMax(array):
    myMax = 0
    returnable = []
    for subArray in array:
        for v in subArray:
            if (len(v) > 0 and len(v[0])> myMax and allCharactersSame(v[1]) == False):
                myMax = len(v[0])
                returnable = v
    return returnable

    

def findRefRepeat(sequence, repeatSize):
    generatedRepeats = genRepeats(repeatSize) 
    result = []
    for s in generatedRepeats:
        reg = "".join(["((", s, ")+)"])
        #must have same 2bp flanking region on both sides
        found = re.findall(reg, str(sequence))
        result.append(found)
    theMax = getMax(result)
    if (len(theMax) <= 1):
        return(None, None, None, None)
    pattern = theMax[1]
    numRepeats = int(len(theMax[0])/repeatSize)
    dotStr = "."*numFlankNucs
    regPart2 = "".join(["(",dotStr,")((", pattern, "){", str(numRepeats), "})(",dotStr,")"])
    foundWflanks = re.findall(regPart2, str(sequence))
    if(len(foundWflanks) == 0):
        return(None, None, None, None)
    flankL = foundWflanks[0][0]
    flankR = foundWflanks[0][3]
    return(pattern, numRepeats, flankL, flankR)

def getMaxLen(array):
    myMax = 0
    for i in array:
        if (len(i[0])> myMax):
            myMax = len(i[0])
    return myMax
            
def findSpecificRepeat(read, repeat, flankL, flankR):
    reg = "".join([flankL,"((", repeat, ")+)", flankR])
    found = re.findall(reg, str(read))
    theMaxLen = getMaxLen(found)
    return theMaxLen/len(repeat)

def getRefSeqPattern(nameOfRead, lengthOfRepeat):
    refSeq = refDict[nameOfRead].seq
    refSeq = str(refSeq)
    refSeqSub = refSeq[flankOffset : -flankOffset]
    refPattern, refnumRepeats, flankL, flankR = findRefRepeat(refSeqSub, lengthOfRepeat)
    return (refPattern, refnumRepeats, flankL, flankR)


def subsetSam(refReadName, samData):
    subSam = samData.loc[samData["RefReadName"] == refReadName]
    subSam["LeftMostPos"] = subSam["LeftMostPos"].astype(int)
    subSam = subSam.loc[subSam["LeftMostPos"].isin(list(range(lowerBound, upperBound)))]
    return subSam

def veiwRefRead(refName):
    print(refDict[refName].seq)

def viewSubset(refName):
    with open('subOut.txt', 'w') as outFile:
        outFile.write(' ')
        outFile.write(str(refDict[refName].seq))
        outFile.write("\n")
        frame = subsetSam(refName)
        reads = frame["ReadSeq"].values
        leftPos = frame["LeftMostPos"].values
        flags = frame["Flag"].values
        
        i = 0
        while i < len(reads):
            finalRead = reads[i]
            outFile.write(' '* leftPos[i])
            outFile.write(str(finalRead))
            outFile.write("\n")
            i+=1
    
def printResults(resultArray):
    uniqueValues = set(resultArray)
    writeOut = []
    for numRepeats in uniqueValues:
        readsWithNumRepeats = resultArray.count(numRepeats)
        if (readsWithNumRepeats >= minNumReads and numRepeats >= minSSRfreq):
            writeOut.append([numRepeats, readsWithNumRepeats])
    if len(writeOut) == 0:
        global stat0
        stat0 +=1
        return "0-0"
        
    elif len(writeOut) == 1:
        global stat1
        stat1 +=1
        numRepeat = writeOut[0][0]
        return (str(int(numRepeat)) +"-"+str(int(numRepeat)))

    elif len(writeOut) == 2:
        allele1Support = writeOut[0][1]
        allele2Support = writeOut[1][1]
        ratio = allele1Support/allele2Support
        if ratio > 1:
            ratio = 1/ratio
        if ratio >= majorMinorRatio:
            global stat2
            stat2 +=1
            return(str(int(writeOut[0][0])) + "-" + str(int(writeOut[1][0])))
        else:
            global stat1_5
            stat1_5 +=1
            return(str(int(writeOut[0][0])) + "-" + str(int(writeOut[0][0])))

    elif len(writeOut) > 2:
        global stat3
        stat3 +=1
        return "0-0"
    return "ERRROR"
    
def writeStats():
    statOut = open("stat.ssr", "w")
    statOut.write("The following are the stats for the sequences in the reference\n")
    statOut.write("Minimum major/minor allele frequency: " + str(majorMinorRatio) + "\n")
    statOut.write("Minimum reads needed for SSR to be considered: " + str(minNumReads) + "\n")
    statOut.write("Minimum SSR frequency: " + str(minSSRfreq) + "\n")
    statOut.write("Minimum refernce SSR frequency: " + str(minRefFreq) + "\n")
    statOut.write("No call: " + str(stat0) + "\n")
    statOut.write("Homozygous: " + str(stat1) + "\n")
    statOut.write("At least homozygous (minimum major/minor frequency not met) : " + str(stat1_5) + "\n")
    statOut.write("Heterozygous: " + str(stat2) + "\n")
    statOut.write("Ambigous: " + str(stat3) + "\n")
    statOut.close()

def findSamReads(subSam, refInput):
    if refInput == 0:
        return "0-0"
    else:
        refPattern = refInput[0]
        refName = refInput[1]
        flankL =refInput[2]
        flankR = refInput[3]
        samRepeatCounts = []
        for index, row in subSam.iterrows():
            readFromSam = row["ReadSeq"]
            appendable = findSpecificRepeat(readFromSam, refPattern, flankL, flankR)
            samRepeatCounts.append(appendable)
        returnable = printResults(samRepeatCounts)
    return returnable



def processSams(refData, outputDict):
    for samFile in samFiles:
        samFile = samFile.rstrip("\n")
        with open(samFile) as samF:
            samData = prepSam(samF)
            samName = samFile[0:7]
            colToAppend = []
            for refName in  refData.keys():
                subSam = subsetSam(refName, samData)
                #refData[refName] is a list
                appendable = findSamReads(subSam, refData[refName])
                colToAppend.append(appendable)
            outputDict[samName] = colToAppend


def searchRef(refDict, outputDict):
    refData = {}
    for name in refDict:
        refName = name
        refSeq = refDict[refName].seq
        refSeq = str(refSeq)
        refPattern2, refNumRepeats2, flankL2, flankR2 = getRefSeqPattern(refName, 2)
        refPattern3, refNumRepeats3, flankL3, flankR3 = getRefSeqPattern(refName, 3)
        refPattern4, refNumRepeats4, flankL4, flankR4 = getRefSeqPattern(refName, 4)
        
        if (refNumRepeats2 == None):
            refNumRepeats2 = 0
        if (refNumRepeats3 == None):
            refNumRepeats3 = 0
        if (refNumRepeats4 == None):
            refNumRepeats4 = 0
        maxFreqArray = [refNumRepeats2, refNumRepeats3, refNumRepeats4]   
        maxFreq = max(maxFreqArray)
        
        if (maxFreq < minRefFreq):
            refData[refName] = 0
        elif(refNumRepeats2 == maxFreq):
            refData[refName] = [refPattern2, refNumRepeats2, flankL2, flankR2]
        elif(refNumRepeats3 == maxFreq):
            refData[refName] = [refPattern3, refNumRepeats3, flankL3, flankR3]
        elif(refPattern4 == maxFreq):
            refData[refName] = [refPattern4, refNumRepeats4, flankL4, flankR4]
        else:
            refData[refName] = 0
        
    #create first 2 columns in outputDict
    outputDict["RefName"] = []
    for i in refData.keys():
        outputDict["RefName"].append(i)
#    outputDict["RefFreq"] = []
#    for r in refData.values():
#        print(r)
#        outputDict["RefFreq"].append(r[1])
    print("processed reference")
    return refData

def main():
    
    startTime = time.time()      
    
    outputDict = {}
    refData = searchRef(refDict, outputDict)
    #process Sams
    processSams(refData, outputDict)
    outputDf = pd.DataFrame(outputDict)
    outputDf.to_csv("SSRoutput.ssr", sep= "\t")
    
    endTime = time.time()
    writeStats()
    runtime = endTime-startTime
    runtime = str(round(runtime/60, 2))
    print("run time:", runtime, "minutes")
        
main()


