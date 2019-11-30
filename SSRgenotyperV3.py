# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 12:42:10 2019

@author: Dan-L
"""


#find fasta repeats
#reverse complelment taken into account (SAM already changes it to reverse comp)
#single letter repeats are removed (i.e. GGGG)
#lowercase seqs will not be considered


import re
import itertools
import pandas as pd
#from collections import Counter
#from Bio.Seq import Seq
import time
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("ReferenceFile", help = "The reference file (FASTA)")
parser.add_argument("SamFiles", help = "Text document with the SAM file names seperated by newline")
parser.add_argument("OutputFile", help = "Output file name ( will end with \".ssr\")")
parser.add_argument("-A","--AlleleRatio", help = "The minmum ratio of major to minor alleles, between 0 and 1 (default = .2)", type=float, default = .2)
parser.add_argument("-N", "--NameSize", help = "The number of characters to be printed from each SAM file name in the output table (default = 100)", type= int, default = 100)
parser.add_argument("-R", "--RefUnits", help = "The minimum number of SSR units in a reference SSR (default = 4)", type=int, default = 4)
parser.add_argument("-P", "--PopUnits", help = "The minimum number of SSR units in an accession SSR (default = 3)", type=int, default = 3)
parser.add_argument("-F", "--FlankSize", help = "The number of flanking bases on each side of the SSR that must match the reference (default = 15)", type=int, default= 15)
parser.add_argument("-S", "--Support", help = "Then minimum number of supporting reads for alleles to be called (default = 3)", type = int, default = 3)
parser.add_argument("-W", "--WindowOffset", help = "Offset on each side of the reference sequence, making a window for searching for the SSR (default = 1)", type=int, default = 1)
parser.add_argument("-r", "--refFilter", help = "If the porportion of accesions that had no call meet this threshhold, then this marker will not be reported, between 0 and 1 (default = 1)", type=float, default = 1)
parser.add_argument("-Q", "--QualityFilter", help = "Reads with quality score below this level will be filtered out (default = 45)", type=int, default=45)
parser.add_argument("-X", "--Xdebug", help = "Provide marker name and SAM file name seperated by ','. This will also be the output file name (default = '')", type=str, default = "")
parser.add_argument("-M", "--Map", help = "Output a table showing relation to two parents. Make sure the first 2 SAM file names are the parents (default = False)", type=bool, default = False)

args = parser.parse_args()


fastaRef = args.ReferenceFile
#mapping done with BWA (no options selected)
#fasta file that made sam (second file) file MUST BE MAPPED to the reference file (first file)
#there are no optical duplicates
#sam files must be filtered for quality and headers removed (samtools view -q 45 fileName)
#then only grab first 10 columns (cut -f1-10 fileName)
#ls *.f1 > files.txt

samsDoc = args.SamFiles
outFile = args.OutputFile
samFiles = []

with open(samsDoc) as sd:
    input = sd.readlines()
    for i in input:
        samFiles.append(i)

majorMinorRatio = args.AlleleRatio

flankOffset = args.WindowOffset
minRefFreq = args.RefUnits
minSSRfreq = args.PopUnits 
minNumReads = args.Support 
numFlankNucs = args.FlankSize
nameSize = args.NameSize
refFilter = args.refFilter
qualityFilter = args.QualityFilter
debugName = args.Xdebug
doMap = args.Map

noSSRinRef = 0 #-1
noReadsMapped =0 #-2
notEnoughCov = 0 #-4
homo = 0 
hetero = 0
alleleFreqNotMet =0
ambiguous = 0 #-3


refDict = SeqIO.to_dict(SeqIO.parse(fastaRef, "fasta"))

def prepSam(samFile, includeNames = False):
    #samFile is string
    samData = {}
    with open(samFile, 'r') as f:
        for line in f:
            if line.startswith("@") == False:
                splitLine = line.split('\t')
                if len(splitLine) > 9:
                    refName = splitLine[2]
                    matchRead = splitLine[9]
                    quality = int(splitLine[4])
                    readID = splitLine[0]
                    if quality >= qualityFilter:
                        if refName not in samData:
                            samData[refName] = []
                        if includeNames:
                            samData[refName].append([matchRead, readID])
                        else:
                            samData[refName].append(matchRead)
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
    ##filter repeat here
    if theMaxLen/len(repeat) < minSSRfreq:
        return None
    return theMaxLen/len(repeat)

def getRefSeqPattern(nameOfRead, lengthOfRepeat):
    refSeq = refDict[nameOfRead].seq
    refSeq = str(refSeq)
    refSeqSub = refSeq[flankOffset : -flankOffset]
    refPattern, refnumRepeats, flankL, flankR = findRefRepeat(refSeqSub, lengthOfRepeat)
    return (refPattern, refnumRepeats, flankL, flankR)


def veiwRefRead(refName):
    print(refDict[refName].seq)

    
def printResults(resultArray):
    if (len(resultArray) < minNumReads) or (len(resultArray) <= 0):
        #also not enough SSR freq
        global notEnoughCov
        notEnoughCov += 1
        return "0,-4"
    uniqueValues = {}
    for i in resultArray:
        if int(i) not in uniqueValues:
            uniqueValues[int(i)] = 1 
        else:
            uniqueValues[int(i)] += 1
    alleleData = sorted(uniqueValues.items(), key=lambda x: x[1], reverse = True)
    if len(alleleData) == 0:
        return "ERROR"
    
    elif len(alleleData) == 1:
        global homo
        homo +=1
        return (str(alleleData[0][0]) +","+ str(alleleData[0][0]))

    elif len(uniqueValues) == 2:
        allele1Support = alleleData[0][1]
        allele2Support = alleleData[1][1]
        ratio = allele2Support/allele1Support
        if ratio >= majorMinorRatio:
            global hetero
            hetero +=1
            return(str(alleleData[0][0]) + "," + str(alleleData[1][0]))
        else:
            global alleleFreqNotMet
            #reported as homo in table and in stats
            alleleFreqNotMet +=1
            homo +=1
            return(str(alleleData[0][0]) + "," + str(alleleData[0][0]))

    elif len(uniqueValues) > 2:
        global ambiguous
        ambiguous +=1
        return "0,-3"
    return "ERRROR"
    
def writeStats(runtime, refProcessTime):
    statOut = open(outFile + ".ssrstat", "w")
    statOut.write("The following are the stats for the sequences in the reference\n")
    statOut.write("Reference processed in: " + refProcessTime + "\n")
    statOut.write("Minimum major/minor allele frequency: " + str(majorMinorRatio) + "\n")
    statOut.write("Minimum number reads needed for SSR to be considered: " + str(minNumReads) + "\n")
    statOut.write("Minimum SSR unit frequency for population: " + str(minSSRfreq) + "\n")
    statOut.write("Minimum reference SSR unit frequency: " + str(minRefFreq) + "\n")
    statOut.write("No SSR found in reference (counted for each SAM file): " + str(noSSRinRef) + "\n")
    statOut.write("No reads mapped to marker: " + str(noReadsMapped) + "\n")
    statOut.write("Not enough coverage to call: " + str(notEnoughCov) + "\n")
    statOut.write("Minimum allele ratio not met, reported as homozygote in table and in stats: " + str(alleleFreqNotMet) + "\n")
    statOut.write("Homozygous: " + str(homo) + "\n")
    statOut.write("Heterozygous: " + str(hetero) + "\n")
    statOut.write("Ambigous (more than 2 alleles found): " + str(ambiguous) + "\n")
    statOut.write("RunTime: " + str(runtime) + " minutes" + "\n")
    statOut.close()

def findSamReads(subSam, refInput):
    refPattern = refInput[0]
    #refName = refInput[1]
    flankL =refInput[2]
    flankR = refInput[3]
    samRepeatCounts = []
    for r in subSam:
        appendable = findSpecificRepeat(r, refPattern, flankL, flankR)
        if appendable != None:
            samRepeatCounts.append(appendable)
    returnable = printResults(samRepeatCounts)
    return returnable



def processSams(refData, outputDict, inSamFiles):
    #inSamFiles is list
    #split inSamFiles accross nodes
    for samFile in inSamFiles:
        samFile = samFile.rstrip("\n")
        samData = prepSam(samFile)
        samName = samFile[0:nameSize]
        colToAppend = []
        for refName in  refData.keys():
            if refData[refName] == 0:
                appendable= "0,-1"
                global noSSRinRef
                noSSRinRef +=1
            elif refName not in samData:
                appendable = "0,-2"
                global noReadsMapped
                noReadsMapped +=1
            else:
                subSam = samData[refName]
                #refData[refName] is a list
                appendable = findSamReads(subSam, refData[refName])
            colToAppend.append(appendable)
        outputDict[samName] = colToAppend
        print("processed SAM file:", samFile)


def searchRef(refDict, outputDict):
    refData = {}
    for name in refDict:
        refName = name.rstrip()
        refSeq = refDict[refName].seq
        refSeq = str(refSeq)
        refPattern2, refNumRepeats2, flankL2, flankR2 = getRefSeqPattern(refName, 2)
        refPattern3, refNumRepeats3, flankL3, flankR3 = getRefSeqPattern(refName, 3)
        
        if (refNumRepeats2 == None):
            refNumRepeats2 = 0
        if (refNumRepeats3 == None):
            refNumRepeats3 = 0
        maxFreqArray = [refNumRepeats2, refNumRepeats3]   
        maxFreq = max(maxFreqArray)
        
        if (maxFreq < minRefFreq):
            refData[refName] = 0
        elif(refNumRepeats2 == maxFreq):
            refData[refName] = [refPattern2, refNumRepeats2, flankL2, flankR2]
        elif(refNumRepeats3 == maxFreq):
            refData[refName] = [refPattern3, refNumRepeats3, flankL3, flankR3]
        else:
            refData[refName] = 0
        
    #create first 2 columns in outputDict
    outputDict["RefName"] = []
    outputDict["SSRPattern"] = []
    outputDict["RefGenotype"] = []
    for i in refData:
        outputDict["RefName"].append(i)
        if refData[i] != 0:
            outputDict["RefGenotype"].append(refData[i][1])
            outputDict["SSRPattern"].append(refData[i][0])
        else:
            outputDict["RefGenotype"].append(0)
            outputDict["SSRPattern"].append("---")
    return refData

def filterTable(outputDf):
    indexToRemove = []
    for index, row in outputDf.iterrows():
        toBeScored = row[3:]
        score = 0
        for i in toBeScored:
            if i.startswith("0",0,1):
                score +=1
        if score/len(toBeScored) >= refFilter:
            indexToRemove.append(index)
    newTable = outputDf.drop(outputDf.index[indexToRemove])
    return newTable
 
def debug(debugName):
    output="ERROR"
    pattern = ""
    split = debugName.split(",")
    markerName = split[0]
    samFileName = split[1]
    outputDict = {}
    refData = searchRef(refDict, outputDict)
    if refData[markerName] != 0:
        pattern = refData[markerName][0]
    samFile = samFileName.rstrip("\n")
    samData = prepSam(samFile, True)
    if markerName not in refData:
        output = "marker not in reference"
    elif markerName not in samData:
        output = "no reads mapped to this marker"
    else:
        output = "Marker: " + markerName + "\n" + "SAM file: " + samFileName + "\n" + "Pattern: " + pattern + "\n"
        refSeq = str(refDict[markerName].seq)
        flankRlocationRef = re.search(refData[markerName][3], refSeq)
        flankLlocationRef = re.search(refData[markerName][2], refSeq)
        refSeq2 = refSeq[:flankLlocationRef.start() + numFlankNucs] + "   " + refSeq[flankLlocationRef.start() + numFlankNucs:flankRlocationRef.start()] + "   " + refSeq[flankRlocationRef.start():]
        output += "\nReference Sequence: " + refSeq2 + "\n\n" 
        subSam = samData[markerName]
        goodReads = []
        badReads = []
        for sub in subSam:
            s = sub[0]
            subName = sub[1]
            flankRlocation = re.search(refData[markerName][3], s)
            flankLlocation = re.search(refData[markerName][2], s)
            if flankLlocation == None or flankRlocation == None:
                badReads.append(s + "  ---" + subName + "---")
            else:
                offset = " " * (flankLlocationRef.start() - flankLlocation.start() + 20) 
                newS = offset + s[:flankLlocation.start() + numFlankNucs] + "   " + s[flankLlocation.start() + numFlankNucs:flankRlocation.start()] + "  "+ s[flankRlocation.start():]
                goodReads.append(newS + "  ---" + subName + "---")
        for g in goodReads:
            g2 = g.rstrip("\n")
            output += (g2 +"\n")
        output += "\n"
        for b in badReads:
            b2 = b.rstrip("\n")
            output += (b2 + "\n")
    with open("debug.txt", "w" ) as w:
        w.write(output)
        
def makeMap(outputDf):
    #make inferences if a parent is missing in separate program? 
    #iterate through row
    if outputDf.shape[1] < 5:
        print("ERROR: provide atleast 2 SAM files (2 parents)")
        return
    newTableAsList = []
    print("making map")
    for index, row in outputDf.iterrows():
        newRow = []
        parent1 = row[3]
        parent2 = row[4]
        p1 = 'A'
        p2 = 'B'
        if ',' in parent1:
            continue
        if ',' in parent2:
            continue
        #skip non informative markers, just where parent1 and parent2 are diff and not hetero
        if parent1 == parent2:
            continue
        if 0 in parent1:
            p1 = 'U'
        if 0 in parent2:
            p2 = 'U'
        newRow.append(row[0])
        newRow.append(p1)
        newRow.append(p2)
        for r in row[5:]:
            if r == parent1:
                newRow.append(p1)
            elif r == parent2:
                newRow.append(p2)
            elif ',' in r:
                # if partent hetero then won't match anyways
                rSplit = r.split(',')
                allele1 = rSplit[0]
                allele2 = rSplit[1]
                if allele1 == parent1 and allele2 == parent2:
                    newRow.append('H')
                elif allele1 == parent2 and allele2 == parent1:
                    newRow.append('H')
            else:
                newRow.append('U')
        newTableAsList.append(newRow)        
    #transform newTableAsList to newTable
    newDf = pd.DataFrame(newTableAsList)
        #get headers from old table
    headers = outputDf.columns.values.tolist()
    newDf.columns = headers[3:]
    
    newDf.to_csv(outFile + ".map", sep= "\t")
    # headers and stuff to newDf and some other formating

def main():
    if debugName != "":
        debug(debugName)
        return
    
    startTime = time.time()      
    
    outputDict = {}
    refData = searchRef(refDict, outputDict)
    refProcessTime = str(round((time.time() - startTime)/60, 2))
    print("processed ref in:", refProcessTime, "minutes")
    #process Sams
    processSams(refData, outputDict, samFiles)
    outputDf = pd.DataFrame(outputDict)
    if refFilter != 1:
        outputDf = filterTable(outputDf)
    outputDf.to_csv(outFile + ".ssr", sep= "\t")
    
    if doMap == True:
        makeMap(outputDf)
    
    endTime = time.time()
    runtime = endTime-startTime
    runtime = str(round(runtime/60, 2))
    print("run time:", runtime, "minutes")
    writeStats(runtime, refProcessTime)

        
main()


