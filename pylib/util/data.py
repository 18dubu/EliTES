__author__ = 'mahandong'

import os
import re
from file import *
from log import strd_logger
import numpy


def combineFiles(listOfPathAndFiles,outFilePathAndName, combine=1, header=0):
    try:
        if os.path.exists(outFilePathAndName):
            os.remove(outFilePathAndName)
            print ("file removed:%s" % outFilePathAndName)
        fhOut = open(outFilePathAndName,'w')
    except IOError:
        return False
    allSet = set()
    allList = []
    for file in listOfPathAndFiles:
        if os.path.exists(file):
            try:
                fhIn = open(file, 'r')
                if combine == 1:
                    if header == 1:
                        fhIn.readline()
                        header = 0
                    for lines in fhIn:
                        allSet.add(lines.strip())
                if combine == 0:
                    if header == 1:
                        fhIn.readline()
                        header=0
                    for lines in fhIn:
                        allList.append(lines.strip())
                elif combine != 1 & combine != 0:
                    print ("combine or not, input:%s\n" % combine)
                    return False
            except Exception as e:
                log.error(e)
                return False
    if combine == 1:
        for classes in sorted(allSet):
            fhOut.write('%s\n' % (classes))  #            fhOut.write(classes+','+fileNames+"\n")
    elif combine == 0:
        for numbers in sorted(allList):
            fhOut.write('%s\n' % (numbers))  #            fhOut.write(numbers+','+fileNames+'\n')
    fhOut.close()

# def columns2profileMatrix(rowCandiList, colCandiList, outPathAndFile,rowNameList=0, colNameList=0):
#     if len(rowCandiList)==len(colCandiList):
#         dataDict = dict(zip(rowCandiList,colCandiList))
#     else:
#         print('Input length of candiList of row and column won\'t match\n')
#         return False
#     fhOut = open(outPathAndFile,'w')
#     fhRow = open(outPathAndFile+'_rowName','w')
# #    fhCol = open(outPathAndFile+'_colName','w')
#     if rowNameList == 0:
#         rowNameList = rowCandiList
#     if colNameList == 0:
#         colNameList = colCandiList
#     for currentCol in colNameList:
#         fhOut.write("%s\t" % currentCol)
#     fhOut.write("\n")
#     allElement = [[0 for x in range(len(colNameList))] for y in range(len(rowNameList))]
#     for currentRow in rowNameList:
#         fhRow.write('%s\n' % currentRow)
#         for checkKeys in dataDict.keys():







def dirname2list(dir, out):

    diseaseList = os.listdir(dir)
    try:
        fhOut = open(out,'w')
        number = 0
        for i in range(1,len(diseaseList)):
            if os.path.isdir(dir + '/' + diseaseList[i]):
                if re.findall('^.', diseaseList[i]):
                    next
                number += 1
                fhOut.write('%s\t%s\n' % (number, diseaseList[i]))
    except Exception as e:
        log.error (e)
        print("ERROR in function dirname2list:%s" % e)
        return None

def sepCol2List(dirName, targetFilename, columnNum, outFilePathandName, combine=1, header=1):
    """
    This function gets all features in different diseases(separated dir, same file name)
    input1: dir name
    input2: target file name
    input3: number of column want to catch (start from 1)
    input4: out file name
    input5: (optional) whether to combine (0,1)
    output:data/preparation/allCDE
    """

    print("Cautious: read CSV files as default! (if not, check ',' separator)\n")
    if re.findall('/$',dirName):
        dir = dirName
    else:
        dir = dirName + '/'
    target = targetFilename#

    fileNames = os.listdir(dir)
    allClass = set()
    allNumber = []
    for i in fileNames:
        if os.path.isfile(dir.join(fileNames)):
            print("there are unexpected files in data folder: omitted\n")
            next
        skip = header #header
        newPath = os.path.join(dir, i)
        newTarget = os.path.join(newPath,target)
        try:
            if os.path.exists(newTarget):
                fhIn = read_csv(newTarget)
                for line in fhIn:
                    if skip == 1:
                        skip = 0
                        continue
                    col = columnNum - 1
                    if combine == 0:
                        allNumber.append(line[col])
                    elif combine == 1:
                        allClass.add(line[col])
                    else:
                        print("error: combine or not?")

        except Exception as e:
            log.error(e)
            print("ERROR in function mapDiseaseName:%s" % e)
            return False
    if combine == 1:
        print('Combine or not: %s\t%s number in total: %s ' % (combine,outFilePathandName, len(allClass)))
    elif combine == 0:
        print('Combine or not: %s\t%s number in total: %s ' % (combine,outFilePathandName,len(allNumber)))
#output file
    outFileName = outFilePathandName
    try:
        fhOut = open(outFileName, 'w')
    except:
        print "open outfile error in get_sep_col_2_list\n"
        return False
    if combine == 1:
        for classes in sorted(allClass):
            fhOut.write('%s\n' % (classes))  #            fhOut.write(classes+','+fileNames+"\n")
    elif combine == 0:
        for numbers in sorted(allNumber):
            fhOut.write('%s\n' % (numbers))  #            fhOut.write(numbers+','+fileNames+'\n')
    fhOut.close()
#sepCol2List(dataRootPath,targetFile,2,"../result/alLFreq", 0)

def integrateSepFiles(dirName, targetFileName, outFilePath):
    print("Cautious: read CSV files as default! (if not, check ',' separator)\n")
    id = 0
    head = ""
    firstFile = True
    fhOut = open(outFilePath,'w')
    if re.findall('/$',dirName):
        dir = dirName
    else:
        dir = dirName + '/'
    for current in os.listdir(dirName):
        pathAna = dir + current + '/'
        if os.path.isdir(pathAna):
            if os.path.exists(pathAna+targetFileName):
                id += 1
                try:
                    fhIn = readCSV(pathAna + targetFileName, "F")
                    if not firstFile:
                        test = fhIn.pop(0)
                        if not test == head:
                            print("header in files disagree:%s!=%s\nfiles:%s\n" % (test, head,pathAna+targetFileName))
                            return False
                    if firstFile:
                        head = fhIn.pop(0)
                        fhOut.write("id\tname")
                        for j in head:
                            fhOut.write("\t%s" % j)
                        fhOut.write('\n')
                        firstFile = 0
                    for line in fhIn:
                        fhOut.write("%s\t%s" % (id,current))
                        for i in line:
                            fhOut.write("\t%s" % (i))
                        fhOut.write('\n')

                except IOError:
                    print("error in integrateSepFiles")
                    return False
    fhOut.close()
#integrateSepFiles(dataRootPath,targetFile,"../result/"+group+"_table")

def sepLine2List(file2sep, sep="-"):
    """

    :param file:
    default:combine
    """
    try:
        fhIn = open(file2sep,'r')
        all = set()
        for line in fhIn:
            if line == '':
                next()
            line = line.replace(sep,'\n').rstrip()
            line = line.split('\n')
            for i in line:
                all.add(i)
        fhIn.close()

        fhOut = open(file2sep+'_sep','w')
        for i in all:
            fhOut.write(i+'\n')
        fhOut.close()
        print("Entry number left: %s" % len(all))
    except:
        print "error in sepLine2List"
        return False
#sepLine2List("../result/allST"," - ")
