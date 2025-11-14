import coolpuppy.lib.io
import glob

import os
import sys

mypath = "./"
inFiles = glob.glob("hic*.clpy")
print(inFiles)

window = int(sys.argv[1]) # Get the windowxwindow squared map
print(window)

for inFile in inFiles:

    print(inFile)
    inFileName = os.path.basename(inFile)
    inFileName = os.path.splitext(inFileName)[0]
    print(inFileName)

    pileup_df = coolpuppy.lib.io.load_pileup_df(inFile)

    # I use the 'resolution' field to get the resolution of the map
    resolution = pileup_df['resolution'][0]
    print("Coolpuppy pile-up matrix at",resolution,"bp resolution")

    # I use the 'data' field, which contains a NxN, where N is the size of the resulting pileup squared matrix
    # In this case the matrix is 21x21, where 21 comes from a region of (100kb (left-flanking) +
    # 10kb (region of interest) + 100kb (right-flanking)) at 10kb resolution
    # n is the counter of the Coolpuppy pile-up matrices contained in inFile
    n = 1
    for matrix in pileup_df['data']:
        #print(type(matrix))
        
        size = len(matrix)
        print("Size of the Coolpuppy pile-up matrix",size)
        print("Matrix indexes go from 0 to",size-1)
        
        print("Central part")
        outFile   = inFileName + "_central_%d.txt" % window 
        fpoutFile = open(outFile,"w") 
        iStart = int((size-1)/2)-int(window/2)
        iEnd   = int(iStart + window)-1
        print("The",window,"x",window,"sub-matrix of interest goes from index",iStart,"to index",iEnd,"included")

        for i in range(iStart,iEnd+1):
            for j in range(iStart,iEnd+1):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))
        fpoutFile.close()

                
        print("Top-left part")
        outFile = inFileName + "_topLeft_%d.txt" % window
        fpoutFile = open(outFile,"w")         
        yiStart = int(size-window)                
        yiEnd   = int(size-1)        
        xiStart = 0
        xiEnd   = int(xiStart + window)-1
        print("The",window,"x",window,"sub-matrix of interest goes from index",xiStart,"to index",xiEnd,"included in x")
        print("The",window,"x",window,"sub-matrix of interest goes from index",yiStart,"to index",yiEnd,"included in y")        

        for i in range(xiStart,xiEnd+1):
            for j in range(yiStart,yiEnd+1):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))
        fpoutFile.close()

              
        print("Bottom-left part")
        outFile = inFileName + "_bottomLeft_%d.txt" % window
        fpoutFile = open(outFile,"w")                 
        yiStart = 0
        yiEnd   = int(yiStart + window)-1
        xiStart = 0
        xiEnd   = int(xiStart + window)-1
        print("The",window,"x",window,"sub-matrix of interest goes from index",xiStart,"to index",xiEnd,"included in x")
        print("The",window,"x",window,"sub-matrix of interest goes from index",yiStart,"to index",yiEnd,"included in y")        

        for i in range(xiStart,xiEnd+1):
            for j in range(yiStart,yiEnd+1):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))                
                #print(inFile,n,i,j,matrix[i][j])
        fpoutFile.close()
        
        print("Top-right part")
        outFile = inFileName + "_topRight_%d.txt" % window
        fpoutFile = open(outFile,"w")                 
        yiStart = int(size-window)                
        yiEnd   = int(size-1)        
        xiStart = int(size-window)
        xiEnd   = int(size-1)
        print("The",window,"x",window,"sub-matrix of interest goes from index",xiStart,"to index",xiEnd,"included in x")
        print("The",window,"x",window,"sub-matrix of interest goes from index",yiStart,"to index",yiEnd,"included in y")        

        for i in range(xiStart,xiEnd+1):
            for j in range(yiStart,yiEnd+1):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))                
        fpoutFile.close()
                
        print("Bottom-right part")
        outFile = inFileName + "_bottomRight_%d.txt" % window
        fpoutFile = open(outFile,"w")                 
        yiStart = 0
        yiEnd   = int(yiStart + window)-1
        xiStart = int(size-window)
        xiEnd   = int(size-1)
        print("The",window,"x",window,"sub-matrix of interest goes from index",xiStart,"to index",xiEnd,"included in x")
        print("The",window,"x",window,"sub-matrix of interest goes from index",yiStart,"to index",yiEnd,"included in y")        

        for i in range(xiStart,xiEnd+1):
            for j in range(yiStart,yiEnd+1):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))                
        fpoutFile.close()
