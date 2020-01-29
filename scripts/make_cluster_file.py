import sys

#function to make cluster file
def makeOrthClusters(proteinClusters):    
    clusterCount = 0
    clusters = open(proteinClusters, 'r')
    outfile = open("ProteinOrthoClusters_reformatted.tsv","w+")
    for line in clusters:
        if(line.startswith("#")):
            header = line.split()[4:]
        else:
            clusterLine = line.split("\t")[3:]
            clusterLine[-1] = clusterLine[-1].strip()
            combinedClusterRow = list()
            for lineLoc in range(0,len(clusterLine)):
                clusterRow = header[lineLoc] +":"+ clusterLine[lineLoc]
                if "*" not in clusterRow:
                    combinedClusterRow.append(clusterRow)
            outfile.write(str(clusterCount)+ "\t" + "\t".join(combinedClusterRow) + "\n")
            clusterCount += 1
    clusters.close()
    outfile.close()

makeOrthClusters(snakemake.input[0])