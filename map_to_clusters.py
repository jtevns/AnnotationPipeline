def MapToClusters(clusters,annotation):
    with open(clusters) as clusters:
        clusterDict = {x.split("\t")[0]:x.split("\t")[1:] for x in clusters }

    outFileName = "MappedToClusters_"+annotation[:-3]+"csv"

    with open(annotation) as annotations, open(outFileName,"w+") as outFile:
        for line in annotations:
            if "#" not in line:
                clusterNum = line.split(" ")[0].split("_cluster_")[1]
                foundClust = clusterDict[clusterNum]
                for geneList in foundClust:
                    if len(geneList.split(",")) > 1:
                        geneList = geneList.split(",")
                        prefix = geneList[0].split(":")[0]
                        genes =  [":".join([prefix,x]) for x in geneList[1:]]
                        for gene in genes:
                            outLine = ",".join([gene.strip()]+line.split()[1:])
                            outFile.write(outLine + "\n")
                    else:
                        outLine = ",".join([geneList.strip()]+line.split()[1:])
                        outFile.write(outLine+"\n")

MapToClusters(sys.argv[1],sys.argv[2])