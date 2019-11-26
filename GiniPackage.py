#!/usr/bin/env python3
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import sys
import csv


#Create GiniObject
#Optional parameters: 
###### genes-cell sample file (required for further manipulations). m x n matrix: Genes as rows and cell samples as columns
###### t-SNE coordinates matrix (Assumes: row names, x coordinates, y coordinates corresponding to 1, 2nd and 3rd columns)
###### clusters file (1 column file, n rows corresponding to n cell samples)
class GiniObject:
    def __init__(self, fname=None, tsnefname=None, clusterfname=None):
        self.counts = []
        self._TSNEcoords = np.array([])
        self.clusters = np.array([])
        self.thresval = []
        self.allgenes = []
        self.allgini = []
        self.thresval = 0
        self.numCluster = 0
        if len(sys.argv) != 1 and len(sys.argv) != 3:
            print("Please provide either just the counts table or counts table, t-SNE Coordinates and the cluster information for each sample")
            sys.exit()
        # User provides name/path of the csv file that contains gene-to-sample data
        if fname:
            self.compile(fname)
        if tsnefname:
            self.setTSNE(tsnefname)
        # User provides name/path of csv file with cluster information for each sample
        if clusterfname:
            self.setClusters(clusterfname)

    #Set counts matrix; Only takes in numeric table with no gene or header information. If you have csv file with extra information, please use Initialization method
    #Alternatively user can use self.compile method for csv file with gene and header information
    #parameters: numpy m x n matrix; m genes and n cell samples
    def setCounts(self, counts):
        self.counts = counts

    #Set t-SNE coordinates (fname: csv file with rownames, x and y coordinates as 1, 2nd and 3rd columns of the file)
    #parameters: csv file name (with path if not in the same directory)
    def setTSNE(self, fname):
        self._TSNEcoords = np.loadtxt(fname, delimiter=",", skiprows= 1, usecols=(1,2))
    
    #Set t-SNE coordinates (fname: csv file with x and y coordinates as 1 and 2nd columns of the file)
    #parameters: csv file name (with path if not in the same directory)
    def setTSNE2(self, fname):
        self._TSNEcoords = np.loadtxt(fname, delimiter=",", skiprows=1)

    #Set list of genes manually (genes: array of genes)
    #parameters: python list of genes
    def setGenes(self, genes):
        self.allgenes =  genes
    
    #Set clusters (fname: csv file with cluster information; 1 column file with n rows, corresponding to n cell samples)
    #parameters: csv file name (with path if not in the same directory)
    def setClusters(self, fname):
        self.clusters = np.loadtxt(fname, delimiter=",")

    #Set counts alternate method (csv file with gene and header information). Method automatically fetches gene information (assumed to be 1st column)
    #parameters: csv file name (with path if not in the same directory)
    def compile(self, fname):
        f = open(fname, 'r')
        f.readline()
        self.allgenes = [row.split()[0].split(',')[0] for row in f]
        f.close()
        with open(fname) as f:
            ncols = len(f.readline().split(','))
        f.close()
        self.counts = np.loadtxt(fname, dtype ='uint32', delimiter=",", skiprows=1, usecols=range(1,ncols))


    #Call runTSNE if you do not have t-SNE coordinates.
    #Optional parameters: dims (for # of dimensions to be returned from PCA).
    def runTSNE(self, dims=None):
        if dims is None:
            dims = 0
        xpca = PCA(n_components=dims).fit_transform(np.transpose(self.counts))
        self._TSNEcoords = TSNE(learning_rate=500).fit_transform(xpca)

    #Checks if all parameters are available for calculating gini coeff and plotting 2D Gini
    def type_check(self, genelist):
        if (self._TSNEcoords.size == 0):
            print(
                "Please either provide t-SNE coordinates (setTSNE) or run t-SNE (runTSNE)")
            sys.exit()
        if (not genelist):
            print(
                "Please provide a valid gene(s) list; Ex: [\"ACTC1\", \"MYH7\"]")
            sys.exit()
        if (self.clusters.size == 0):
            print("Please provide clusterID information; You can use setClusterID method to do so.")
            sys.exit()

    #DrawGini method calls make2DGini to draw 2D Gini plot
    #parameters: Optional parameter to draw and save the 2D Gini Index plots
    def DrawGini(self, draw=False):
        genelist =  self.allgenes
        self.type_check(genelist)
        clusters = len(np.unique(self.clusters))
        self.numCluster = clusters
        self.allgini = np.zeros((len(genelist), clusters))
        for i in range(0, len(genelist)):
            for j in range(0, clusters):
                coords = self._TSNEcoords[self.clusters == j,:] # TSNE coordinates belonging to cluster k
                geneind = self.allgenes.index(genelist[i]) # index of the gene
                clusID = np.zeros((len(coords),)) #Initialize GINI clusters (0 representing all and 1 representing expressed)
                for k in range(0, len(coords)):
                    if (self.counts[geneind, k] > self.thresval):
                        clusID[k] = 1
                if sum(clusID) > 0:
                    self.allgini[i, j] = self.make2DGini(coords, clusID, genelist[i], j, draw)
                else:
                    self.allgini[i, j] = -1
        self.SaveGini()

    #Save all the Gini values in a CSV file in the directory of the GiniPackage python script
    def SaveGini(self):
        df = pd.DataFrame(self.allgini, index=self.allgenes)
        a = [i for i in range(0, self.numCluster)]
        df.to_csv(r'allgini.csv', sep='\t', encoding='utf-8', header=a)

    #Method to set clusterIDs (2D Gini) and to draw 2D Gini curve. Calls computeGini for calculating gini coeff
    #Not intended for direct usage. Please use DrawGini method
    def make2DGini(self, coords, clusID, gene, cluster, draw):
        anglelist = np.arange(0, 360, 360/1000)
        numCluster = len(np.unique(clusID))
        # x is your matrix (in this case, tSNE coordinates), numCluster is a 1D matrix/vector indicating which clusterID each coordinate belongs to
        gini = np.zeros((anglelist.shape[0],2))
        for i in range(0, numCluster):
            if i == 1:
                coordinate = coords[clusID == i,:]
            else:
                coordinate = coords
            for j in range(0, anglelist.shape[0]):
                angle = anglelist[j]
                anglearr = [math.cos(angle), math.sin(angle)]
                value = np.dot(coordinate, anglearr)
                valarr = np.subtract(value, min(value))
                valarr = np.asarray(valarr)
                gini[j, i] = self.computeGini(np.ones((valarr.shape[0],1)), valarr)
            plt.polar(anglelist, gini, linewidth=0.25)
        if draw:
            plt.savefig(gene+"_Cluster"+str(cluster)+".png", bbox_inches='tight')
            plt.close()
        ginival = math.sqrt(sum((gini[:,0] - gini[:,1])**2)/1001)
        return ginival


    #Calculates gini coefficient
    #Not intended for direct usage. Please use DrawGini method
    def computeGini(self, pop, val):
        pop = [0] + pop
        val = [0] + val
        #assert (pop>=0).all(), 'Pop matrix: Some values are less than zero; check if you passed in the right values'
        #assert (val>=0).all(), 'Value matrix: Some values are less than zero; check if you passed in the right values'
        #z = pop*val
        z = val
        ord = np.argsort(val)
        pop = pop[ord]
        z = z[ord]
        pop = np.cumsum(pop)
        z = np.cumsum(z)
        relpop = pop/pop[-1]
        relz = z/z[-1]
        g = 1-np.sum((relz[0:-1]+relz[1:]) * np.diff(relpop))
        return g

