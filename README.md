# 2D Gini Index

Calculates Gini RMSD values for each gene/cluster.
Special thanks to Dr. Jake Chen, Dr. Thanh Nguyen, Dr. Xu Nuo and Zongliang Yue.

### Required packages:

    numpy, pandas, mathplotlib, sklearn.
    Install these packages before using pip3 command (ex: sudo pip3 install numpy) or run setup.sh (on linux) to autoinstall them.

## Class Gini Object:
    
    Allows user to create a new instance of Gini Object. Each Gini Object
    corresponds to a new set of data.
    
    Subsequent methods allow user to calculate Gini RMSD values for every
    gene in each cluster. If a gene is non-expressed in one cluster, then
    value of -1 is returned.
    
### Class Parameters:

    fname:
        Specify a genes-sample matrix file. Assumes the first row corresponds 
        to cell-sample information and the first column corresponds to list of
        genes. If not, refer to setCounts and setGenes method to set genes-sample
        matrix and genes list.
        (path to file)
    
    tsnefname:
        Specify t-SNE coordinates corresponding to genes-sample matrix file
        (path to file)
    
    clusterfname:
        Specify a cluster information file corresponding to each cell sample
        (path to file)
        
    At minimum, genes-sample matrix file is required. If users specify a t-SNE
    Coordinate file, they must also specify a cluster information file. Based on
    a clustering algorithm. Ex: k-means).
    

### Methods:

    setCounts(self, counts)
        Desc: Set counts matrix; Only takes in numeric table with no gene or header information. If you have csv file with extra information, please use Initialization method
            Alternatively user can use self.compile method for csv file with gene and header information
        Parameters: numpy m x n matrix; m genes and n cell samples
    
    setTSNE(self, fname)
        Desc: Set t-SNE coordinates (fname: csv file with rownames, x and y coordinates as 1, 2nd and 3rd columns of the file)
        Parameters: csv file name (with path if not in the same directory)
    
    setTSNE2(self, fname)
        Desc: Alternative method; Set coordinates (fname: csv file with x and y coordinates as 1 and 2nd columns of the file)
        Parameters: csv file name (with path if not in the same directory)
    
    setGenes(self, genes)
        Desc: Set list of genes manually (genes: array of genes)
        Parameters: list of genes
        
    setClusters(self, fname)
        Desc: Set clusters (fname: csv file with cluster information; 1 column file with n rows, corresponding to n cell samples)
        Parameters: csv file name (with path if not in the same directory)
    
    compile(self, fname)    
        Desc: Set counts alternate method (csv file with gene and header information). Method automatically fetches gene information (assumed to be 1st column)
        Parameters: csv file name (with path if not in the same directory)
    
    runTSNE(self, dims=None)    
        Desc: runTSNE calculates tSNE coordinates. Run it if you do not have t-SNE coordinates.
        Parameters (optional): dims (for # of dimensions to be returned from PCA). Default is 0.
    
    type_check(self, genelist)
        Desc: Checks if all parameters are available for calculating gini coeff and plotting 2D Gini. Not intended for direct usage.
        Paramers: List of genes and a list of threshold values.
    
    DrawGini(self)
        Desc: DrawGini method calls make2DGini to draw 2D Gini plot.
        
    SaveGini(self)
        Desc: Saves Gini RMSD Values for each gene (row) and cluster (column) in a CSV file.
    
    make2DGini(self, coords, clusID, gene, cluster)
        Desc: Method to set clusterIDs (2D Gini) and to draw 2D Gini curve. Calls computeGini for calculating gini coeff
        Not intended for direct usage. Please use DrawGini method
        Parameters: List of coordinates, cluster IDs, gene and cluster
        
        
    computeGini(self, pop, val)
        Desc: Calculates gini coefficient
        Not intended for direct usage. Please use DrawGini method

###    Default Usage:
    
    Instantiante a new Gini Object with Genes-Cell Sample Matrix, tSNE Coordinates and Cluster IDs.
    If genes-cell sample matrix does not contain gene/cell sample information, use setCounts, setGenes,
    setTSNE and setCluster to set all of them manually.
    
    Call DrawGini method.
    
    Refer to tester.py for sample usage. Sample data files are available under Data folder.

