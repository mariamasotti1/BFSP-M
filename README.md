# BFSP-M 
## Bayesian Functional Spatial Partitioning for Multiple Lesion Detection 

Example is included as "example_data.csv". The coordinates are saved as the first 2 columns of the matrix. The last column contains the mcoord voxel-wise probabilities. 

The file "bfsp_m.R" contains functions, libraries to run the BFSP-M algorithm.  

The file "run_bfsp_m.R" generates lesion masks for example data. 

### bd_bfsp()

Estimates boundaries and gaussian spatial processes in each region via MCMC. 

#### Arguments

x: vector of length n giving voxel spatial coordinate x

y: vector of length n giving voxel spatial coordinate y

value: vector of voxel-wise values (mcoord)

iterations: total number of MCMC iterations 

burn: number of MCMC iterations to discard (must be less than iterations)

max_lesion: the maximum possible number of lesions contained in the image  

spatial: boolean, if T: models autocorrelation with exponential gaussian process, if F: assumes data follow independent normal distributions

#### Value

list containing MCMC samples of estimated parameters  

### post_process_bfspm()

#### Arguments

result: output from bd_bfsp

prob: posterior lesion probability cutoff (default is .5 meaning that the lesion was present in at least half of MCMC samples) 

#### Value

a dataframe with x,y,value, cluster, cluster_uncert, average_voxel_value. 
cluster: a mask indicating which cluster each voxel belongs to (0 is the healthy region)
cluster_uncert: the posterior lesion probabiliy specific to each lesion
average_voxel_value: the average value of mcoord in each cluster 
