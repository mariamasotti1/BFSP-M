source("bfsp_m.R")
data<-read.csv("example_data.csv") #read in data
set.seed(1) #set random seed for reproducible results 
result<-bfspm(data$x,data$y,data$mcoord,iterations=1000,burn=100,max_lesion=4,spatial=T) #run bfsp-m algorithm
mask<-post_process_bfspm(result,prob = .5) #post-process 
write.csv(mask,"bfsp_m_results.csv") #save results as .csv file 
