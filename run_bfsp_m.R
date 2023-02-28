data<-read.csv("example_data.csv") #read in data
set.seed(1) #set random seed for reproducible results 
result<-bd_bfsp(data$x,data$y,data$mcoord,iterations=5000,burn=2000,max_lesion=4,spatial=F) #run bfsp-m algorithm
mask<-post_process_bfspm(result,prob = .5) #post-process 
write.csv(result,"bfsp_m_results.csv") #save results as .csv file 
