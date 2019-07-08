#title: roadmap_filt_dtr
#Author: Sunho Lee (sunho2519@gmail.com)
#date created: July.1.2019 (Albert Einstein College of Medicine)
#purpose: merge roadmap downloaded data to create a single table file
  #RAM limitation problem improvement from Seungsoo Kim's roadmap_filt.r;
    #resource expected: Free RAM of at least 3,800 MB  (3.8min)


  #################
  ####detour fx####
  #################
#function: reads rds file of a given path to filter and return its result
#return: a list, with its 1st index having the filtered table, and 2nd index having either 1 or 0 to indicate the exitence of the file
roadmap_filt_dtr = function(path) {

  rtn_li = list()     #returning object that contains the result and an indicator of file existence
  temp_enh = data.frame()   #oject to save the result

  road_df = try(readRDS(path))
  if('try-error' %in% class(road_df)) {     # if file not found
    cat(paste0('\n', path, " - file not found. Moving on to the next file"))
    rtn_li[[2]] = 1
  }
  else {        #if file found
    temp_enh = subset( road_df, name %in% c("13_EnhA1","14_EnhA2","15_EnhAF","16_EnhW1","17_EnhW2","18_EnhAc") )
    rtn_li[[2]] = 0
  }
  rtn_li[[1]] = temp_enh    #filtered or empty file saved
  return (rtn_li)
}


    ###################
    ######main fx######
    ###################
library(plyr)
library(data.table)
library(numbers)

args = commandArgs(trailingOnly=T)
hmsg = 'Rscript roadmap_filt_lp.r does not require any argument.'
if(length(args)>0) stop(hmsg)

source('src/pdtime.r')
t0 = Sys.time()

err_cnt = 0         #counts the number of file absent in the folder
no_file = 129       #total number of files
road_li = list()    #list object to save the read tables
dir = 'db/roadmap/'   #directory where downloaded roadmap data are saved
f_name = 'db/roadmap_enh.bed'    #name of the file to be written
cid = as.character(formatC(c(1:no_file),width=3,flag='0'))    #'001'~'129'

  #read road files and save them in a list
cat("\n>>> start reading files <<<\n")
for (i in 1:length(cid)) {
  if( mod(i, 10) == 0 ) cat( paste0('\n', i, '/', no_file, ' being processed') )
  path = paste0(dir,'E', cid[i],'_25_imputed12marks_hg38lift_dense.bed.rds')
  result = roadmap_filt_dtr(path)          #roadmap_filt_dtr returns the result as list. the 1st index contains the filtered table, 2nd contains 1 if error in reading file else 0 if no-error
  road_li[[i]] = result[[1]]                #filtered file is saved in the list
  if (result[[2]] == 1) {                   #if the file is not found, increase err_cnt by 1
    err_cnt = err_cnt + 1
  }
}
cat(paste0("\n\n>>> finished reading ", no_file, 'files <<<'))
road_enh = ldply(road_li, data.frame)      #the tables saved in the list is merged to form a single dataframe,
rm(road_li)

cat(paste0("\n\n>>> filtering in progress <<<"))
road_enh = road_enh[, c(1:3,6)]   #only col 1,2,3,6 are required
write.table(road_enh,f_name,row.names=F,col.names=F,quote=F,sep='\t')


cat(paste0('\n\n', err_cnt, ' files were not found \n'))
cat(paste0('File written: ',f_name,'\n'))
cat(paste0('>> ',pdtime(t0,1),'\n'))
