#get the file name without path information
getFileFullNameNopath=function(fullfnames){
   res = NULL
   for(each in fullfnames) {
     splitted=unlist( strsplit(each, "/") )
     res= c(res, splitted[length(splitted) ])
   }
   return (res)
}