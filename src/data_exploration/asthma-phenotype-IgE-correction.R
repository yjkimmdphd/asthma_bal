################################################################
# revise the asthma phenotype data frame
# replace values that have "<" or ">"with an actual numeric value
# replace nonsensical values
################################################################

library(tidyverse)
# asthma biomarker phenotype file, nasal, saved in  'phenotype'
filename2<-("C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/revised_data/asthma-phenotype-filtered.csv")
file.exists(filename2)
phenotype<-read.csv(filename2)

r.neg.titers<-runif(2243,min=0.001, max=0.09)%>%round(2)
phenotype.dup<-phenotype # duplicated phenotype data
c<-grep("<",phenotype.dup) # which col has a data containing "<"
d<-grep(">",phenotype.dup) # which col has a data containing ">"
phenotype.dup[,c]
phenotype.dup[,d]

# fx that finds whether the data contains "<"
findLess<-function(data){
  grepl("<",data)
}
# fx that finds whether the data contains ">"
findMore<-function(data){
  grepl(">",data)
}

# ask if each cell in each column of phenotype.dup has "<", and save TRUE/False in each cell. Save the resulting data.frame as neg.titers 
# also make a list of row numbers after finding which row of each column has "<"
neg.titers<-sapply(phenotype.dup[,c],findLess)
neg.titers.col<-colnames(neg.titers)
neg.titers.row<-apply(neg.titers,2,function(r){which(r==TRUE)})

# same as above, but with "<"
high.titers<-sapply(phenotype.dup[,d],findMore)
high.titers.col<-colnames(high.titers)
high.titers.row<-apply(high.titers,2,function(r){which(r==TRUE)})

# some of the titers have weird entry. define a function 'findVal' that looks for each specific weird entries

findVal<-function(val){
  e<-apply(phenotype.dup,2,function(r){r==val})
  f<-colSums(e,na.rm=TRUE)
  which(f>0)
}

# which one has ">0.35". Then replace that cell with a random value <0.09
e.c<-findVal(">0.35")%>%as.numeric # which column has ">0.35"
e.r<-which(phenotype.dup[,e.c]==">0.35")
print(e.c)
print(e.r)
phenotype.dup[e.r,e.c]<-runif(n=1,min=0.01, max=0.09)%>%round(2)

# which cells have ">5000", then replace each of these cells with 5000
f.c<-findVal(">5000")%>%as.numeric # which column has ">5000"
f.r<-sapply(phenotype.dup[,f.c],function(r){which(r==">5000")})
print(f.c)
print(f.r)
for(i in 1:length(f.c)){
  l<-f.c[i]
  k<-f.r[[i]]
  phenotype.dup[k,l]<-5000
}


# which cells have ">27.5", then replace each of these cells with 27.5
g.c<-findVal(">27.5")%>%as.numeric
g.r<-which(phenotype.dup[,g.c]==">27.5")
phenotype.dup[g.r,g.c]<-27.5

for(i in 1:length(neg.titers.col)){
  k<-neg.titers.col[[i]]
  l<-neg.titers.row[[k]]
  phenotype.dup[l,k]%>%print
  phenotype.dup[l,k]<-sample(r.neg.titers, 1)
  phenotype.dup[l,k]%>%print
}

phenotype.dup[,neg.titers.col]

# check again which ones have ">" in the values. The remaining such entries should be ">100". Replace with 101
high.titers<-sapply(phenotype.dup[,d],findMore)
high.titers.col<-colnames(high.titers)
high.titers.row<-apply(high.titers,2,function(r){which(r==TRUE)})

for(i in 1:length(high.titers.col)){
  k<-high.titers.col[[i]]
  l<-high.titers.row[[k]]
  phenotype.dup[l,k]%>%print
  phenotype.dup[l,k]<-101
  phenotype.dup[l,k]%>%print
}

# some IgE testing columns test "Pos" or "Neg" allergic to certain allergens
# replace any nonsensical NA entry with a true 'NA'
# factorize 'Pos' and 'Neg'
col.initial.titer.pn<-c(104:113)
phenotype.dup[,col.initial.titer.pn]
h.c<-findVal("N/A")%>%as.numeric
h.r<-sapply(phenotype.dup[,h.c],function(r){which(r=="N/A")})
for(i in 1:length(h.c)){
  l<-h.c[i]
  k<-h.r[[i]]
  phenotype.dup[k,l]<-NA
}
kk<-sapply(phenotype.dup[,col.initial.titer.pn],function(r){factor(r,levels=c("Pos","Neg"),exclude=NULL)})
phenotype.dup[,col.initial.titer.pn]<-kk

# cols with actual IgE titers are currently character vectors. change to numeric
col.initial.titer<-c(114:118)
col.repeat.titer<-c(131:144)
phenotype.dup[,col.initial.titer]<-sapply(phenotype.dup[,col.initial.titer],as.numeric)
phenotype.dup[,col.repeat.titer]<-sapply(phenotype.dup[,col.repeat.titer],as.numeric)

write.csv(phenotype.dup,"C:/Users/kimyo/Dropbox/Research/asthma-allergy-bunyavanich/data/revised_data/asthma-phenotype-filtered-revised.csv")


