library(gkmSVM,readtext)
require(dplyr)

posfn= 'positives.fa' 
negfn= 'negatives.fa' 
kernelfn= 'test_kernel.txt'
setwd("/Users/Peter/PycharmProjects/RNALibrary/")
fileName <- 'InverseSequenceOutput.txt'
a=readChar(fileName, file.info(fileName)$size)
b=strsplit(a,"\n")
listofseqs<-c()
k<-1
for (i in 4:length(b[[1]])){ #sample the sequences in the format predefined by the testrun.txt seq output
  if((i-4)%%5==0){
    listofseqs[[k]]<-((b[[1]])[[i]])
    k<-k+1
  }
}

negative_sequence_names<-c() #Negative Seq FASTA labels
k<-1
for (i in 2:length(listofseqs)-1){
  negative_sequence_names[[k]]=(paste("negativesequence",i,sep=""))
  k<-k+1
}

positive_sequence_name="positive sequence" #Positive Seq FASTA label
positive_file_name="positiveseq.fa"
negative_file_name="negativeseq.fa"

writeFasta<-function(data, filename){ #c/o https://bootstrappers.umassmed.edu/guides/main/r_writeFasta.html referenced by Peter Bowman-Davis on June 4, 2020
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}
diversity_score<-list() #Ordered list of sequence diversity score computed by the average of the gkm-svm kernel scores (first column of the triangular matrix slice as described https://cran.r-project.org/web/packages/gkmSVM/gkmSVM.pdf )
for (i in 1:length(listofseqs)){ #triple switch: first seq, last seq, middle n-2 seqs because I am bad at R
  if(i==1){
    candidateseq=listofseqs[i]
    negativeseqs=listofseqs[2:length(listofseqs)]
    negseqdata = dplyr::data_frame(name = negative_sequence_names,seq = negativeseqs)
    writeFasta(negseqdata, negative_file_name)
    posseqdata = dplyr::data_frame(name = positive_sequence_name,seq = candidateseq)
    writeFasta(posseqdata, positive_file_name)
    gkmsvm_kernel(positive_file_name, negative_file_name, paste("test_kernel",".txt",sep=""));
    #Generate Kernel Diversity Score
    a=readChar("test_kernel.txt", file.info(fileName)$size)
    b=strsplit(a,"\t")
    total<-0
    for (c in 1:length((b[[1]]))){
      if(substr(toString(((b[[1]])[c])), 1,1)=="\n"){
        k=substr(toString(((b[[1]])[c])),2,nchar(toString(((b[[1]])[c]))))
        if(is.na(as.double(k))==FALSE){
          total=total+as.double(k)
        }
      }
    }
    diversity_score[i]=total/length((b[[1]]))
    #print(diversity_score)
    
  }
  else if(i==length(listofseqs)){
    candidateseq=listofseqs[i]
    negativeseqs=listofseqs[1:i-1]
    negseqdata = dplyr::data_frame(name = negative_sequence_names,seq = negativeseqs)
    writeFasta(negseqdata, negative_file_name)
    posseqdata = dplyr::data_frame(name = positive_sequence_name,seq = candidateseq)
    writeFasta(posseqdata, positive_file_name)
    gkmsvm_kernel(positive_file_name, negative_file_name, paste("test_kernel",".txt",sep=""));
    #Generate Kernel Diversity Score
    a=readChar("test_kernel.txt", file.info(fileName)$size)
    b=strsplit(a,"\t")
    total<-0
    for (c in 1:length((b[[1]]))){
      if(substr(toString(((b[[1]])[c])), 1,1)=="\n"){
        k=substr(toString(((b[[1]])[c])),2,nchar(toString(((b[[1]])[c]))))
        if(is.na(as.double(k))==FALSE){
          total=total+as.double(k)
        }
      }
    }
    diversity_score[i]=total/length((b[[1]]))
    #print(diversity_score)
    
  }
  else {
    candidateseq=listofseqs[i]
    negativeseq1=listofseqs[1:i-1]
    negativeseq2=listofseqs[(i+1):length(listofseqs)-i]
    negativeseqs=c(negativeseq1,negativeseq2)
    print(length(negativeseqs))
    negseqdata = dplyr::data_frame(name = negative_sequence_names,seq = negativeseqs)
    writeFasta(negseqdata, negative_file_name)
    posseqdata = dplyr::data_frame(name = positive_sequence_name,seq = candidateseq)
    writeFasta(posseqdata, positive_file_name)
    gkmsvm_kernel(positive_file_name, negative_file_name, paste("test_kernel",".txt",sep=""));
    #Generate Kernel Diversity Score
    a=readChar("test_kernel.txt", file.info(fileName)$size)
    b=strsplit(a,"\t")
    total<-0
    for (c in 1:length((b[[1]]))){
      if(substr(toString(((b[[1]])[c])), 1,1)=="\n"){
        k=substr(toString(((b[[1]])[c])),2,nchar(toString(((b[[1]])[c]))))
        if(is.na(as.double(k))==FALSE){
          total=total+as.double(k)
        }
      }
    }
    diversity_score[i]=total/length((b[[1]]))
  }
}
#file cleanup to prevent confusion/ambiguity in filesystem
if (file.exists("positiveseq.fa")) 
  file.remove("positiveseq.fa")
if (file.exists("negativeseq.fa")) 
  file.remove("negativeseq.fa")
if (file.exists("test_kernel.txt")) 
  file.remove("test_kernel.txt")
#dump the diversity score list to a text file for later parsing in python
capture.output(diversity_score, file = "gkm_svm_output.txt")

