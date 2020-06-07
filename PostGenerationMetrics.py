import operator
def Reverse(lst): 
    return [ele for ele in reversed(lst)] 
def file_len(fname):
    with open(fname) as g:
        for i, l in enumerate(g):
            pass
    return i + 1
def pick_top_n_diversity(n):#fully integrated function to generate the top n number of sequences by minimizing normality score within the domain and ranking them. Returns sorted dict least to greatest normality score
	v = open("Metrics.txt", "r")
	start_position_1=11
	start_position_2=2
	listofseqs=[]
	listofnormvals=[]
	triplelist=[]
	a=0
	for x in range(file_len("Metrics.txt")):

		temp=v.readline()
		if(temp[0]=="M"):
			triplelist.append(temp[7:8])

		if(start_position_1%12==0): #the sequence
			triplelist.append(temp[0:-1]) 

		if(start_position_2%12==0):#the normality
			triplelist.append(temp[11:-1]) 

		start_position_1+=1
		start_position_2+=1
	i=0
	approved=0
	for item in triplelist:
		if(i%3==0):
			if(triplelist[i+1]=="V"):
				approved+=1
				listofseqs.append(triplelist[i])
				listofnormvals.append(triplelist[i+2])
		i+=1
	dictionary = dict(zip(listofseqs, listofnormvals))
	sorted_x = sorted(dictionary.items(), key=operator.itemgetter(1))
	print(str(approved)+" of the "+str(len(triplelist)/3) +" sequences were approved by refolding")
	v.close()
	return sorted_x[0:n]

def pick_top_n_editdist(n):#fully integrated function to generate the top n number of sequences by minimizing normality score within the domain and ranking them. Returns sorted dict least to greatest normality score
	v = open("Metrics.txt", "r")
	start_position_1=11
	start_position_2=3
	listofseqs=[]
	listofnormvals=[]
	triplelist=[]
	for x in range(file_len("Metrics.txt")):
		temp=v.readline()
		if(temp[0]=="M"):
			triplelist.append(temp[7:8])
		if(start_position_1%12==0): #the sequence
			triplelist.append(temp[0:-1]) 
		if(start_position_2%12==0):#the normality
			triplelist.append(float(temp[33:-1])) 
		start_position_1+=1
		start_position_2+=1
	i=0
	approved=0
	for item in triplelist:
		if(i%3==0):
			if(triplelist[i+1]=="V"):
				approved+=1
				listofseqs.append(triplelist[i])
				listofnormvals.append(triplelist[i+2])
		i+=1
	dictionary = dict(zip(listofseqs, listofnormvals))
	sorted_x = sorted(dictionary.items(), key=operator.itemgetter(1))
	v.close()
	reversedlist= Reverse(sorted_x)
	return reversedlist[0:n]

def comparediversityratings(diversity,editdist): #compares the Levenshteinian diversity ranking to the SVM generated diversity rankings as a percent shared 
	iterant=0
	a=0
	for diversitycandidateindex in range(len(diversity)):
		for editdistcandidateindex in range(len(editdist)):
			if((diversity[diversitycandidateindex])[0]==(editdist[editdistcandidateindex])[0]):
				a+=1
	print(str(a/float(str(len(diversity))+'.0')*100)+ "% "+"of sequences in the SVM diversity ranking within the range were shared by the Levenshtein diversity ranking")
		
comparediversityratings(pick_top_n_diversity(5),pick_top_n_editdist(5))







