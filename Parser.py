import operator
def file_len(fname):
    with open(fname) as g:
        for i, l in enumerate(g):
            pass
    return i + 1
z=[]
b=[]
def pick_top_n_diversity(n):
	v = open("SequenceMetrics.txt", "r")
	c=0
	o=0
	it=11
	it2=2
	t=0
	listofseqs=[]
	listofnormvals=[]
	for x in range(file_len("Metrics.txt")):
		temp=v.readline()
		if(temp[0]=="M"):
			if(temp[7:8]=="V"):
				c+=1
			o+=1
		if(it%12==0):
			listofseqs.append(temp[0:-1])
		if(it2%12==0):
			listofnormvals.append(float(temp[11:-1]))
		it+=1
		it2+=1
	dictionary = dict(zip(listofseqs, listofnormvals))
	sorted_x = sorted(dictionary.items(), key=operator.itemgetter(1))
	v.close()
	return sorted_x[0:n]

for x in pick_top_n_diversity(20):
	print(x)
# print(t)
# v.close()
# def topseq(number_of_sequences):

#print(str(float(c*100/o))+"% verification success yielding "+str(c)+" qualified sequences")








