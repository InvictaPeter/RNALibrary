import operator
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
structlist=[]
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
	global structlist
	v = open("Metrics.txt", "r")
	start_position_1=11
	start_position_2=3
	listofseqs=[]
	listofnormvals=[]
	triplelist=[]

	
	for x in range(file_len("Metrics.txt")):
		temp=v.readline()
		if(temp[0]=="." or temp[0]=="("):
			structlist.append(temp[0:-1])
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
	#print(triplelist)
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
structlib=[]
def get_data():#fully integrated function to generate the top n number of sequences by minimizing normality score within the domain and ranking them. Returns sorted dict least to greatest normality score
	v = open("Metrics.txt", "r")
	start_position_1=11
	start_position_2=2
	start_position_3=3
	listofseqs=[]
	listofnormvals=[]
	triplelist=[]
	structlib=[]
	a=0
	for x in range(file_len("Metrics.txt")):

		temp=v.readline()

		if(temp[0]=="." or temp[0]=="("):
			structlib.append(temp[0:-1])

		if(start_position_1%12==0): #the sequence
			triplelist.append(temp[0:-1]) 

		if(start_position_2%12==0):#the normality
			triplelist.append(float(temp[11:-1]))
		if(start_position_3%12==0):#the normality
			triplelist.append(float(temp[33:-1])) 

		start_position_1+=1
		start_position_2+=1
		start_position_3+=1
	return structlib, triplelist[0::3], triplelist[1::3],triplelist[2::3]
structlib,seqs,edit,diversity=get_data()

structlookupdictionary = dict(zip(structlib, seqs))
# print(diversity)
# print(edit)

EditDict = {} 
for i in range(len(seqs)):
	EditDict[seqs[i]]=edit[i]
sortedEditDict = sorted(EditDict.items(), key=operator.itemgetter(1))
DivDict = {} 
for i in range(len(seqs)):
	DivDict[seqs[i]]=diversity[i]
sortedDivDict = sorted(DivDict.items(), key=operator.itemgetter(1))

sortedbyedit=[]
for i in range(len(seqs)):
	sortedbyedit.append(sortedEditDict[i][0])
	sortedbyedit.append(sortedEditDict[i][1])
	sortedbyedit.append(DivDict[sortedEditDict[i][0]])
sbe_seq=sortedbyedit[0::3]
sbe_div=sortedbyedit[2::3]
sbe_edi=sortedbyedit[1::3]

sortedbydiversity=[]
for i in range(len(seqs)):
	sortedbydiversity.append(sortedDivDict[i][0])
	sortedbydiversity.append(sortedDivDict[i][1])
	sortedbydiversity.append(EditDict[sortedDivDict[i][0]])
sbd_seq=sortedbydiversity[0::3]
sbd_div=sortedbydiversity[1::3]
sbd_edi=sortedbydiversity[2::3]


N = len(diversity)
x = range(1,len(seqs)+1)
# y=diversity
# y = sbd_div
# edit = sbd_edi


# print(sbd_div)
# print(sbe_div)
#edit= [i**2 for i in edit]
normalizededit = [number / (max(edit)) for number in edit] #normalize the edit dist vals so as to map it to the space of opacity

alphas=normalizededit

#print(pick_top_n_diversity(50))
comparediversityratings(pick_top_n_diversity(45),pick_top_n_editdist(45))

rgba_colors = np.zeros((len(alphas),4))
# for red the first column needs to be one
rgba_colors[:,0] = 0.75
# the fourth column needs to be your alphas
rgba_colors[:, 3] = alphas

def best_fit(X, Y):

    xbar = sum(X)/len(X)
    ybar = sum(Y)/len(Y)
    n = len(X) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(X, Y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in X]) - n * xbar**2

    b = numer / denum
    a = ybar - b * xbar

    print('best fit line:\ny = {:.2f} + {:.2f}x'.format(a, b))

    return a, b


# print(pick_top_n_editdist(25))

# print(pick_top_n_diversity(25))
y = sbe_div
x = sbe_edi
#plt.plot(np.poly1d(np.polyfit(x, y, 1)))
#plt.scatter(x,y)
a, b = best_fit(x, y)
#plt.scatter(a,b)
yfit = [a + b * xi for xi in x]
#plt.plot(x, yfit,'b')

mymodel = np.poly1d(np.polyfit(x, y, 4))
mymodel2 = np.poly1d(np.polyfit(x, y, 3))
mymodel3 = np.poly1d(np.polyfit(x, y, 2))
myline = np.linspace(min(x), max(x), 100)

fig, axs = plt.subplots(3, 2, sharey=False, tight_layout=True)

# We can set the number of bins with the `bins` kwarg

# plt.plot(myline, mymodel(myline),'r')
# plt.plot(myline, mymodel2(myline),'g')
# plt.plot(myline, mymodel3(myline),'y')
def calcstructmetrics(structlist):
	structdict={}
	lengthdict={}

	for structseq in structlist:
		k=0
		for structureatpoint in structseq:
			if(structureatpoint=="("):
				structdict.update({structseq : k})
				break
			k+=1

	for structseq in structlist:
		k=0
		stemlen=0
		for structureatpoint in structseq:
			if(structureatpoint=="("):
				break
			k+=1
		i=0
		for structureatpoint in structseq:
			try:
				if(structseq[i+1]=="."):
					temp=True
				else:
					temp=False
			except IndexError:
				temp=True
			if(structureatpoint==")" and temp):
				break
			i+=1
		lengthdict.update({structseq: i-k+1 })

	print("the range of structure length values ranged from: "+ str(min(lengthdict.values())) + " to " + str(max(lengthdict.values())))
	
	ranges=axs[0][0].hist(structdict.values(),ls='dotted',fc='#e8505b',ec='black',bins=range(min(structdict.values()),max(structdict.values())+2,1))
	axs[0][0].set_ylabel('Number Of Sequences')
	axs[0][0].set_xlabel("Structure Starting Position (histogram)")
	
	
	

	bins=axs[0][1].hist(lengthdict.values(),ls='dotted',fc='#f9d56e',ec='black',bins=range(min(lengthdict.values()),max(lengthdict.values())+2,1))
	axs[0][1].set_ylabel('Number Of Sequences')
	axs[0][1].set_xlabel("Structure length (histogram)")
	# sortedlengthtuplelist = sorted(lengthdict.items(), key=operator.itemgetter(1))
	sortedbins=[]
	binDict={}
	nvals=[]

	#Gets binning from the histogram itself
	binning=[]
	for index in range(len(bins[1])):
		try:
			nextdif=(bins[1][index+1])-(bins[1][index])
			tmp=[]
			for dx in range(0,nextdif):
				tmp.append(bins[1][index]+dx)
			binning.append(tmp)
		except IndexError:
			pass

	#print(binning)
	
	#collects seqs of same feature number into a list of lists. [[n,seq1,seq2],[n2,seq3,seq4]...etc] where n = feature size as a label
	for n in range(100):
		tmplst=[]
		tmplst.append(n)
		for item in range(len(lengthdict)):
			if(lengthdict.values()[item]==n):
				tmplst.append(lengthdict.keys()[item])
		nvals.append(tmplst)
	newnvals=[]
	for item in nvals:
		if (len(item)!=1):
			newnvals.append(item)
	#print(binning)
	#print(newnvals)
	#combine bins operation
	assem=[]
	for bingroups in binning:
		tmpgroup=[]
		for binnum in bingroups:
			for binseqgroup in newnvals:
				if(binseqgroup[0]==binnum):
					tmpgroup.append(binseqgroup)
		mastergroup=[]
		try:
			barxagg=0
			for i in range(len(tmpgroup)):
				barxagg+=(tmpgroup[i][0])
			barxagg=barxagg/float(str(len(tmpgroup))+".0")
			mastergroup.append(barxagg)
			for i in range(len(tmpgroup)):
				for item in tmpgroup[i][1:]:
					mastergroup.append(item)
			assem.append(mastergroup)
		except:
			pass
		
	newnvals=assem
	
	#gets the bin numbers (x axis)
	xval=[]
	for lengthgroups in newnvals:
		xval.append(lengthgroups[0])
	#print(xval)
	yval=[]
	#generates the gkm-div scores for each bin (avg)
	for lengthgroups in newnvals:
		aggregate=0
		for seq in lengthgroups[1:]:
			
			aggregate+=DivDict.get(structlookupdictionary.get(seq))
		yval.append(aggregate/(len(lengthgroups)-1))

	#print(yval)
	axs[1][1].bar(xval,yval,color="#5e6f64")
	axs[1][1].set_title("Bin gkm-svm Diversity of Above")
	axs[1][1].set_ylim(min(yval)-0.5*min(yval),max(yval)+0.5*max(yval))
	










	sortedbins=[]
	binDict={}
	nvals=[]

	#Gets binning from the histogram itself
	binning=[]
	for index in range(len(ranges[1])):
		try:
			nextdif=(ranges[1][index+1])-(ranges[1][index])
			tmp=[]
			for dx in range(0,nextdif):
				tmp.append(ranges[1][index]+dx)
			binning.append(tmp)
		except IndexError:
			pass
	
	
	#collects seqs of same feature number into a list of lists. [[n,seq1,seq2],[n2,seq3,seq4]...etc] where n = feature size as a label
	for n in range(100):
		tmplst=[]
		tmplst.append(n)
		for item in range(len(structdict)):
			if(structdict.values()[item]==n):
				tmplst.append(structdict.keys()[item])
		nvals.append(tmplst)
	newnvals=[]
	for item in nvals:
		if (len(item)!=1):
			newnvals.append(item)
	#print(binning)
	#print(nvals)
	#print(newnvals)
	#combine bins operation
	assem=[]
	for bingroups in binning:
		tmpgroup=[]
		for binnum in bingroups:
			for binseqgroup in newnvals:
				if(binseqgroup[0]==binnum):
					
					tmpgroup.append(binseqgroup)
		print("break")
		print(tmpgroup)
		print("break")
		mastergroup=[]
		barxagg=0
		for i in range(len(tmpgroup)):
			barxagg+=(tmpgroup[i][0])
		try:
			barxagg=barxagg/float(str(len(tmpgroup))+".0")
			mastergroup.append(barxagg)
			for i in range(len(tmpgroup)):
				for item in tmpgroup[i][1:]:
					mastergroup.append(item)
			assem.append(mastergroup)
		except:
			pass
		
	newnvals=assem

	#print(newnvals)
	#gets the bin numbers (x axis)
	xval=[]
	for lengthgroups in newnvals:
		xval.append(lengthgroups[0])
	#print(xval)
	yval=[]
	#generates the gkm-div scores for each bin (avg)
	for lengthgroups in newnvals:
		aggregate=0
		for seq in lengthgroups[1:]:
			
			aggregate+=DivDict.get(structlookupdictionary.get(seq))
		yval.append(aggregate/(len(lengthgroups)-1))

	#print(yval)
	axs[1][0].bar(xval,yval,color="#14b1ab")
	axs[1][0].set_title("Bin gkm-svm Diversity of Above")
	axs[1][0].set_ylim(min(yval)-0.5*min(yval),max(yval)+0.5*max(yval))












	sortedbins=[]
	binDict={}
	nvals=[]

	#Gets binning from the histogram itself
	binning=[]
	for index in range(len(ranges[1])):
		try:
			nextdif=(ranges[1][index+1])-(ranges[1][index])
			tmp=[]
			for dx in range(0,nextdif):
				tmp.append(ranges[1][index]+dx)
			binning.append(tmp)
		except IndexError:
			pass
	
	
	#collects seqs of same feature number into a list of lists. [[n,seq1,seq2],[n2,seq3,seq4]...etc] where n = feature size as a label
	for n in range(100):
		tmplst=[]
		tmplst.append(n)
		for item in range(len(structdict)):
			if(structdict.values()[item]==n):
				tmplst.append(structdict.keys()[item])
		nvals.append(tmplst)
	newnvals=[]
	for item in nvals:
		if (len(item)!=1):
			newnvals.append(item)
	#print(binning)
	#print(nvals)
	#print(newnvals)
	#combine bins operation
	assem=[]
	for bingroups in binning:
		tmpgroup=[]
		for binnum in bingroups:
			for binseqgroup in newnvals:
				if(binseqgroup[0]==binnum):
					
					tmpgroup.append(binseqgroup)
		print("break")
		print(tmpgroup)
		print("break")
		mastergroup=[]
		barxagg=0
		for i in range(len(tmpgroup)):
			barxagg+=(tmpgroup[i][0])
		try:
			barxagg=barxagg/float(str(len(tmpgroup))+".0")
			mastergroup.append(barxagg)
			for i in range(len(tmpgroup)):
				for item in tmpgroup[i][1:]:
					mastergroup.append(item)
			assem.append(mastergroup)
		except:
			pass
		
	newnvals=assem

	#print(newnvals)
	#gets the bin numbers (x axis)
	xval=[]
	for lengthgroups in newnvals:
		xval.append(lengthgroups[0])
	#print(xval)
	yval=[]
	#generates the gkm-div scores for each bin (avg)
	for lengthgroups in newnvals:
		aggregate=0
		for seq in lengthgroups[1:]:
			
			aggregate+=EditDict.get(structlookupdictionary.get(seq))
		yval.append(aggregate/(len(lengthgroups)-1))

	#print(yval)
	axs[2][0].bar(xval,yval,color="#14b1ab")
	axs[2][0].set_title("Bin edit distance of above")
	axs[2][0].set_ylim(min(yval)-0.5*min(yval),max(yval)+0.5*max(yval))










	sortedbins=[]
	binDict={}
	nvals=[]

	#Gets binning from the histogram itself
	binning=[]
	for index in range(len(bins[1])):
		try:
			nextdif=(bins[1][index+1])-(bins[1][index])
			tmp=[]
			for dx in range(0,nextdif):
				tmp.append(bins[1][index]+dx)
			binning.append(tmp)
		except IndexError:
			pass

	#print(binning)
	
	#collects seqs of same feature number into a list of lists. [[n,seq1,seq2],[n2,seq3,seq4]...etc] where n = feature size as a label
	for n in range(100):
		tmplst=[]
		tmplst.append(n)
		for item in range(len(lengthdict)):
			if(lengthdict.values()[item]==n):
				tmplst.append(lengthdict.keys()[item])
		nvals.append(tmplst)
	newnvals=[]
	for item in nvals:
		if (len(item)!=1):
			newnvals.append(item)
	#print(binning)
	#print(newnvals)
	#combine bins operation
	assem=[]
	for bingroups in binning:
		tmpgroup=[]
		for binnum in bingroups:
			for binseqgroup in newnvals:
				if(binseqgroup[0]==binnum):
					tmpgroup.append(binseqgroup)
		mastergroup=[]
		try:
			barxagg=0
			for i in range(len(tmpgroup)):
				barxagg+=(tmpgroup[i][0])
			barxagg=barxagg/float(str(len(tmpgroup))+".0")
			mastergroup.append(barxagg)
			for i in range(len(tmpgroup)):
				for item in tmpgroup[i][1:]:
					mastergroup.append(item)
			assem.append(mastergroup)
		except:
			pass
		
	newnvals=assem
	
	#gets the bin numbers (x axis)
	xval=[]
	for lengthgroups in newnvals:
		xval.append(lengthgroups[0])
	#print(xval)
	yval=[]
	#generates the gkm-div scores for each bin (avg)
	for lengthgroups in newnvals:
		aggregate=0
		for seq in lengthgroups[1:]:
			
			aggregate+=EditDict.get(structlookupdictionary.get(seq))
		yval.append(aggregate/(len(lengthgroups)-1))

	#print(yval)
	axs[2][1].bar(xval,yval,color="#5e6f64")
	axs[2][1].set_ylim(min(yval)-0.5*min(yval),max(yval)+0.5*max(yval))
	axs[2][1].set_title("Bin edit Diversity of Above")


	#print(structdict.values())
	# for index in range(len(bins[1])):
	# 	try:
			
	# 	except:
	# 		pass
	# 	print(bins[1][index])


	# sortedlengthdict=dict(sortedlengthtuplelist)
	
	# binvals=[]
	# print(bins[1])
	# for o in range(100):
	# 	if(len([k for k,v in sortedlengthdict.items() if v == o]) != 0):
	# 		binvals.append(o)
	# 		sortedbins.append([k for k,v in sortedlengthdict.items() if v == o])



	# # print(sortedbins)
	# aggvals=[]
	# binlens=[]
	# binavg=[]

	# for bins in sortedbins:
	# 	aggval=0
	# 	binlens.append(len(bins))
	# 	for indivstruct in bins:
	# 		aggval+=(DivDict.get(structlookupdictionary.get(indivstruct)))
	# 	aggvals.append(aggval)
	
	#print(binlens)
	# for item in range(len(aggvals)):
	# 	binavg.append(aggvals[item]/binlens[item])
	# # print(len(binvals))
	# # print(len(binavg))
	# #plt.bar(binvals,binavg,color="green")
	# axs[1][1].bar(binvals,binavg,color="green")
	# axs[1][1].set_title("Avg. Bin Diversity (bar graph)")
	#print(sortedlengthtuplelist)


	#print(sortedlengthdict)
	#print(dict.keys(sortedlengthdict))
	#print(sortedlengthdict)
	#for item in sortedlengthdict
	# print(len(lengthdict.values()))
	# print(len(structdict.values()))
calcstructmetrics(structlist)
print(len(structlist))

plt.tight_layout()
#plt.scatter(x, y, color=rgba_colors,s=10)
# plt.ylim(min(diversity)-0.25*(min(diversity)), max(diversity)+0.25*(max(diversity)))
# plt.xlabel("edit distance")
# plt.ylabel("gkm-svm normality")
# plt.title("Sequences Edit Dist. vs Diversity (Opacity~Edit Dist) Sorted by edit dist increasing")
plt.show()






