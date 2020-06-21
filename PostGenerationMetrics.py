import operator
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

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
def get_data():#fully integrated function to generate the top n number of sequences by minimizing normality score within the domain and ranking them. Returns sorted dict least to greatest normality score
	v = open("Metrics.txt", "r")
	start_position_1=11
	start_position_2=2
	start_position_3=3
	listofseqs=[]
	listofnormvals=[]
	triplelist=[]
	a=0
	for x in range(file_len("Metrics.txt")):

		temp=v.readline()

		if(start_position_1%12==0): #the sequence
			triplelist.append(temp[0:-1]) 

		if(start_position_2%12==0):#the normality
			triplelist.append(float(temp[11:-1]))
		if(start_position_3%12==0):#the normality
			triplelist.append(float(temp[33:-1])) 

		start_position_1+=1
		start_position_2+=1
		start_position_3+=1
	return triplelist[0::3], triplelist[1::3],triplelist[2::3]
seqs,edit,diversity=get_data()
#print(seqs)
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
y = sbd_div
edit = sbd_edi

print(sbd_div)
print(sbe_div)
edit= [i**2 for i in edit]
normalizededit = [number / (max(edit)) for number in edit] #normalize the edit dist vals so as to map it to the space of opacity

alphas=normalizededit


rgba_colors = np.zeros((len(alphas),4))
# for red the first column needs to be one
rgba_colors[:,0] = 0.75
# the fourth column needs to be your alphas
rgba_colors[:, 3] = alphas

plt.scatter(x, y, color=rgba_colors,s=10)
plt.ylim(min(diversity)-0.25*(min(diversity)), max(diversity)+0.25*(max(diversity)))
plt.xlabel("sequence ID")
plt.ylabel("gkm-svm normality")
plt.title("Sequences Edit Dist. vs Diversity (Opacity~Edit Dist) Sorted by normality increasing, squared to emphasize outliers")
plt.show()






