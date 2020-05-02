import time
f = open("FinalBatch.txt", "r")

def file_len(fname):
    with open(fname) as g:
        for i, l in enumerate(g):
            pass
    return i + 1
z=[]
print("total lines: "+str(file_len("FinalBatch.txt")))
for x in range(file_len("FinalBatch.txt")):
	temp=f.readline()
	if(temp[0]=="."):
		z.append(temp)


for item in z[0:]:
	#print(type(item))
	z.append(len(item))
globe=0.0
indx=0

for element in z[0:64222]:
	#print(indx)
	globe+=len(element)
	indx+=1
globe=float(globe/(len(z)))
print(globe)




