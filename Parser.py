import time
f = open("testrun.txt", "r")

def file_len(fname):
    with open(fname) as g:
        for i, l in enumerate(g):
            pass
    return i + 1
z=[]
b=[]
print("total lines: "+str(file_len("testrun.txt")))
for x in range(file_len("testrun.txt")):
	temp=f.readline()
	if(temp!=''):
		z.append(temp)
	f.readline()
	f.readline()
	temp=f.readline()
	if(temp!=''):
		z.append(temp)
	f.readline()
for x in z[0::2]:
	for y in z[1::2]:
		b.append(x+y[0:-2])


def comparestring(string1,string2):
	counter=0
	similarity=0
	for letter in string1:
		if letter == string2[counter]:
			similarity+=1
		counter+=1
	return similarity


def calcdiff(inlist):
	Dict={}
	enum=0
	#print(inlist)
	for basestr in inlist:
		listexceptstring=inlist[:enum]+inlist[enum+1:]
		#print(listexceptstring)
		a=0
		for elem in listexceptstring:
			a+=comparestring(basestr,elem)
		Dict.update( {basestr : a} )
		enum+=1
	return Dict
y=calcdiff(b)
print(y)







