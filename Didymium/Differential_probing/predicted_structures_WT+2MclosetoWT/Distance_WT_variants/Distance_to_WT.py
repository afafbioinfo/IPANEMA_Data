import os, numpy as np
from os.path import isfile, join
from collections import OrderedDict
import math
def GetListFile(PathFile, FileExtension):
    return [os.path.splitext(f)[0] for f in os.listdir(PathFile) if isfile(join(PathFile, f)) and os.path.splitext(f)[1] == '.' + FileExtension]

def Parsefile(Path):
    reactivities=[]
    with open(Path) as infile:
    	for line in infile:
		reactivities.append(float(line.split("\t")[1][:-1]))
    return np.array(reactivities)
def Dist(list1,list2):
        diff=[]
	for elem1,elem2 in zip(list1,list2):
		diff.append(math.abs(elem1-elem2,2))
	return sum(diff)
ReactivityWT=Parsefile("didyDMS0Shape.txt")
#print ReactivityWT
DistancetoWT=OrderedDict()
for filz in GetListFile("Probing", "txt"):
        reactivity_variant=Parsefile(os.path.join("Probing", filz  + ".txt"))
	DistancetoWT[filz]=Dist(ReactivityWT,reactivity_variant)
	#print filz,DistancetoWT[filz]

print sorted(DistancetoWT.items(),key=lambda x:(x[1]))[:20]
#C192GShape', 0.5913185200000002), ('A238UShape', 0.5938680499999998), ('C147GShape', 0.8682818600000006), ('G109CShape', 0.94572358), ('G86CShape', 1.2722753199999994), ('G183CShape', 1.3308919299999995), ('C251GShape', 1.3675083799999999), ('U160AShape', 1.3937060699999992), ('U218AShape', 1.4823087700000004), ('G216CShape', 1.51366348), ('A173UShape', 1.5372835900000006), ('A175UShape', 1.6542395299999995), ('G65CShape', 1.6733013999999995), ('A145UShape', 1.7075460399999982), ('A130UShape', 1.7160515399999996), ('G235CShape', 1.7624028400000005), ('G101CShape', 1.800334529999998), ('A77UShape', 1.8369160700000002), ('A89UShape', 1.8807724799999983), ('C240GShape', 1.8883353400000005)]

