#!/usr/bin/env python
#!!!!!!!!!!!!!!!!Libraries*************************
# -*- coding: utf-8 -*-
import time
import shutil
import subprocess as sub
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import scipy
import pylab
from itertools import islice
from sklearn import cluster, datasets
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import load_svmlight_file
import cubehelix 
import seaborn as sns
#sns.set()

from matplotlib.colors import ListedColormap
#np.random.seed(0)
import collections 
from collections import defaultdict
import clusters_draw as Drawdistance
import varna_draw as Varna
#*************************!!!!!!!!!!!!!Containers*************************!!!
MatDist=collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
origine=collections.defaultdict(lambda: collections.defaultdict(lambda:0))
occ=collections.defaultdict(lambda:collections.defaultdict(lambda: 0))
Epsilon=collections.defaultdict(lambda:0)
#*************************!!!!!!!! *************************!!!!!!!!!!!!!!!!!!!
#create folder if it doesn't exist
def CreateFold(dir):
 	try:
    		os.stat(dir)
	except:
         	os.mkdir(dir)
#get all files with a specific extension from a specific path 
def GetListFile(mypath,extension):
	from os import listdir
	from os.path import isfile, join
	return  [f for f in listdir(mypath) if isfile(join(mypath, f)) and os.path.splitext(f)[1]== extension ]
#parse an input file, return  content line by line
def Parsefile2(Path):  
	lines=[] 
	with open(Path) as f:
    		for line in  islice(Path, 1, None):
        		lines.append([n for n in line.strip().split(' ')[0]])
    	return lines 

def Parsefile(Path):   
    	fileIn = open(Path,"r")
    	lines = fileIn.readlines()
    	fileIn.close()
    	return lines
def parseReactivityfile(fileinput):
	Reactvities=[]
	lines=Parsefile(fileinput)
	for it in range(len(lines)):
                Reactvities.append(lines[it].split("\t")[1])
	return Reactvities
def GetlinefromFile(Path,Linenumber):   
    	return Parsefile(Path)[Linenumber]

def MergeFiles(Path,fileslist,output):
	filenames = [Path+'/'+i for i in fileslist]
	with open(output, 'w') as outfile:
    		for fname in filenames:
        		with open(fname) as infile:
            			for line in islice(infile, 1, None):
                			outfile.write(line)
	return output
#*************************************Dealing with probabilities from dot plot ****************************
# generate dot plot base pairs
def DotplotRnaFold(dir,PathConstrainteFile,PathConstrainteFileShape):
         CreateFold(dir)
	 print PathConstrainteFileShape,"jj"
	 for filename in GetListFile(PathConstrainteFileShape,'.fa'):
		print filename,"lll"
		name=os.path.splitext(filename)[0]
		Input= PathConstrainteFileShape+"/"+filename
                ShapeFile=PathConstrainteFileShape+"/"+name+'Shape.txt'
		#print Input, ShapeFile,'hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh'
                output=dir+'/'+name
                command='RNAfold -p -d2 --noLP  ' + '<' +Input+  ' --shape '+ ShapeFile + '>'+ output
		os.system(command)
		#print command
         	shutil.move(name+'_dp.ps', dir+"/"+name+'_dp.ps')
		shutil.move(name+'_ss.ps', dir+"/"+name+'_ss.ps')


#def extract matrix values
def Writeproba( dir,Matrixproba,constraintes,rna):  
	CreateFold(Matrixproba)   
	for file in constraintes:
		PSPath =dir+"/"+file+"_dp.ps"
    		bpm = loadDotPlotPS(PSPath,"RNAfold")
    		dp = DotPlot(rna,bpm)
        	with open(Matrixproba +file.split('.')[0]+".proba","w") as o:
			for (i,j)in bpm.keys():
				o.write("%i\t%i\t%.6f\n"%(i,j,bpm[(i,j)]))
        	o.close()

# extract probabiliy values and build up the probability matrix 
B=collections.defaultdict(lambda: collections.defaultdict(lambda:collections.defaultdict(lambda: 0)))
from os import listdir
from os.path import isfile
from os.path import join
from sklearn import manifold
def Load_Probabilities(mypath):
    Dic={}       
    Beta = [f for f in listdir(mypath) if isfile(join(mypath, f)) and f.endswith('.proba')]
    for i,file in enumerate(Beta):
        Dic[i]=file.split(".proba")[0]
	lines=Parsefile(mypath+"/"+file)
	for it in range(len(lines)):
                lines[it]=lines[it].split("\t")
		B[i][int(lines[it][0])][int(lines[it][1])]=float(lines[it][2])	
		
    return Dic,B
# Calclate eucledian distance between 2 probability dot plot matrix
def Eucledian_distance(B,lenseq):
	dist=scipy.zeros([len(B),len(B)]) 
	for i in range(len(B)):
		for j in range(i+1,len(B)):
				for x in range(lenseq):
					for y in range(x+1,lenseq):
						#Just upper right half of the dotplot should be considered
						dist[i][j]+=math.pow(B[i][x][y]-B[j][x][y],2)
                                dist[j][i]=dist[i][j]
				#print "i",i,"j",j,math.sqrt(dist[i][j])
	return dist

import numpy as np
from matplotlib import pyplot as plt

from sklearn.datasets import make_biclusters
from sklearn.datasets import samples_generator as sg
from sklearn.cluster.bicluster import SpectralCoclustering
from sklearn.metrics import consensus_score


from sklearn.cluster.bicluster import SpectralBiclustering
def a():
    return []
def plotClusteringDistribution(lenconstraint,Folder_name,Lenrna):
        D=scipy.zeros([lenconstraint,lenconstraint])
	Dic,B=Load_Probabilities(Folder_name)
	clusters=defaultdict(a)
	#print Dic.values(),"ddddddddddddddddddddd"
        for element in Dic.keys():
		print Dic[element]
	#calculate the Eucledian distance between different entries
        D=Eucledian_distance(B,Lenrna)
        data= np.array(D)
        #rows=[Dic[element][-4:] for element in Dic.keys()]
        #columns=[Dic[element][4:] for element in Dic.keys()]
	rows=[Dic[element] for element in Dic.keys()]
        columns=[Dic[element] for element in Dic.keys()]
        #print "Clustering with kmeans for n=3"#n=8
        # by looking at the reuslt of spectralbilustering, it seeems that a subdivision of 8 is the appropriate one
        algorithm = cluster.MiniBatchKMeans(n_clusters=8)# 8
        #algorithm=cluster.AffinityPropagation(damping=.9, preference=None)
        algorithm.fit(data)
    	if hasattr(algorithm, 'labels_'):
        	y_pred = algorithm.labels_.astype(np.int)
    	else:
        	y_pred = algorithm.predict(data)
    	for i in range(len(y_pred)):
        	clusters[y_pred[i]].append(Dic[i])
	print clusters
        # score = consensus_score(model.biclusters_,(rows[:, row_idx], columns[:, col_idx]))
        # order ['NMIAMg', '1M7ILU', 'DMSMg', 'NO', 'NMIA','Nai','1M7ILUMg', '1M7ILU3Mg', 'CMCTMg','NaiMg', 'NMIAMgCE', 'BzCNMg','1M7', '1M7Mg' ,'1M7ILU3']
        model = SpectralBiclustering( random_state=0)
        model.fit(data)
        fit_data = data[np.argsort(model.row_labels_)]
        fit_data = fit_data[:, np.argsort(model.column_labels_)]
	#fit_data = data[np.array([6,0,8,13,11,4,10,2,14,1,12,3,9,7,5])]
        #fit_data = fit_data[:, np.array([6,0,8,13,11,4,10,2,14,1,12,3,9,7,5])]
        fig = plt.figure()
        ax = fig.add_subplot(111)
        orig=ax.matshow(data, cmap=plt.cm.Blues)
        fig.colorbar(orig)
        #plt.title("Base pairs eucledian distance between probing conditions") 
        ax.set_xticklabels(rows)
        ax.set_xticks(np.arange(len(rows)))# to show all labels
	ax.set_yticklabels(columns)
        ax.set_yticks(np.arange(len(columns)))
        #plt.show()
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        cx2 = cubehelix.cmap(reverse=True)#, "#FC8D62" ,"#8DA0CB", "#E78AC3", "#A6D854"
        #flatui = ["#66C2A5", "#95a5a6","#E78AC3", "#A6D854", "#FC8D62"]
        flatui=["#66C2A5","#457d6b","#365d50","#263e36","#17221e"]
        mycmap = ListedColormap(sns.color_palette(flatui).as_hex())

        b2=ax.matshow(fit_data, cmap=mycmap)
        #b2.ax.tick_params(labelsize=18) 
        cbar=fig.colorbar(b2,ticks=[0,50])
        cbar.set_label(label='Euclidean distance',size=12)

	for font_objects in cbar.ax.yaxis.get_ticklabels():
    		font_objects.set_size(20)
        #plt.title("Eucledian Distance matrix between conditions after bi-clustering") 
        #rows=[Dic[label][4:] for label in np.argsort(model.row_labels_)]
        #columns=[Dic[label][4:] for label in np.argsort(model.column_labels_)]
	rows=[Dic[label][4:] for label in np.argsort(model.row_labels_)]
        columns=[Dic[label][4:] for label in np.argsort(model.column_labels_)]
        #print rows
        ax.set_xticklabels(rows ,rotation=90)
        ax.set_xticks(np.arange(len(rows)))
	ax.set_yticklabels(columns, rotation_mode="anchor")
        ax.set_yticks(np.arange(len(columns)))
        ax.grid(False)
	# Turns off grid on the secondary (right) Axis.
	#ax.right_ax(False)
        plt.tick_params(axis='both', which='major', labelsize=13)
        plt.tight_layout() 
        plt.savefig("bi_clustering.pdf",format="pdf", dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None,transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None, metadata=None)
	
        
        for i in range(len(B)):
		
		print i,Dic[i], np.mean([ D[i][j] for j in range(len(B)) if i!=j]), np.min([ D[i][j] for j in range(len(B)) if i!=j])
		#%,[ D[i][j] for j in range(len(B))]
                print '\n'
	#print "Distance" , D
        # Clustering process with th plot
        adist = np.array(D)
	amax = np.amax(adist)
	adist /= amax
	mds = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=6)
	results = mds.fit(adist)
	coords = results.embedding_
	# plot results 
	fig = plt.figure()
	plt.subplots_adjust(bottom = 0.1)
	plt.scatter(coords[:, 0], coords[:, 1], marker = 'o',s=100,c="#66C2A5")# s for marker size, c for color
	k=0
	##print Dic.values()
        listoriented=['NaiMg', 'NMIAMgCE', 'BzCNMg', 'Nai', '1M7ILUMg', '1M7ILU3Mg', '1M7ILU3', '1M7', '1M7Mg', 'CMCTMg', 'NMIAMg', '1M7ILU', 'DMSMg', 'NMIA']
	for label, x, y in zip(listoriented, coords[:, 0], coords[:, 1]):
		k+=1
                Pos=(4,6)
                topbottom='bottom'
		if k%2==0:
			state="right"
		else:
			state="left"
		#if(label=="didyNaiMg"):
		#	state="left"
		#if(label in ["didyNMIAMgCE","didy(-)","didy1M7ILU3Mg","didy1M7ILU","didyBzCNMg"]):
			Pos=(4,-4)
                	topbottom='top'
    		plt.annotate(
			label[4:],
			xy = (x, y), xytext = Pos,
			textcoords = 'offset points', ha = state, va = topbottom, fontsize=11,
			bbox = dict(boxstyle = 'round,pad=0.3', fc = "#66C2A5", alpha = 0.3))
			#arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        plt.xlim(-0.4, 0.4)
   	plt.ylim(-0.7, 0.6)
        plt.tick_params(axis='both', which='major', labelsize=11)

	#plt.show()
	fig.savefig("Euclidean_distance_dot_plot_Matrix.pdf",format="pdf", dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None,transparent=False, bbox_inches=None, pad_inches=0.1,frameon=None, metadata=None)
#*************************************Base pairs Handling!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def ListBasePairsFromStruct(Struct):   #return dic={structure:[liste de pairs de base ],....}
    	lista=[]
        stack=[]
	for i in range (len(Struct)):# sequence length
		if Struct[i] == '(':# Opening base pair
			stack.append(i)
		elif Struct[i] == ')':# Closing base pair
			k=stack.pop()
		        lista.append((k, i))	    
    	return lista

# Parse an RNAsubopt file to extract Base pairs 
def GetBasePairsFromStructFile(faPath):   #return dic={structure:[liste de pairs de base ],....}
	#print faPath	
	DicStruct={}  
    	lines = Parsefile(faPath)
	#print lines 
    	SeqLen=len(lines[1])-1 
	#print SeqLen,"seq length"     
    	for j in range(len(lines)):
		DicStruct[j]=ListBasePairsFromStruct(lines[j].strip().split(' ')[0])
    	return len(lines),DicStruct 

def SortDictionary(Dico):
	 for key, value in sorted(Bo.iteritems(), key=lambda (k,v): (v,k)):
		print "%s: %s" % (key, value)
#Calculate distance between each 2 structures and generate a matrix 2 by 2 , stored in the file SVMLFile
def DistanceTwoBPlist(Struct1,Struct2):
	return len(set(Struct1).symmetric_difference(set(Struct2)) )
	#return len(set(Struct1).intersection(set(Struct2)) )
# we define redondante structutre as a global variable

Redondantestructure=collections.defaultdict(lambda:collections.defaultdict(lambda:False))
def DistanceStruct(StructFile,SVMlFile,numberofsruct,MFESnbrstruct, constrainte): 
	
        Redondantestructure1=[]   	
	#print 'on est la lololo'
    	DicStruct={}
	Dicnumberofsruct={}
	print constrainte
	for i in range(len(constrainte)-1):
		Dicnumberofsruct[constrainte[i]]=numberofsruct
	Dicnumberofsruct[constrainte[len(constrainte)-1]]=MFESnbrstruct
	#print "*************", Dicnumberofsruct
    	nb,DicStruct=GetBasePairsFromStructFile(StructFile)
	#print 'llll',   nb,DicStruct
    	for i in range(0,nb) :
		for j in range(i+1,nb):
			MatDist[i][j]= DistanceTwoBPlist(DicStruct[i],DicStruct[j])# difference symetrique entre deux structure i et j where each structure is defined as the number of pair de base in the list of pairs        
			if MatDist[i][j]==0:
				#print "hhhh", i,j 
				#print int(j/numberofsruct)
				if j not in Redondantestructure1:
					if j >numberofsruct*(len(constrainte)-1):
						Dicnumberofsruct[constrainte[len(constrainte)-1]]-=1		
					else:
						Dicnumberofsruct[constrainte[int(j/numberofsruct)]]-=1
					Redondantestructure1.append(j)
				
    			MatDist[j][i]=MatDist[i][j]
	print Redondantestructure1
	for elem in  Redondantestructure1:
		print "sos",elem
		if elem <numberofsruct*(len(constrainte)-1):
	 		ConditionNumber=int((elem)/numberofsruct)         		
		else:
			ConditionNumber=len(constrainte)-1
		StructureNumber=elem-ConditionNumber* numberofsruct	
                #print ConditionNumber,StructureNumber,
		Redondantestructure[constrainte[ConditionNumber]][StructureNumber]=True	
	# strore the distance matrix in the file SVMLFile
    	o=open (SVMlFile,"w") 
	for i in range(len(MatDist)):
		#if  i not in Redondantestructure1:
                o.write("%i\t"%(i+1))
		for j in range(len(MatDist)):
			if (i!=j) :
				o.write("%i:%.4f\t"%(j+1,MatDist[i][j]))	
		o.write("\n")
	o.close()
	
	#fig.savefig('Eucledian_distance_dot_plot_Matrix.png')
	print "Warning! redondante structures"
	#Redondantestructure1
	print "The sample size for each condition", Dicnumberofsruct
	
	return Redondantestructure1,Dicnumberofsruct,MatDist
	#print "order of the structure(s) with assigned constraint",Redondantestructure






#import matplotlib.colors as colors
#!!!!!!!!!!!!!!!!!!!plot3d after clustering!!!


	#figx = pickle.load(file('FigureObject.fig.pickle'))

	#figx.show() 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# To execute command line in the python script 

def Process(Command):
	p = sub.Popen(Command,stdout=sub.PIPE,stderr=sub.PIPE)
    	output, errors = p.communicate()
    	#print output
# Position correponding to the maximum value in a list
def MaxListIndex(Lista):
	return [i for i, j in enumerate(Lista) if j == max(Lista)]

# retrun cumulated value for a specific number of elements
def GetSumList(lista, nbr):
	CUMUL=np.cumsum(sorted(np.array(lista)))
        return CUMUL[nbr]
# return key for specific value in a dictionary
def keys_of_value(dic, value):
    	for key in dic:
        	if dic[key] == value:
       			return key
import operator
def Key_max_value(dico):
	return max(dico.iteritems(), key=operator.itemgetter(1))[0]
#*************************!!!!Structure handling*************************!!!!!!!!!!!!!

# The file containing the structure for different conditions is arranged as follow: for each condition in the consraint list, a number of structure is added to this file with repect to the order of the condition in the constraint list. Thus for a sampling of 20, the file will contain (20* nbr of conitions), ordred by the condition position in the lst of constraints.

# for a given structure characterized by a 'StructureNumber' return the condition represented by this structure
def GetOrigineofStructure(StructureNumber,numberofsruct):
	if (StructureNumber % numberofsruct!=0):
		return int(StructureNumber/numberofsruct)+1
        else:
        	return int(StructureNumber/numberofsruct)

def GetRandomValuesfromList(lista, maxLen):
	import random
	while len(lista)>maxLen:
  		lista.remove(random.choice(lista))
	return lista

#************************* Sampling step**********************************************
def StructSampling(PathConstrainteFile,PathConstrainteFileShape,numberofsruct,Tmpr):
	 dir='OutputSamples'+str(numberofsruct)
         CreateFold(dir)
         
	 for filename in GetListFile(PathConstrainteFile,'.fa'):
		#print filename
		Input= PathConstrainteFile+"/"+filename
		output=dir+'/'+os.path.splitext(filename)[0]
                os.system('RNAsubopt --noLP -p ' + str(numberofsruct) +' -s -T' + str(Tmpr)+ ' -C  <'+Input+  '>'+ output)
                #Command=['RNAsubopt', '--noLP', '-p', str(numberofsruct),' -s', '-C <'+Input,  '>'+ output]
		#Process(Command)
         
	 for filename in GetListFile(PathConstrainteFileShape,'.fa'):
		#print filename
		Input= PathConstrainteFileShape+"/"+filename
                ShapeFile=PathConstrainteFileShape+"/"+os.path.splitext(filename)[0]+'Shape.txt'
                output=dir+'/'+os.path.splitext(filename)[0]
		os.system('RNAsubopt --noLP -p ' + str(numberofsruct) +' -s -T' + str(Tmpr)+ ' --shape '+ShapeFile + '<' +Input+  '>'+ output)
	 return dir

def StructSamplingdelta(PathConstrainteFile,PathConstrainteFileShape,numberofsruct):
	 dir='OutputSamples'+str(numberofsruct)
         CreateFold(dir)
         
	 for filename in GetListFile(PathConstrainteFile,'.fa'):
		#print filename
		Input= PathConstrainteFile+"/"+filename
		output=dir+'/'+os.path.splitext(filename)[0]
                nbe=int(numberofsruct)/100
		#print nbe,"nombre supos d estructures"
                os.system('RNAsubopt --noLP -e ' + str(nbe) +' -s -C  <'+Input+  '>'+ output)
                #Command=['RNAsubopt', '--noLP', '-p', str(numberofsruct),' -s', '-C <'+Input,  '>'+ output]
		#Process(Command)
         
	 for filename in GetListFile(PathConstrainteFileShape,'.fa'):
		#print filename
		Input= PathConstrainteFileShape+"/"+filename
                ShapeFile=PathConstrainteFileShape+"/"+os.path.splitext(filename)[0]+'Shape.txt'
                output=dir+'/'+os.path.splitext(filename)[0]
		os.system('RNAsubopt --noLP -e ' + str(int(numberofsruct)/100) +' -s --shape '+ShapeFile + '<' +Input+  '>'+ output)
	 return dir
#***************************************Clustering*****************************

#***************************************Write clusters in an csv file******************
def ClustersDistributions(clusters, Filename,constraintes,numberofsruct,numberssamples):
	origine=collections.defaultdict(lambda: collections.defaultdict(lambda:0))
        it=collections.defaultdict(lambda: 0)
        o=open (Filename,"w")
        o.write("Cluster \t structures \t")
        for j in range(0,numberssamples):
        	o.write("constraint %i = %s \t"%(j+1, constraintes[j]))
        o.write("Number of structures \t Number of groups\n") 
	occ=SamplePresentInClusters(clusters,numberofsruct)
		 
        for elem in clusters:
                o.write("%i \t %s \t"%( elem+1,clusters[elem]))
		for j in range(1,numberssamples+1):
                         if (occ[elem][j]!=0):
				it[elem]+=1	
			 o.write("%i\t"%(occ[elem][j]))	
                o.write("%i\t%i\t"%(len(clusters[elem]),it[elem]))
                o.write("\n")
        
        o.write("Cluster(s) with  high number of  present conditions is(are) : %s"%([ v+1 for v in it.keys() if it[v]==max(it.values())]))
       
#*************************!!!! Get structure's energy values by calling RNAeval from Viennapackage!!!!!!!!!!!!!!!!!
def FromStructFiletoRNAEvalInput(StructFile,InputRNAeval,rna):
	
    	lines=Parsefile(StructFile)
        
        #print len(lino)
	
	#print len(lines), 'wht happend'
    	o=open (InputRNAeval,"w") # geneate intermediate file with sequence+strcuture , seq+strcture .... as the input format  to use RNAeval   
        # print "sdfqspkojr",len(lines)
    	for i in range(1,len(lines)):
	 	o.write("%s%s\t"%(rna,lines[i]))
   	o.close()

def RNAEvalCommand(InputRNAeval,OutputRNAeval): 
	#argListDict= ['RNAeval','<',InputRNAeval, '>',OutputRNAeval]
    	#Process(argListDict)
	os.system('RNAeval <'+InputRNAeval+ '>'+OutputRNAeval)
#StructFile contains the RNA sequence in the first line and list of correponding structures by line
def ENERGY_VALUES_STRUCTURES(StructFile,rna):
    	Energy=[]
    	#generate the rnaeval input file
    	FromStructFiletoRNAEvalInput(StructFile,"InputRNAeval",rna)
   	 # launch the RNaeval command
    	RNAEvalCommand("InputRNAeval","energyvalues")
    	# Parse the RNAevaloutput to extract energy values
    	lines=Parsefile("energyvalues")
    	for i in xrange(1,len(lines),2):
        	# i is the stucture number and 'lines[i].split(" ")[1][1:-2]' is  the  corresponding  energy value 		
		#print 'holla',(lines[i].split(" ")[1][1:-2])
		Energy.append(lines[i].split(" ")[1][1:-2])
    	return Energy

##Boltzmman energy from thermodyamic energy e  according to the formula B=exp^\frac{-e}{RT}
def BoltzmaanEnergy(Energy):
        T=37+273.15
        R=0.0019872370936902486
	return np.exp(-float(Energy)/float(100.*R*T))

#***********************************************************

#************************epsilon by condition
def EpsilonbyCondition(numberofsruct,constraintes,ConditionalBoltzman,percent):
	
	for ConditionNumber in range(len(constraintes)):
		Epsilon[ConditionNumber]=GetSumList(ConditionalBoltzman[constraintes[ConditionNumber]], numberofsruct[constraintes[ConditionNumber]]/percent)
		
	return Epsilon

def EpsilonbyCondition(constraintes,ZConditions,Percentage):
	
	for ConditionNumber in range(len(constraintes)):		
		#Epsilon[ConditionNumber]=Percentage*ZConditions[constraintes[ConditionNumber]]
		Epsilon[ConditionNumber]=Percentage
	print "Epsilon",Epsilon
	print "Zconditions",ZConditions
	return Epsilon
#*************************Get cardinal for coditions verifying Epsilon test
def GetCardinalConditions( clusters,ConditionalBoltzman,constraintes,numberofsruct,Epsilon):
	EnergiesbyCluster={}
        CardinalConditions={}
	
	for ClusterNumber in clusters:                
		# for each cluster, a sum overall boltzmaan energies for a given condition is calculated
		# a condition counts for the cardinal of existing constrainte if the sum overall its strcutures is greater than epsilon[condition]
		l=collections.defaultdict(lambda:0)
                dic={}
		for structure in clusters[ClusterNumber]:
                        temp=0
                        ConditionNumber=int((structure-1)/numberofsruct)
                        StructureNumber=(structure-1)-ConditionNumber* numberofsruct
			temp=ConditionalBoltzman[constraintes[ConditionNumber]][StructureNumber]
			l[ConditionNumber]+=temp# sum of B energies for a given condition
                        dic[structure]=temp
                EnergiesbyCluster[ClusterNumber]=dic
		for ConditionNumber in range(len(constraintes)):
			if l[ConditionNumber]< Epsilon[ConditionNumber]:
				del l[ConditionNumber]
		CardinalConditions[ClusterNumber]=len(l)
	return CardinalConditions,EnergiesbyCluster
#****************************************************************
def  CumulatedBoltzmanEnergiesbyCluster(clusters,ConditionalBoltzman,numberofsruct,constraintes): 
	cBE={}        
	for ClusterNumber in clusters:
                l=0.
		for structure in clusters[ClusterNumber]:
                        ConditionNumber=int((structure-1)/numberofsruct)
                        StructureNumber=(structure-1)-ConditionNumber* numberofsruct
			l+=ConditionalBoltzman[constraintes[ConditionNumber]][StructureNumber]
		cBE[ClusterNumber]=l
	return cBE
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!MEA structure
def MEA_algo(P_ss,BPp,rna):
	startime1=time.time()
	startime=time.time()
	nb=len(rna)
	W=collections.defaultdict(lambda: collections.defaultdict(lambda:0))
	for i in range(0,nb):
		W[i][i]=P_ss[i]
 	for d in range(1,nb):
  		for j in range(d,nb):
   			i=j-d;
			#print "ok++++++++",i,j,k
			
   			Res1=P_ss[i]+W[i+1][j]
   			Res2=P_ss[j]+W[i][j-1]

			# condition hairpin should be added?!!!
			#print CanonicalPair(rna[i],rna[j])
   			Res3=2*BPp[i][j]+W[i+1][j-1]
    			lista=[]
    			for k in range(i,j):
				#print "ok",i,j,k
     				lista.append(W[i][k]+W[k+1][j])
			#print "lista",lista
    			Res4=np.max(lista)
			
    			W[i][j]=np.max([Res1,Res2,Res3,Res4])
			
			#W[i][j]=np.max([P_ss[i]+W[i+1][j],P_ss[j]+W[i][j-1],2*BPp[i][j]+W[i+1][j-1],np.max([W[i][k]+W[k+1][j] for k in range(i,j)])])
		endtime=time.time()
		#print("End of 1 iteration %53f\t"%(endtime-startime))
		startime=endtime
	endtime1=time.time()
	#print("End of W matrix fill in %53f\t"%(endtime1-startime1))
	return W


def BackTracing(W,BPp,rna,P_ss,i,j,pair):
	
	#print i, j ,"Brffff what is the problem???"
	'''
        if W[i][j]==(2*BPp[i][j]+W[i+1][j-1])*CanonicalPair(rna[i],rna[j]) and W[i][j]!=0:
		print i,j, "those should be bas epairs!!"
	'''
	if i<j:
		#print i, j ,"what is the problem???",P_ss[j],W[i][j-1]
  		if W[i][j]==P_ss[i]+W[i+1][j]:
			#print "SOS1",i,j,W[i][j],P_ss[i],W[i+1][j]
   			BackTracing(W,BPp,rna,P_ss,i+1,j,pair)
                
  		elif  W[i][j]==P_ss[j]+W[i][j-1]:
			#print "SOS2",i,j,W[i][j],P_ss[j],W[i][j-1]
   			BackTracing(W,BPp,rna,P_ss,i,j-1,pair)
                
  		elif W[i][j]==(2*BPp[i][j]+W[i+1][j-1]):		
   			pair.append((i,j))
   			BackTracing(W,BPp,rna,P_ss,i+1,j-1,pair)
  		else:
   			for k in range(i,j):
    				if W[i][j]==W[i][k]+W[k+1][j]:
     					BackTracing(W,BPp,rna,P_ss,i,k,pair)
    					BackTracing(W,BPp,rna,P_ss,k+1,j,pair)
					break
	
 	return pair
def fromPairsToStruct(rna, Pairs):
	structure=["." for i in range(len(rna)-1)]
	for (i,j) in Pairs:
		structure[i]='('
		structure[j]=')'
	return "".join(structure)

def MEA(BPp,rna,pair):
	startime=time.time()
	P_ss=collections.defaultdict(lambda:0)
	W=collections.defaultdict(lambda: collections.defaultdict(lambda:0))
	for i in range(len(rna)):
		 P_ss[i]=1-sum([BPp[min(i,j)][max(i,j)] for j in range(len(rna))])	
	W=MEA_algo(P_ss,BPp,rna)
	endtime=time.time()
	print("End of MEA %53f\t"%(endtime-startime))
	#print BackTracing(W,BPp,rna,P_ss,0,len(rna)-1)
        #print fromPairsToStruct(rna, BackTracing(W,BPp,rna,P_ss,0,len(rna)-1,pair))
	return fromPairsToStruct(rna, BackTracing(W,BPp,rna,P_ss,0,len(rna)-1,pair)),BackTracing(W,BPp,rna,P_ss,0,len(rna)-1,pair)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def  BasePairsbyCluster(clusters,Structurefile,Boltzmann,numberofsruct,constrainte,cutoff): 
	ListBPbyCluster=collections.defaultdict() 
	ListBPbyStruct=collections.defaultdict(lambda:collections.defaultdict()) 
	Zcluster=collections.defaultdict(lambda:collections.defaultdict()) 
	BoltzmannOverPairsbyCluster=collections.defaultdict(lambda: collections.defaultdict(lambda:collections.defaultdict(lambda: 0)))
	Fiftypercent=collections.defaultdict()
	#print BoltzmannProba
        
	for ClusterNumber in clusters:
		liste=[]
		for structure in clusters[ClusterNumber]:
			
			ConditionNumber=int((structure-1)/numberofsruct)
		        StructureNumber=(structure-1)-ConditionNumber* numberofsruct
			liste.append(Boltzmann[constrainte[ConditionNumber]][StructureNumber])
                Zcluster[ClusterNumber]=sum(liste)
		
		#print 'cluster',ClusterNumber

		list1=[]
		for structure in clusters[ClusterNumber]:
			
			ConditionNumber=int((structure-1)/numberofsruct)
                        StructureNumber=(structure-1)-ConditionNumber* numberofsruct
			for (i,j) in ListBasePairsFromStruct(GetlinefromFile(Structurefile,structure-1)):
				BoltzmannOverPairsbyCluster[ClusterNumber][i][j] += Boltzmann[constrainte[ConditionNumber]][StructureNumber]/Zcluster[ClusterNumber]
			'''  
			print 'herin',ListBasePairsFromStruct(GetlinefromFile(Structurefile,structure-1)) 
			'''                       	 
			ListBPbyStruct[ClusterNumber][structure]=ListBasePairsFromStruct(GetlinefromFile(Structurefile,structure-1))
                        '''
			list1+=ListBPbyStruct[ClusterNumber][structure]
			'''
			

                #print "how does it look like", sorted(BoltzmannOverPairsbyCluster[ClusterNumber].values())
		'''
                with open("ProbaPairs_"+str(ClusterNumber), 'w') as outfile:
			for j in BoltzmannOverPairsbyCluster[ClusterNumber].keys():
				for i in BoltzmannOverPairsbyCluster[ClusterNumber].keys():
	             			if BoltzmannOverPairsbyCluster[ClusterNumber][i][j]>0:
		                     		#print "Base pairs for a cutoff of",cutoff , i,j,BoltzmannOverPairsbyCluster[ClusterNumber][i][j]
						outfile.write("%i \t %i \t %f \n"%( i,j,BoltzmannOverPairsbyCluster[ClusterNumber][i][j]))
		'''
                ListBPbyCluster[ClusterNumber]=list1
		Fiftypercent[ClusterNumber]=len(clusters)/2
	#print "pairs pairs",BoltzmannOverPairsbyCluster
	#print "Zcluster",Zcluster
	return ListBPbyStruct,ListBPbyCluster,Fiftypercent, BoltzmannOverPairsbyCluster
# count the occurence of each pairs, then filter them n a way to keep only those with an accurence above the fixed Threshold
from collections import Counter
# Counter (ListPairs) contains base pairs as key and its occurence as value
# .OccurenceBasePairs(args) return only base paairs with frequency above the Threshold
def OccurenceBasePairs(ListPairs, Threshold):
	return [(elem[0],elem[1],Counter(ListPairs)[elem]) for elem in Counter(ListPairs)  if Counter(ListPairs)[elem]>=Threshold]

#Cluster Diameters calculation
def ClustersDiameter(clusters,ListBPbystrcut):
	print "CLusters Diameters:"
	lista=[]
	for ClusterNumber in clusters:
		d=max([DistanceTwoBPlist(ListBPbystrcut[ClusterNumber][structure1],ListBPbystrcut[ClusterNumber][structure2])  for structure1  in clusters[ClusterNumber] for structure2 in clusters[ClusterNumber] ])
		print "Cluster number",ClusterNumber, "diameter", d
		lista.append(d)
	print "Average distance",np.mean(lista)
	return lista
def ClustersDistances(clusters,Boltzmanprobabilities,ListBPbystrcut,numberofsruct,condition):
	E=collections.defaultdict() 
	#print "clusters",clusters
	#print "CLusters Distances:"
	for ClusterNumber in clusters:
		liste=[]
		#print "boltzman proba",Boltzmanprobabilities
		for structure1  in clusters[ClusterNumber]:
			for structure2 in clusters[ClusterNumber]:
				
				# TODO Change the conversion to  be considered as fucntion!!! 
				ConditionNumber1=int((structure1-1)/numberofsruct)
         			StructureNumber1=(structure1-1) -ConditionNumber1* numberofsruct
				ConditionNumber2=int((structure2-1)/numberofsruct)
         			StructureNumber2=(structure2-1) -ConditionNumber2* numberofsruct
					
				#print structure1,structure2, "ffff",ConditionNumber1,StructureNumber1,ConditionNumber2,StructureNumber2
				#print Boltzmanprobabilities[condition[ConditionNumber1]][StructureNumber1]
				#print Boltzmanprobabilities[condition[ConditionNumber2]][StructureNumber2]
				#print DistanceTwoBPlist(ListBPbystrcut[ClusterNumber][structure1],ListBPbystrcut[ClusterNumber][structure2])
					
			
				liste.append(Boltzmanprobabilities[condition[ConditionNumber1]][StructureNumber1]*Boltzmanprobabilities[condition[ConditionNumber2]][StructureNumber2]*DistanceTwoBPlist(ListBPbystrcut[ClusterNumber][structure1],ListBPbystrcut[ClusterNumber][structure2]) )
		E[ClusterNumber]=2*sum(liste)/float(len(clusters[ClusterNumber])*(len(clusters[ClusterNumber])-1))
	print "Mean distance in the clusters:",E
	return E
def CentroidBycluster(clusters,StructFile,Boltzmann,numberofsruct,constrainte,cutoff,rna):
	'''
	ListBP=collections.defaultdict() 
	Fiftypercent=collections.defaultdict()
        CentroidBP=collections.defaultdict()
	listCentroidStructure=collections.defaultdict()
	ListBPbystrcut=collections.defaultdict(lambda:collections.defaultdict()) 
	'''
	E=collections.defaultdict() 
	mycentroid=collections.defaultdict() 
	listpairscentroid=collections.defaultdict(lambda: collections.defaultdict()) 
	Myproba=collections.defaultdict(lambda: collections.defaultdict(lambda:collections.defaultdict(lambda: 0)))
	ListBPbystrcut,ListBP,Fiftypercent,Myproba=BasePairsbyCluster(clusters,StructFile,Boltzmann,numberofsruct,constrainte,cutoff)
	#print Diameters
	ListDiameters=ClustersDiameter(clusters,ListBPbystrcut)
	E=ClustersDistances(clusters,Boltzmann,ListBPbystrcut,numberofsruct,constrainte)
	for ClusterNumber in clusters:
		pair=[]
		
		print "MEA algorithm is running for ClusterNumber",ClusterNumber
		mycentroid[ClusterNumber],listpairscentroid[ClusterNumber]=MEA(Myproba[ClusterNumber],rna,pair)
		
	
		#print "cetroids are calculated here"
		'''
		CentroidBP[ClusterNumber]=OccurenceBasePairs(ListBP[ClusterNumber], Fiftypercent[ClusterNumber])
		print 	CentroidBP[ClusterNumber]
                MINdist=min([DistanceTwoBPlist(CentroidBP[ClusterNumber],ListBPbystrcut[ClusterNumber][structure]) for structure in clusters[ClusterNumber] ])
		#print MINdist
		for structure in clusters[ClusterNumber]:
			print "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh",structure
		#print Myproba[ClusterNumber]
		
				listCentroidStructure[ClusterNumber]=[structure for structure in clusters[ClusterNumber] if  DistanceTwoBPlist(CentroidBP[ClusterNumber],ListBPbystrcut[ClusterNumber][structure])== MINdist]
		'''
	print "Centroids calculation is done"
	#herein we add teh dstance between clusters:
	print "Cluster distances"
	MatriceDistanceCentroids=scipy.zeros([len(clusters),len(clusters)])
	MatriceDistanceClustersEucld=scipy.zeros([len(clusters),len(clusters)])
	for ClusterNumber in clusters:
		for ClusterNumber2 in clusters:
			if ClusterNumber2 > ClusterNumber:
				
				l= DistanceTwoBPlist(listpairscentroid[ClusterNumber],listpairscentroid[ClusterNumber2])  
				#print "distance between clusters comparing the centroide's distances",l
				MatriceDistanceCentroids[ClusterNumber][ClusterNumber2]=l
				MatriceDistanceCentroids[ClusterNumber2][ClusterNumber]=l
				#print "distance between clusters comparing the means distances", ClusterNumber, ClusterNumber2, np.abs(E[ClusterNumber]-E[ClusterNumber2]),np.sqrt(abs(pow(E[ClusterNumber],2)-pow(E[ClusterNumber2],2)))
				l=np.sqrt(abs(pow(E[ClusterNumber],2)-pow(E[ClusterNumber2],2)))
				MatriceDistanceClustersEucld[ClusterNumber][ClusterNumber2]=l
				MatriceDistanceClustersEucld[ClusterNumber2][ClusterNumber]=l
				#print "distance between clusters compring the centroide's distances", ClusterNumber, ClusterNumber2, DistanceTwoBPlist(ListBPbystrcut[ClusterNumber][listCentroidStructure[ClusterNumber][0]],ListBPbystrcut[ClusterNumber2][listCentroidStructure[ClusterNumber2][0]])  
	Drawdistance.plotDisatnce(MatriceDistanceCentroids,clusters,"blue")
	Drawdistance.plotDisatnce(MatriceDistanceClustersEucld,clusters,"red")
	return mycentroid,E,MatriceDistanceCentroids,ListDiameters
#halte!!!!!!!!!!!!!!!!!!!!!!
def BoltzmanEnergies(constraintes,StructfileRepos,numberofsruct,MFESnbrstruct,rna):
        Energy=collections.defaultdict(lambda:collections.defaultdict())
        Boltzman=collections.defaultdict(lambda:collections.defaultdict())
        ConditionalBoltzman=collections.defaultdict(lambda:collections.defaultdict())
	#PS_ij=collections.defaultdict(lambda:collections.defaultdict())
        ZBolzman=collections.defaultdict(lambda:collections.defaultdict())
        
        for Condition in constraintes:
		
		FileStructure=StructfileRepos+'/'+Condition
        	Energy[Condition]=ENERGY_VALUES_STRUCTURES(FileStructure,rna)# list of energy values for the structures present in the Condition
		#print Condition
	#print Energy
        for Condition in constraintes:
		
 		if Condition=="MFES":
			Boltzman[Condition]=[BoltzmaanEnergy(Energy[Condition][i]) for i in range(MFESnbrstruct)]
		else:	
			listawithoutRedonddnace=[]		
			for i in range(numberofsruct):
				Boltzman[Condition][i]=BoltzmaanEnergy(Energy[Condition][i])
				if not(Redondantestructure[Condition][i]):# if the structure is not redundunt
					listawithoutRedonddnace.append(BoltzmaanEnergy(Energy[Condition][i]))	
			
                ZBolzman[Condition]=sum(listawithoutRedonddnace)# Partition function 
	import pickle
	#pickle.dump(Boltzman, file('Boltzmann.pickle', 'wb'))       	       
	# Relative Boltzmaan energy for each strcture/Condition 
        listall=[]
        for Condition in constraintes[:-1]: #to not count MFES
                lista=[]		
		for i in range(numberofsruct):
			if not(Redondantestructure[Condition][i]):
				lista.append(BoltzmaanEnergy(Energy[Condition][i])/ZBolzman[Condition])
			else:
				lista.append(0)# to solve the problem of the number of structure variation
		listall+=lista			
		ConditionalBoltzman[Condition]=lista
        #sumoverallconditions=sum(listall)
        '''
	for Condition in constraintes:
		for element in range(numberofsruct) :
			print Condition, element,ConditionalBoltzman[Condition][element],"ghd"
			PS_ij[Condition][element]=float(ConditionalBoltzman[Condition][element]/ sumoverallconditions)
	'''
	#print "hfghfdgdgfgdgdhfh",listall
	#print ConditionalBoltzman
	#print "when dividing", float(ConditionalBoltzman/sumoverallconditions)
	return Boltzman,ConditionalBoltzman,ZBolzman



#*************************!!!!!!!!!!!!!!!!clusters handling*************************!!!!!!!!!!
def FilterCluster(cluster,redundant):#(liststruct, listredndantstruct)
	return [elem for elem in cluster if elem not in redundant]
# count the occurence of present conditions in a given cluster
def SamplePresentInClusters(clusters,numberofsruct): 
        for ClusterNumber in clusters:		
		for StructureNumber in clusters[ClusterNumber]:			
			origine[ClusterNumber][StructureNumber]= GetOrigineofStructure(StructureNumber,numberofsruct)
        #calculate occurence of conditions present whithin a cluster               
        for ClusterNumber in clusters:
        	for ConditionInCluster in origine[ClusterNumber].values():
			occ[ClusterNumber][ConditionInCluster]=origine[ClusterNumber].values().count(ConditionInCluster)

 	return occ


#*************************!!!!!!!!!!!!!Pareto front*************************!!!!!!!!!!!!!!!!!!
def dominates(row, rowCandidate):
    	return all(r >= rc for r, rc in zip(row, rowCandidate))
# Return the dominating cluster with the 
# the input is a dicionary : data= { cluster number:[V1 value, V2 value],cluster number:[V1 value, V2 value]...}
def Pareto(Dico):
    	cleareCluster=[]
    	remaining = Dico.keys()

    	while remaining:
    		candidateKey = remaining.pop()
		candidate = Dico[candidateKey]

        	new_remaining = []
        	for other in remaining:
			if not dominates(candidate, Dico[other]):
				new_remaining.append(other)	   

        	if not any(dominates(Dico[other], candidate) for other in new_remaining):
            		cleareCluster.append(candidateKey)
        	remaining = new_remaining
		print len(remaining) 
    	#return cleareCluster,cleared 
	return cleareCluster

def dominant( c,d):
	return (c[0]>d[0] and c[1]>=d[1] and c[2]>=d[2] ) or (c[0]>=d[0] and c[1]>d[1] and c[2]>=d[2] ) or (c[0]>=d[0] and c[1]>=d[1] and c[2]>d[2] )

def dominee( c,d):
	return(c[0]<d[0] and c[1]<=d[1] and c[2]<=d[2] ) or (c[0]<=d[0] and c[1]<d[1] and c[2]<=d[2] ) or (c[0]<=d[0] and c[1]<=d[1] and c[2]<d[2] )
	
# Return the dominating cluster with the 
# the input is a dicionary : data= { cluster number:[V1 value, V2 value],cluster number:[V1 value, V2 value]...}
def Pareto2(Dico):
	Cmax=[]
	C= Dico.keys()
	while C:
		
		c=C.pop()
		#print c
		dominated=False
		for d in Cmax :
			#print Dico[c],Dico[d],dom(Dico[c],Dico[d])
			if (dominant(Dico[c],Dico[d])):
				Cmax.remove(d)
			elif (dominee(Dico[c],Dico[d])) :
			        dominated=True
		if not dominated:
			Cmax.append(c)
				
	return Cmax
#******************************************plot Clusters in function of (Boltzmaan cumulated Energy) & (Cardinal energies)
# something not normal is taking place here, the plot function confuses X and Y axis
def plotClustercBECard(clusternumber,cBE,cardinal,xlabelo,ylabelo,output):

	Labels=clusternumber
        X=cBE
	Y=cardinal
 	fig = plt.figure()
	fig.suptitle('Pareto front for clusters', fontsize=14, fontweight='bold')
	ax = fig.add_subplot(111)
	fig.subplots_adjust(top=0.85)
	ax.set_title('Cluster distribution')
	ax.set_xlabel(xlabelo)
	ax.set_ylabel(ylabelo)
        ax.set_ylim(bottom=0)	
        plt.axis([0,max(X)+1 ,0, max(Y)+np.mean(Y)])		
        ax.grid(True)
	for i in Labels:
		ax.text(X[i]+0.2,Y[i], i+1,fontsize=10,horizontalalignment='center',color='b')
        plt.plot(X,Y, 'r*')
        fig.savefig(output)

#**********************************************************Dot plot class
def loadDotPlotPS(path,tag):
    positions=[]
    if(tag=="RNAfold"):
	positions=[0,1,2,3]
    if(tag=="RNAalifold"):
	positions=[3,4,5,6]
    res = {}
    outTmp = open(path)
    for l in outTmp:
        data = l[:-1].split()
        if len(data) == positions[3]+1 and data[positions[3]]=="ubox":
            i = int(data[positions[0]])-1
            j = int(data[positions[1]])-1
            p = math.pow(float(data[positions[2]]),2.)
            res[i,j] = p
    outTmp.close()
    return res

class DotPlot:
    """Class for holding/producing base-pair probability matrices"""
    def __init__(self, rna , bpm = None):
        self.rna = rna[:]
        if bpm is None:
            # we will avoid this case to be sure that rnafold from the min works well
            self.bpm = self.runRNAFold()
        else:
            self.bpm = bpm
    def getSeq(self):
        return self.rna

    def getBPProb(self,i,j):
        if (i,j) in self.bpm:
            return self.bpm[i,j]
        else:
            return 0.

    def getUnpairedProbs(self):
        res = [1. for i in self.rna]
        for i,j in self.bpm:
            res[i] -= self.bpm[i,j]
            res[j] -= self.bpm[i,j]
        return res 

def DrawvarnaANdall(File1,Listclusters,CentroidStructure,numberofsruct,rna,Centroids_Energies,PathConstrainteFileShape):	
	c = canvas.canvas()	
	i=0
	with open (File1,"w") as OPt:
		lista=[]	
		for elem in Listclusters:
			i=i+1
			#print "The corresponding Boltzmann Energy structure", EnergiesbyCluster[elem]
			#Optimalstructure= mylibrary.Key_max_value(EnergiesbyCluster[elem])
			filename=File1+numberofsruct+"_"+str(elem)
		
			OPt.write("%s\n"%(CentroidStructure[elem]))
			#Heatmap
			lista+=ListBasePairsFromStruct(CentroidStructure[elem])

			# Varna call
			with open (filename,"w")as outt:
				outt.write(">HIV\n%s\n%s"%(rna,CentroidStructure[elem]))
			print " HIV  structure Centroid of the ",File1,elem,"with the energy value of",Centroids_Energies[elem-1], CentroidStructure[elem]
		        
			Varna.drawStructure(filename,PathConstrainteFileShape+"/HIV1M7Shape.txt",filename+".eps")
			c.insert(epsfile.epsfile(0, i*40, filename+".eps"))
	return c,lista	
