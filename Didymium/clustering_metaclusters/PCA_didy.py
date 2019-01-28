import ThreeDplot
import Functions as mylibrary
import varna_draw as Varna
import heatmap as HMP
import plotstruct as PSt
from ConfigParser import SafeConfigParser 
import time 
import collections 
from collections import defaultdict
#import Diana as Dn
#import draw3dPareto
#rna=Functions.ParseFasta(rnafile,"RNAfold")
#numberssamples=len(constraintes)

if __name__ == '__main__':
    #startime=time.time()
    #Get parameters from the config file
    config= SafeConfigParser()
    config.read("RNAstruct.Config")
    rnafile=config.get("Paths","rnafile")
    PathConstrainteFile=config.get("Paths","PathConstrainteFile")
    PathConstrainteFileShape=config.get("Paths","PathConstrainteFileShape")
    numberofsruct=config.get("Conditions","numberofsruct")
    Temperature=config.get("Conditions","Temperature")
    percent=int(config.get("Pareto","percent"))
    CutoffZcondition=float(config.get("Pareto","CutoffZcondition"))
    cutoff=config.get("Conditions","cutoffBasePairs")

    Psdotpath=config.get("Paths","Psdotpath")
    Matrixproba=config.get("Paths","Matrixproba")
    #constraintes=config.get("Conditions","Constraintes")
    Algorithm=config.get("Clustering","Algorithm")
    #constraintes=["didyNO","didy1M7", "didy1M7Mg", "didyNMIA", "didyNMIAMg", "didyBzCNMg", "didyCMCTMg",  "didyNaiMg","didyDMSMg","didyNai", "didyNMIAMgCE","didy1M7ILU3","didy1M7ILU","didy1M7ILU3Mg","didy1M7ILUMg"]
    constraintes=["didyNO","didy1M7", "didy1M7Mg", "didyNMIA", "didyNMIAMg", "didyBzCNMg", "didyCMCTMg",  "didyNaiMg","didyDMSMg","didyNai", "didyNMIAMgCE","didy1M7ILU3","didy1M7ILU","didy1M7ILU3Mg","didy1M7ILUMg"]


 
    #constraintes=["a83uDMS","a152uDMS","a223uDMS","c102gDMS","c219gDMS","c247gDMS","didyDMS","g101cDMS","g134cDMS","u123aDMS","u163aDMS"]
    rna=mylibrary.Parsefile(rnafile)[1]
    startimebig=time.time()

    print '**************Calculation of Eucledian distance between different BP dot plot conditions**********'
    #mylibrary.plotClusteringDistribution(numberofsruct,Matrixproba,len(rna))

    #mylibrary.DotplotRnaFold(Psdotpath,PathConstrainteFile,PathConstrainteFileShape)

    #mylibrary.Writeproba(Psdotpath,Matrixproba,constraintes,rna)

    mylibrary.plotClusteringDistribution(int(numberofsruct),Matrixproba,len(rna))
    print "Eucledian distance done"
