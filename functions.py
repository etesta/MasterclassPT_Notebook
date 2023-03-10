import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import warnings
warnings.filterwarnings('ignore')

def extract(filename, keys = ['PTV', 'CTV', '-GTV', 'Kidney_R', 'Stomach', 'SmallBowel',
                              'LargeBowel', 'Celiac', 'Liver', 'Heart', 'SpinalCord', 
                              'duodenum', 'Skin']):
    
    data = {} 
    MyKeys = []

    #known keys
    #keys = set(['PTV',                
            # 'CTV',
            # '-GTV',
            # 'Kidney_R',
            # 'Kidney_L',
            # 'Stomach',
            # 'SmallBowel',
            # 'LargeBowel',
            # 'Celiac',
            # 'Liver',
            # 'Heart',
            # 'SpinalCord',
            # 'duodenum',
            # 'Skin'])

    mat = scipy.io.loadmat(filename) #load file

    #dive into the main structure
    deep = mat['hgS_070000'][0, 0] # a Matlab 7 figure

    #All wanted data are deep inside the mat structure
    #for each organ data is inside a list of arrays with first key 'graph2d.lineseries'
    #[3:full array 'axis'][0:dive][0:dive][3:graph2d.lineseries list]
    graph2d = deep[3][0][0][3]

    #get number of keys in this list
    m = 0
    for key in graph2d:
        m+=1

    #iterate over each dataseries using k index and dive to retrieve x and y data
    k = 0
    while(k < m):
        #[i:select graph2d.lineseries][0:dive][0:dive][2:data for this graph2d.lineseries][0: dive][0: dive]
        arr = graph2d[k][0][2][0][0] 
        i = 0 #index of key in list, used to retreive the array of data which is right after the array containing the key
        for key in arr:
            #print('--- ', key)
            try: #array is a mess and sometimes key might be missing, the try statement prevent crashing
                if str(key[0]) in keys: #we often need to dive using [0] because each variable of the structure is actually an array
                    #print('--- ', key)
                    xdata = arr[i+1][0]
                    ydata = arr[i+2][0]
                    #print('xdata len : ', len(xdata))
                    #print('ydata len : ', len(ydata))
                    if(len(xdata) <= 1 or len(ydata) <= 1): #array is empty, dose is null 
                        #print('')#'no data')
                        break
                    else:
                        name = str(key[0])
                        if (name == 'CTV' or name == 'GTV' or name == 'PTV'):
                            name = '-' + name 
                        MyKeys.append(name) 
                        data[name] = {'x': xdata, 'y': ydata}
                i+=1
            except Exception as e:
                print(e)
        k+=1
    return MyKeys, data

def plot_AllHDV(data, keys, xmax, file):

    fig, ax = plt.subplots(figsize=(5,5))
    #fig, ax= plt.figure(figsize=(15,5))


    labels = []
    for key  in keys:
        plt.plot(data[key]['x'], data[key]['y'], label = key)
        
    #plt.legend(loc='center left', bbox_to_anchor=(1, .5))
    plt.legend()
    plt.xlabel('Dose (Gy)')
    plt.ylabel('Fraction of the volume (%)')
    plt.xlim(0,xmax)
    plt.ylim(0,105)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    plt.grid()
    
    plt.rc('legend',fontsize=15)
    fontsize = 15
    parameters = {'axes.labelsize': fontsize,
          'axes.titlesize': fontsize, 'xtick.labelsize':fontsize, 'ytick.labelsize':fontsize}
    plt.rcParams.update(parameters)   
    plt.savefig(file+'.png',bbox_inches='tight')
    plt.savefig(file+'.pdf',bbox_inches='tight')
    
    
def plot_HDV(data, organ):
    fig= plt.figure(figsize=(15,5))
    #fig, axs = plt.subplots(1, 2, figsize=(15, 5))
    plt.rc('legend',fontsize=15)
    
    parameters = {'axes.labelsize': 15,
          'axes.titlesize': 15, 'xtick.labelsize':15, 'ytick.labelsize':15}
    plt.rcParams.update(parameters)  
    plt.plot(data[organ]['x'], data[organ]['y'])
    

def IdealTumorHDV(data, PrescribedDose):
    ICurve=[]
    j=0
    length = len(data['-GTV']['x'])
    #length = get_size(data)
    #print("length : " , length)
    while(j<length):
        if(data['-GTV']['x'][j]<PrescribedDose):  
            ICurve.append(100)
        else:
            ICurve.append(0)
        j+=1
        #print(j)
    return ICurve

def get_TumorDoseDeviation(data, PrescribedDose):
    i=0
    count=0
    deviation=0
    tmp=0
    length = len(data['-GTV']['x'])
    ICurve = IdealTumorHDV(data, PrescribedDose)
    while(i<length):
        tmp=data['-GTV']['y'][i]-ICurve[i]
        if(np.abs(tmp) > 1):
            print(tmp)
            deviation+= tmp*tmp
            count+=1
        i+=1
    deviation=np.sqrt(deviation)/(count) 
    return deviation
    
def OAR_DoseMaxVolume(data, organ, Vmaxvalue):
    length = len(data[organ]['x'])
    i=0
    while(i<length-1 and data[organ]['x'][i]<Vmaxvalue): 
        i+=1
    return data[organ]['y'][i]

def Tumor_DoseMinVolume(data, organ, Vminvalue):
    length = len(data[organ]['x'])
    i=0
    while(i<length-1 and data[organ]['x'][i]<Vminvalue): 
        i+=1
    return 100-data[organ]['y'][i]