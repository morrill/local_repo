# import libraries
import numpy as np
import scipy.integrate
import csv
import itertools
import matplotlib.pyplot as plt
import os.path
#%% user functions 
# electrical measurement data functions
def read_in_data(filename):
    # read in .csv file and extract voltage and current data
    if os.path.isfile(filename):
        infile = open(filename)
        volt = np.zeros((100,1))
        current = np.zeros((100,1))
        f = csv.reader(infile)
        count = 0
        for row in itertools.islice(f,7,108):
            volt[count] = row[1]
            current[count] = row[13]
            count += 1
        infile.close()
        # caclulate slope of data set
        dI = np.gradient(current[:,0])
        slopeV = np.gradient(volt[:,0],dI)
        Res = np.mean(slopeV)
        return volt, current, slopeV, Res
    else:
        #print "File: "+ file[59:] +" not found in directory"
        volt = []
        current = []
        slopeV = []
        Res = []
        return volt, current, slopeV, Res
        
def plotmeas(plot_sample):
# Plots all IV measurements in the dictionary 'plot_sample'
    num_line = range(len(plot_sample))
    color_idx = np.linspace(0,1,len(num_line))
    keylst = plot_sample.keys()
    keylst.sort()
    title = keylst[0][:3]
    plt.figure()
    for i in num_line:
        if not any(plot_sample[keylst[i]]['I']):
            print("No data for sample: "+str(keylst[i]))
            continue
        xdata = plot_sample[keylst[i]]['I']
        ydata = plot_sample[keylst[i]]['V']
        plt.plot(xdata, ydata, color = plt.cm.gist_earth(color_idx[i]), linewidth=1)
    plt.xlabel("Current [A]")
    plt.ylabel("Voltage [V]")
    plt.title(title)
    plt.autoscale(tight=True)#plt.axis([0,max(xdata),0,1.1*max(ydata)])
    plt.grid()
       
 # profile functions       
def read_profile(filename):
    # reads in height profile data
    infile = open(filename)
    num_rows = sum(1 for row in infile)
    pos = np.zeros((num_rows,1))
    height = np.zeros((num_rows,1))
    infile = open(filename)
    f = csv.reader(infile)
    count = 0
    for row in f:
        pos[count] = row[0]
        height[count] = row[1]
        count += 1
    return height, pos 
    
def get_area(height,pos):
    # calculates profile data
    area = scipy.integrate.trapz(height,pos,axis=0) # Area in [nm^2]
    #print 'Area = : %.5f nm^2' %area
    return area
    
#%% plotting profile, MUST BE IN SAME WORKING DIRECTORY AS PROFILE FILES
plt.figure()
# Input the file names you want to plot, must be in current working directory
filenames = ['profile_13_5_1.csv','profile_13_5_2.csv','profile_13_5_3.csv']
linename = ['end 1','middle','end 2']
color_ind = np.linspace(0,1,3)
j = 0
for i in filenames:
    [height, pos] = read_profile(i)
    plt.plot(pos*1e-3,height*1e-3,color = plt.cm.brg(color_ind[j]),label = linename[j])
    j +=1
plt.legend(loc='best', bbox_to_anchor=(1,1))
#plt.autoscale() 
plt.ylim([0,7])
#plt.xlim([0,300])  
plt.title('Sample 13: 175 C for 15 min, 5 layers')
plt.xlabel("Cross line length [um]")
plt.ylabel("Height [um]")  
# get area
areas = np.zeros((3,1))
print('Area:')
j = 0
for i in filenames:
    [height, pos] = read_profile(i)
    areas[j] = get_area(height,pos)*(1e-3)**2
    print(linename[j] + ' = %.3f [um^2]' %areas[j])
    j += 1
ave_area = np.mean(areas)
var_area = np.var(areas)
std_area = np.std(areas)
print('Average for ' + filenames[0][8:12] + ' = %.3f [um^2]' %ave_area)
print('Standard deviation for ' + filenames[0][8:12] + ' = %.3f [um^2]' %std_area)
  
#%% Run to get all raw electrical measurement data, MUST BE IN RIGHT WORKING DIRECTORY
# generate list of sample ID numbers for file name
# current working directory must be where the data is
  
# for samples with multiple different layers
#samplePath = []
#samplenum = []
#sampleabrv = []
#for samp in np.arange(4,5):       # INPUT sample numbers, does not include upper bound
#    for lay in np.arange(4,5):      # layer numbers, does not include upper bound   
#        for meas in np.arange(1,6): # measurement number, does not include upper bound
#            path = os.path.join(os.getcwd(),str(samp)+"_"+str(lay)+"_"+str(meas)+".csv")            
#            samplePath.append(path)
#            #samplePath.append('G:\\Users\\Michael Orrill\\Documents\\python_notebooks\\2nd_run\\' + str(samp)+"_"+str(lay)+"_"+str(meas)+".csv")
#            samplenum.append(str(samp)+"_"+str(lay)+"_"+str(meas))
#        sampleabrv.append(str(samp)+"_"+str(lay))    
        
# for samples with no variation in layers
samplePath = []
samplenum = []
sampleabrv = []
heat = 'ht'
noheat = 'nht'
for samp in np.arange(1,11):    # INPUT sample numbers, does not include upper bound
    for meas in np.arange(1,4): # measurement number, does not include upper bound
        path = os.path.join(os.getcwd(),heat+str(samp)+"_"+str(meas)+".csv")            
        samplePath.append(path)
        #samplePath.append('G:\\Users\\Michael Orrill\\Documents\\python_notebooks\\2nd_run\\' + str(samp)+"_"+str(lay)+"_"+str(meas)+".csv")
        samplenum.append(str(samp)+"_"+str(meas))
    sampleabrv.append(str(samp))    
        
# put all data into dictionary called 'data' each measurement is accessable via its sampleID (1_2_3)
datahead = ['V', 'I', 'dV', 'R']
data = {head:{i:[] for i in datahead} for head in samplenum}
for i in range(len(samplenum)):
    [volt,current,slopeV,Res] = read_in_data(samplePath[i])
    data[samplenum[i]]['V'] = volt
    data[samplenum[i]]['I'] = current
    data[samplenum[i]]['dV'] = slopeV
    data[samplenum[i]]['R'] = Res
    # if a sample is not in the folder for read_in_data(), input NaN for all entries
    #if not any(data[samplenum[i]]['I']):
        #print "removing "+str(samplenum[i])+" from Dictionary 'data'"
     #   del data[samplenum[i]]
keylist = data.keys()
keylist.sort()
# create list of numbers for calling out certain samples
start = np.arange(0,len(keylist),5)
end = np.arange(5,len(keylist)+5,5)
# get overall maximum voltage measurement
maxVolt = np.zeros(len(keylist))
count = 0
for i in keylist:
    if not any(data[i]['V']):
        continue    
    maxVolt[count] = max(data[i]['V'])
    count += 1
OverallmaxVoltage = maxVolt.max() # this is the maximum V in the loaded data set
#%% Isolate resistance, calculate statistics and uncertainty, store in stats
currentR = np.zeros((5,1))
stats = {name:{u:[] for u in ['meanR','variance','SD','uR']} for name in sampleabrv}
u0c = 0.5*500e-9 # zero order error for current measurement
# Design stage Current uncertainty    
udc = np.sqrt(u0c**2 + (0.0002*1e-3+1.5e-6)**2) 
for i in range(len(end-1)): # cycle through all samples and making current_sample a dict with all data for all 5 meas
    current_sample = {k:data[k] for k in keylist[start[i]:end[i]]}
    samplist = current_sample.keys()
    samplist.sort()
    for k in [0,1,2,3,4]: # get the 5 resistance values for each measurement
        if not any(current_sample[samplist[k]]['I']):
            continue
        currentR[k] = current_sample[samplist[k]]['R']
        maxV = np.amax(current_sample[samplist[k]]['V'])
    # Zero order uncertainty for voltage depends on measurement range    
    if maxV <= 20e-3:
        u0v = 0.5*10e-9
        v_perc = 0.001
        v_acc = 150e-6
    elif 20e-3 < maxV <= 200e-3:
        u0v = 0.5*100e-9
        v_perc = 0.00012
        v_acc = 200e-6
    elif 200e-3 < maxV <= 2:
        u0v = 0.5*1e-6
        v_perc = 0.00012
        v_acc = 300e-6
    elif 2 < maxV <= 20:
        u0v = 0.5*10e-6
        v_perc = 0.00015
        v_acc = 1e-3
    elif 20 < maxV <= 200:
        u0v = 0.5*100e-6
        v_perc = 0.00015
        v_acc = 10e-3
    # Design stage # Design stage voltage uncertainty       
    udv = np.sqrt(u0v**2 + (v_perc*maxV+v_acc)**2)    
    # Uncertainty propegation to resistance measurement:
    stats[sampleabrv[i]]['uR'] = np.sqrt((-maxV/(1e-3)**2*udc)**2 + (1/1e-3*udv)**2) # combined uncertainty in IV measurement
    stats[sampleabrv[i]]['meanR'] = np.mean(currentR) 
    stats[sampleabrv[i]]['SD'] = np.std(currentR) # Standard deviation
    stats[sampleabrv[i]]['SDu'] = np.std(currentR)/np.sqrt(len(current_sample[samplist[0]]['I'])) # Standard uncertianty
    #stats[sampleabrv[i]]['variance'] = np.var(currentR)
    stats[sampleabrv[i]]['combined'] = np.sqrt(stats[sampleabrv[i]]['uR']**2 + stats[sampleabrv[i]]['SDu']**2)*2# combined standard uncertainty in R at 95% confidence

#%% Old code
''' This is all the old code that was used to build some functions/parts of scripts    
    
    
# Uncertainty analysis
# uncertainty from Keithley datasheet is systematic error
# uncertainty from varience of measured values is random error
# Keithley uncertainty [uk]:
u0c = 0.5*500e-9 # zero order error for current measurement
# Zero order uncertainty for voltage depends on measurement range
if maxVoltage <= 20e-3:
    u0v = 0.5*10e-9
    v_perc = 0.001
    v_acc = 150e-6
elif 20e-3 < maxVoltage <= 200e-3:
    u0v = 0.5*100e-9
    v_perc = 0.00012
    v_acc = 200e-6
elif 200e-3 < maxVoltage <= 2:
    u0v = 0.5*1e-6
    v_perc = 0.00012
    v_acc = 300e-6
elif 2 < maxVoltage <= 20:
    u0v = 0.5*10e-6
    v_perc = 0.00015
    v_acc = 1e-3
elif 20 < maxVoltage <= 200:
    u0v = 0.5*100e-6
    v_perc = 0.00015
    v_acc = 10e-3
# Design stage # Design stage voltage uncertainty       
udv = np.sqrt(u0v**2 + (v_perc*maxVoltage+v_acc)**2)    
# Design stage Current uncertainty    
udc = np.sqrt(u0c**2 + (0.0002*1e-3+1.5e-6)**2) 
# Uncertainty propegation to resistance measurement:
uR = np.sqrt((-maxVoltage/(1e-3)**2*udc)**2 + (1/1e-3*udv)**2)

    
    
    
 '''   
    
    
    
    
    
    