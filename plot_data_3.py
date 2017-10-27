import pandas as pd
import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import math as m

class case_dataset:  
    def __init__(self,inputfile,casename):
        self.inputfile = inputfile
        self.casename = casename
        self.filetag = inputfile.replace(".csv","") 
        self.dataframe = pd.read_csv(inputfile, skiprows=4)
        self.dataframe.columns = self.dataframe.columns.str.replace\
        ('Monitor Point: ','')
       
    def getaverage(self,N):
        self.average = (self.dataframe.tail(N)).mean()
        
    def filterdata(self,tags):
        self.filtered = {}
        for tag in tags:
            self.filtered[tag] = \
            self.average.filter(regex=tag).append(self.average.filter(regex=tag)[:1])   
                
    def plotpolar(self,*data,**options):
        toprint = options.get('save')
        
        # --- Default settings
        plt.figure(figsize=(6,6))
        fig = plt.subplot(111, polar=True)
        list_of_max = []
        print_tags = []
        
        for i in range(0,len(data)):
            fig.plot(theta, self.filtered[data[i]],'-o',\
                    label = data[i]+', '+self.casename)   
            fig.legend(bbox_to_anchor=(1.2, 1.1))
            print_tags.append( data[i]+'-'+self.casename)
            list_of_max.append(max(self.filtered[data[i]].values))   
        
        fig.grid(True)
        fig.set_theta_zero_location("N")
        fig.set_theta_direction(-1)
        fig.set_rlabel_position(70)
        fig.set_xlabel('Wall Shear [Pa]')
        fig.xaxis.set_label_coords(1.1, 0.7)
        fig.set_ylim([0,max(list_of_max)*1.1])
        fig.set_yticks([0,
                        round(max(list_of_max)*1,2)/4,
                        round(max(list_of_max),2)/2,
                        round(max(list_of_max),2)/4*3])  
        
        if (toprint == 1):
            filename = ("_".join(print_tags))+'.png'
            plt.savefig(filename,dpi=300,bbox_inches='tight')
       
    def plotpolarave(self,tag,**options):
        toprint = options.get('save')
        
        # --- Default settings
        plt.figure(figsize=(6,6))
        fig = pl.subplot(111, polar=True)
        list_of_max = []
       
        circle = pl.Circle((0, 0), self.filtered[tag].values.mean(),\
                           transform=fig.transData._b, color='red',\
                           linewidth=2,linestyle='--',fill=False)
        fig.add_artist(circle)
        # --- Plot
        fig.plot(theta, self.filtered[tag],'-o',\
                    label = tag+', '+self.casename)
        
        fig.legend(bbox_to_anchor=(1.2, 1.1))
        print_tags.append( tag+'-'+self.casename)
        list_of_max.append(max(self.filtered[tag].values))   
        
        fig.grid(True)
        fig.set_theta_zero_location("N")
        fig.set_theta_direction(-1)
        fig.set_rlabel_position(70)
        fig.set_xlabel('Wall Shear [Pa]')
        fig.xaxis.set_label_coords(1.1, 0.7)
        fig.set_ylim([0,max(list_of_max)*1.1])
        fig.set_yticks([0,
                        round(max(list_of_max)*1,2)/4,
                        round(max(list_of_max),2)/2,
                        round(max(list_of_max),2)/4*3])  
        
        if (toprint == 1):
            filename = self.casename+'-'+tag+'-polar_average.png'
            plt.savefig(filename,dpi=300,bbox_inches='tight')

    def plotpolarall(self,data,**options):
        toprint = options.get('save')
        
        # --- Default settings
        plt.figure(figsize=(8,8))
        fig = plt.subplot(111, polar=True)
        list_of_max = []
        
        for i in range(0,len(data)):
            fig.plot(theta, self.filtered[data[i]],'-o',\
                    label = data[i]+', '+self.casename)   
            fig.legend(bbox_to_anchor=(1.35, 1.2))
            print_tags.append( data[i]+'-'+self.casename)
            list_of_max.append(max(self.filtered[data[i]].values))   
        
        fig.grid(True)
        fig.set_theta_zero_location("N")
        fig.set_theta_direction(-1)
        fig.set_rlabel_position(70)
        fig.set_xlabel('Wall Shear [Pa]')
        fig.xaxis.set_label_coords(1.1, 0.7)
        fig.set_ylim([0,max(list_of_max)*1.1])
        fig.set_yticks([0,
                        round(max(list_of_max)*1,2)/4,
                        round(max(list_of_max),2)/2,
                        round(max(list_of_max),2)/4*3])  
        
        if (toprint == 1):
            filename = (self.casename)+'_all_polar.png'
            plt.savefig(filename,dpi=300,bbox_inches='tight')
        
    def plotmonitor(self,*tags,**options):
        n = options.get("n")             
        for tag in tags:
            plt.plot(self.dataframe[tag].loc[(self.dataframe.shape[0]-n):],\
                     label=tag)
        plt.xlabel('Iterations')
        plt.ylabel('Wall Shear [Pa]')
        plt.ylim(0,self.dataframe[tag].loc[(self.dataframe.shape[0]-n):].max()*1.1)
        plt.legend()
            
    def getstat(self,N):
        self.stat = pd.DataFrame()
        for tag in tags:
            self.stat[tag] = self.dataframe.tail(N).\
                filter(regex=tag).mean().values
        temp = self.stat.copy()
        self.stat.loc['mean'] = temp.mean()
        self.stat.loc['min'] = temp.min()
        self.stat.loc['max'] = temp.max()
        self.stat.loc['std'] = temp.std()
        self.stat.loc['diff'] = temp.max().values-temp.min().values
    
    def com(self,tags):
        segments = 8
        gamma = np.arange(segments) * 360/segments + 360/segments/2
        gamma = gamma*m.pi/180
        self.distances = pd.DataFrame()
        
        for tag in tags:
            ws = np.array(self.average.filter(like=tag).values)
            x = np.array(np.multiply(ws, np.cos(gamma)))
            y = np.array(np.multiply(ws, np.sin(gamma)))
            xy = np.matrix((x,y))
            xyave = xy.mean(axis=1)         
            self.distances.set_value(tag,'orig',\
                                     m.sqrt(np.sum(np.multiply(xyave,xyave))))
          
        self.distances['norm'] = (\
            (self.distances['orig']-self.distances['orig'].min())/ \
            (self.distances['orig'].max()-self.distances['orig'].min()))
           
        
    def kereszt(self,tags):
       self.kereszt = pd.DataFrame()
       self.keresztnorm = pd.DataFrame()
       self.keresztarany = pd.DataFrame()
       self.keresztaranynorm = pd.DataFrame()
       
       for tag in tags:
           self.kereszt[tag] = abs(self.stat[tag].iloc[0:4].values - \
               self.stat[tag].iloc[4:8].values)
           self.keresztnorm[tag] = (self.kereszt[tag]-self.kereszt.min().min()) \
               / (self.kereszt.max().max()-self.kereszt.min().min())
           for i in range(4):
               if (self.stat[tag].iloc[i] > self.stat[tag].iloc[i+4]):
                   value = self.stat[tag].iloc[i] / self.stat[tag].iloc[i+4]
               else:
                   value = self.stat[tag].iloc[i+4] / self.stat[tag].iloc[i]
               self.keresztarany.set_value(i,tag,value)
           
       self.keresztaranynorm = \
               (self.keresztarany-self.keresztarany.min().min()) / \
               (self.keresztarany.max().max()-self.keresztarany.min().min())
            
       
    def wave(self):
        self.wave = \
        self.keresztnorm.max()*0.5 + self.distances['norm']*0.5
    
    
    
print_tags = []
              
# --- Define tags 
tags = [
        'Felni111',
        'Felni112',
        'Felni121',
        'Felni122',
        'Felni211',
        'Felni212',
        'Felni221',
        'Felni222',
        'Felni311',
        'Felni312',
        'Felni321',
        'Felni322',
        ]

# --- Plot data
theta = np.array([0,45,90,135,180,225,270,315,360])
theta = theta + 22.5
theta = theta*m.pi/180


N = 2    
# --- Define objects    
vr0 = case_dataset("transient_temperature_vol_rate_0.csv",'VR0')
vr0.getaverage(2)
vr0.filterdata(tags)
vr0.getstat(2)
vr0.com(tags)
vr0.kereszt(tags)
vr0.wave()


dfwave = pd.DataFrame()
dfwave['0'] = vr0.wave
'''
plt.figure()
volrate210.wave.plot.bar()
plt.ylabel('Súlyozott átlag [-]')
plt.savefig('VR210_wave.png',dpi=300,bbox_inches='tight')
'''


