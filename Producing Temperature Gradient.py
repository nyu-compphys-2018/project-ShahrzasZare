#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  4 11:40:16 2018

@author: shabani_lab
"""

### libraries




import numpy as np
from random import random
import matplotlib.pyplot as plt




### constants and initializations which we may want to change in differen experiments



    
    ### more general constants

k = 1 
L = [10.0,10.0] #size of the room
number_of_microparticles = 1000
number_of_regions = 10
deltaX = L[0]/number_of_regions
deltaT = 0.005

    ### flags

thermostat_flag = 1 # on (1) or off (0)
temperature_gradient_flag = 1     # temperature gradient flag: 0: uniform ; 1: linear ; 2: parabolic

    ###constants related to the classes

micro_mass = 1
macro_mass = 1
radius_of_the_macroparticle = 1.0

    ### temperature

location = np.zeros(number_of_regions) #location of the reigon along x axis     
Temperature = np.zeros(number_of_regions) #temperature of the reigon along x axis
min_temperature = 20
max_temperature = 200

def temperature_initialization(temperature_gradient_flag,min_temperature,max_temperature):
    for i in range(int(number_of_regions)):
        location[i] = i*deltaX +deltaX/2
    ### temperature gradient flag:
        ### 0 : uniform
        ### 1 : linear
        ### 2 : parabolic
    if (temperature_gradient_flag==0):
    ### uniform temperature gradient
        for i in range(int(number_of_regions)):
            Temperature[i] = (max_temperature+min_temperature)/2
    else:
        if(temperature_gradient_flag==1):
        ### linear temperature gradient
            for i in range(int(number_of_regions/2)):
                Temperature[i] = min_temperature + (max_temperature-min_temperature)*location[i]/(L[0]*0.5)
            j = (number_of_regions/2)-1
            for i in range(int(number_of_regions/2),int(number_of_regions)):
                Temperature[i] = Temperature[int(j)]
                j -= 1
        else:
            if(temperature_gradient_flag==2):
                ### parabolic temperature gradient
                for i in range(int(number_of_regions)):
                    Temperature[i] = ((location[i]-0.5*L[0])**2)*(max_temperature+min_temperature)/(0.5*L[0])**2 + max_temperature   







###functions
                    
                    
                    

    ### thermostat and measurments


def thermostat(microparticles,temperature,density,kinetic_energy):  
    
    for i in range( int(number_of_regions) ):
        ### measurment
        temperature[i] = 0.0
        density[i] = 0.0000001 #0.0000001 is set to prevent deviding by zero
        kinetic_energy[i] = 0.0000001 #0.0000001 is set to prevent deviding by zero
        alpha = 0.0 #the coefficient which correct the velocities
        the_particles = np.zeros(number_of_microparticles) #number of microparticles is the maximum possible number of the particles which are located in the reigen, but the actual number is density-0.0000001
        ### find the microparticles which are located in the riegon
        for j in range(number_of_microparticles):
            if( microparticles[j].x<=(location[i]+deltaX/2) and microparticles[j].x>=(location[i]-deltaX/2) ):
                v2 = microparticles[j].vx**2 + microparticles[j].vy**2
                kinetic_energy[i] += 0.5*microparticles[j].m*v2
                the_particles[int(density[i])] = j
                density[i] += 1
                microparticles[j].region = i
        ### calculate temperatures for the reigon            
        temperature[i] = kinetic_energy[i]/(density[i]*k)
        alpha = np.sqrt(density[i]*k*Temperature[i]/kinetic_energy[i])
        ### corrrect the velocities
        for j in range(int(density[i])):
            microparticles[int(the_particles[j])].vx *= alpha
            microparticles[int(the_particles[j])].vy *= alpha
            

def measurments(microparticles,temperature,density,kinetic_energy):  
    
    for i in range( int(number_of_regions) ):
        ### measurment
        temperature[i] = 0.0
        density[i] = 0.0000001 #0.0000001 is set to prevent deviding by zero
        kinetic_energy[i] = 0.0000001 #0.0000001 is set to prevent deviding by zero
        ### find the microparticles which are located in the riegon
        for j in range(number_of_microparticles):
            if( microparticles[j].x<=(location[i]+deltaX/2) and microparticles[j].x>=(location[i]-deltaX/2) ):
                microparticles[j].region = i
                v2 = microparticles[j].vx**2 + microparticles[j].vy**2
                kinetic_energy[i] += 0.5*microparticles[j].m*v2
                density[i] += 1
        ### calculate temperatures for the reigon            
        temperature[i] = kinetic_energy[i]/(density[i]*k)           
    
        
    ### separate measurments


def calculate_density(density):
    
    for i in range( int(number_of_regions) ):
        density[i] = 0.0 
        ### find the microparticles which are located in the riegon
        for j in range(number_of_microparticles):
            if( microparticles[j].x<=(location[i]+deltaX/2) and microparticles[j].x>=(location[i]-deltaX/2) ):
                density[i] += 1
                
    
def calculate_temperature(temperature):
    
    for i in range( int(number_of_regions) ):
        temperature[i] = 0.0
        density = 0.0000001 #0.0000001 is set to prevent deviding by zero 
        kinetic_energy = 0.0
        ### find the microparticles which are located in the riegon
        for j in range(number_of_microparticles):
            if( microparticles[j].x<=(location[i]+deltaX/2) and microparticles[j].x>=(location[i]-deltaX/2) ):
                v2 = microparticles[j].vx**2 + microparticles[j].vy**2
                kinetic_energy += 0.5*microparticles[j].m*v2
                density += 1
        ### calculate temperatures for the reigon            
        temperature[i] = kinetic_energy/(density*k)
        temperature_error[i] = Temperature[i]-temperature[i] 

        
def calculate_kinetic_energy(kinetic_energy):
    
    for i in range( int(number_of_regions) ): 
        kinetic_energy[i] = 0.0
        ### find the microparticles which are located in the riegon
        for j in range(number_of_microparticles):
            if( microparticles[j].x<=(location[i]+deltaX/2) and microparticles[j].x>=(location[i]-deltaX/2) ):
                v2 = microparticles[j].vx**2 + microparticles[j].vy**2
                kinetic_energy[i] += 0.5*microparticles[j].m*v2        


def calculate_momentumX():
    
    momentumX = 0.0
    for i in range( int(number_of_microparticles) ): 
        momentumX += microparticles[i].m*microparticles[i].vx
    return momentumX


def calculate_momentumY():
    
    momentumY = 0.0
    for i in range( int(number_of_microparticles) ): 
        momentumY += microparticles[i].m*microparticles[i].vy
    return momentumY




    ### update


    
    
def update(microparticles,temperature,density,kinetic_energy):       
    for i in range(number_of_microparticles):    
        microparticles[i].x += deltaT*microparticles[i].vx
        microparticles[i].y += deltaT*microparticles[i].vy
        microparticles[i].check_event()
    if(thermostat_flag==1):
        thermostat(microparticles,temperature,density,kinetic_energy)




### classes
   
    
    
    
class Microparticle:    
    
    x = 0.0 #location of each macroparticle
    y = 0.0 #location of each macroparticle
    vx = 0.0 #velocity of each macroparticle
    vy = 0.0 #velocity of each macroparticle
    fx = 0.0 #force of each macroparticle
    fy = 0.0 #force of each macroparticle
    m = 0.0 #mass        
    location_assigned = 0.0
    region = 0            
    
    def check_event(self):    
        if(self.x<0):
            self.x = L[0]
        if(self.x>L[0]):
            self.x = 0
        if(self.y<0):
            self.y = L[1]
        if(self.y>L[1]):
            self.y = 0                    




### objects and variables
            
            

            
    ### define some quantities which we want to measure for each reigne along x axis
    
temperature = np.zeros(number_of_regions) #actual temperature of each reigen
temperature_error = np.zeros(number_of_regions) #the difference between actual and expected temperature of each reigen
density = np.zeros(number_of_regions) #number of microparticles in each reigon 
kinetic_energy = np.zeros(number_of_regions) #kinetic energy of each reigon

trackX = []
trackY = []
Energy = []
MomentumX = []
MomentumY = []
Time = []
        
    ### mobile objects
            
microparticles = [Microparticle() for i in range(number_of_microparticles)]




### main




#
    ### initialization
    
temperature_gradient_flag = 1
temperature_initialization(temperature_gradient_flag,20,800)
    
        ###location
        
for i in range(number_of_microparticles): 
    microparticles[i].m = micro_mass
    switch = 0
    while (switch==0): # I do not want the microparticles to be in same place at first step
        microparticles[i].x = random()*L[0]
        microparticles[i].y = random()*L[1]
        switch = 1
        for j in range(i):
                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
                    switch = 0        
    for j in range(number_of_regions):
        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
            microparticles[i].region = j
    
    ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls

for i in range(number_of_microparticles):
    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
    microparticles[i].vx = np.random.normal(0,Var)
    microparticles[i].vy = np.random.normal(0,Var)
#    

#fig2 = plt.figure()
#calculate_density(density)
#plt.plot(location,density,'b',label='linear temperature gradient: begining')
#for i in range(1000):
#    update(microparticles,temperature,density,kinetic_energy)
#calculate_density(density)
#plt.plot(location,density,'r',label='linear temperature gradient:after 100 time steps')
#fig2.show()




fig = plt.figure()
for i in range(50):
    for i in range(number_of_microparticles):
        color_code1 = Temperature[int(microparticles[i].region)]/max(Temperature)
        microparticle_show = plt.Circle((microparticles[i].x,microparticles[i].y),0.05,color=(color_code1, 0.0, 0.0))
        fig = plt.gcf()
        f = fig.gca()
        f.add_artist(microparticle_show)
    plt.pause(0.000001)
    fig.clf()    
    plt.xlim(0, L[0])
    plt.ylim(0, L[1])
    plt.gca().set_aspect('equal', adjustable='box')
    update(microparticles,temperature,density,kinetic_energy) 
fig.show()





#plt.figure()
#plt.title('Temperature of 20 Regions for 2000 microparticles')
#plt.xlabel('x')
#plt.ylabel('temperature')


#    ### initialization
#    
#temperature_gradient_flag = 0 
#temperature_initialization(temperature_gradient_flag,20,20)
#plt.plot(location,Temperature,label = 'expected uniform temperature = 20')
#    
#        ###location
#        
#for i in range(number_of_microparticles): 
#    microparticles[i].m = micro_mass
#    switch = 0
#    while (switch==0): # I do not want the microparticles to be in same place at first step
#        microparticles[i].x = random()*L[0]
#        microparticles[i].y = random()*L[1]
#        switch = 1
#        for j in range(i):
#                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
#                    switch = 0        
#    for j in range(number_of_regions):
#        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
#            microparticles[i].region = j
#    
#    ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls
#
#for i in range(number_of_microparticles):
#    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
#    microparticles[i].vx = np.random.normal(0,Var)
#    microparticles[i].vy = np.random.normal(0,Var)
#    
#
#
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0,0,0.4),label='uniform temperature: begining')
#for i in range(100):
#    update(microparticles,temperature,density,kinetic_energy)
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0,0,0.8),label='uniform temperature: after 100 time steps')
#
#
#
#    ### initialization
#    
#temperature_gradient_flag = 0 
#temperature_initialization(temperature_gradient_flag,100,100)
#plt.plot(location,Temperature,label = 'expected uniform temperature = 100')
#    
#        ###location
#        
#for i in range(number_of_microparticles): 
#    microparticles[i].m = micro_mass
#    switch = 0
#    while (switch==0): # I do not want the microparticles to be in same place at first step
#        microparticles[i].x = random()*L[0]
#        microparticles[i].y = random()*L[1]
#        switch = 1
#        for j in range(i):
#                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
#                    switch = 0        
#    for j in range(number_of_regions):
#        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
#            microparticles[i].region = j
#    
#    ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls
#
#for i in range(number_of_microparticles):
#    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
#    microparticles[i].vx = np.random.normal(0,Var)
#    microparticles[i].vy = np.random.normal(0,Var)
#    
#
#
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0.4,0,0),label='uniform temperature: begining')
#for i in range(100):
#    update(microparticles,temperature,density,kinetic_energy)
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0.8,0,0),label='uniform temperature: after 100 time steps')
#



#    ### initialization
#    
#temperature_gradient_flag = 1 
#temperature_initialization(temperature_gradient_flag,min_temperature,max_temperature)
#plt.plot(location,Temperature,label = 'expected linear temperature gradient')
#    
#        ###location
#        
#for i in range(number_of_microparticles): 
#    microparticles[i].m = micro_mass
#    switch = 0
#    while (switch==0): # I do not want the microparticles to be in same place at first step
#        microparticles[i].x = random()*L[0]
#        microparticles[i].y = random()*L[1]
#        switch = 1
#        for j in range(i):
#                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
#                    switch = 0        
#    for j in range(number_of_regions):
#        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
#            microparticles[i].region = j
#    
#    ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls
#
#for i in range(number_of_microparticles):
#    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
#    microparticles[i].vx = np.random.normal(0,Var)
#    microparticles[i].vy = np.random.normal(0,Var)
#    
#
#
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0.4,0,0),label='linear temperature gradient: begining')
#for i in range(500):
#    update(microparticles,temperature,density,kinetic_energy)
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0.8,0,0),label='linear temperature gradient: after 500 time steps')



#
#    ### initialization
#    
#temperature_gradient_flag = 2 
#temperature_initialization(temperature_gradient_flag,min_temperature,max_temperature)
#plt.plot(location,Temperature,label = 'expected parabolic temperature gradient')
#    
#        ###location
#        
#for i in range(number_of_microparticles): 
#    microparticles[i].m = micro_mass
#    switch = 0
#    while (switch==0): # I do not want the microparticles to be in same place at first step
#        microparticles[i].x = random()*L[0]
#        microparticles[i].y = random()*L[1]
#        switch = 1
#        for j in range(i):
#                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
#                    switch = 0        
#    for j in range(number_of_regions):
#        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
#            microparticles[i].region = j
#    
#    ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls
#
#for i in range(number_of_microparticles):
#    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
#    microparticles[i].vx = np.random.normal(0,Var)
#    microparticles[i].vy = np.random.normal(0,Var)
#    
#
#
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0.4,0,0),label='parabolic temperature gradient: begining')
#for i in range(500):
#    update(microparticles,temperature,density,kinetic_energy)
#calculate_temperature(temperature)
#plt.plot(location,temperature,color=(0.8,0,0),label='parabolic temperature gradient: after 500 time steps')




#plt.figure()
#plt.title('Finding suitable number of relaxing')
#plt.xlabel('number of steps')
#plt.ylabel('temperature variance')
#
#
#variance = []
#number = []
#
#for number_of_steps in range(0,10):
#    
#    ### initialization
#    
# temperature_gradient_flag = 1 
# temperature_initialization(temperature_gradient_flag,20,100)
#    
#        ###location
#        
# for i in range(number_of_microparticles): 
#    microparticles[i].m = micro_mass
#    switch = 0
#    while (switch==0): # I do not want the microparticles to be in same place at first step
#        microparticles[i].x = random()*L[0]
#        microparticles[i].y = random()*L[1]
#        switch = 1
#        for j in range(i):
#                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
#                    switch = 0        
#    for j in range(number_of_regions):
#        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
#            microparticles[i].region = j
#    
#    ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls
#
# for i in range(number_of_microparticles):
#    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
#    microparticles[i].vx = np.random.normal(0,Var)
#    microparticles[i].vy = np.random.normal(0,Var)
#    
# number.append(number_of_steps)
#
#
# for i in range(number_of_steps):
#    update(microparticles,temperature,density,kinetic_energy)
# calculate_temperature(temperature)
# variance.append(np.var(temperature_error))
# print(variance)
#
#
#plt.plot(number,variance,color = (0.5,0,0))




##plt.figure()
##plt.title('Temperature of 50 Regions for 5000 microparticles')
##plt.xlabel('x')
##plt.ylabel('temperature')
#




#plt.legend(loc=2, fontsize = 'x-small')
#plt.show()
           
            
            
                   


