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
from scipy import stats




### constants and initializations which we may want to change in differen experiments



    
    ### more general constants

k = 1 
L = [10.0,10.0] #size of the room
number_of_microparticles = 500
number_of_regions = 10
deltaX = L[0]/number_of_regions
deltaT = 0.001
number_of_time_steps = 500
show_step = 5

    ### flags

thermostat_flag = 1 # on (1) or off (0)
temperature_gradient_flag = 1  # temperature gradient flag: 0: uniform ; 1: linear ; 2: parabolic
surface_temperature_flag = 0 # 0: no specific temperature for macroparticle  1: each side of macroparticle has a temperature 
microparticles_show_flag = 1 # on (1) or off (0)

    ###constants related to the classes

micro_mass = 1
macro_mass = 1
radius_of_the_macroparticle = 0.5

    ### temperature

location = np.zeros(number_of_regions) #location of the reigon along x axis     
Temperature = np.zeros(number_of_regions) #temperature of the reigon along x axis
min_temperature = 10
max_temperature = 100
left_surface_temperature = 20
right_surface_temperature = 20

def temperature_initialization(temperature_gradient_flag,min_temperature,max_temperature):
    for i in range(int(number_of_regions)):
        location[i] = i*deltaX +deltaX/2
    ### temperature gradient flag:
        ### 0 : uniform
        ### 1 : linear
        ### 2 : styloid
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
                v2 = microparticles[j].vx**2 + microparticles[j].vy**2
                kinetic_energy[i] += 0.5*microparticles[j].m*v2
                density[i] += 1
                microparticles[j].region = i
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
                microparticles[j].region = i
                
    
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
    momentumX += particle.m*particle.vx
    return momentumX


def calculate_momentumY():
    
    momentumY = 0.0
    for i in range( int(number_of_microparticles) ): 
        momentumY += microparticles[i].m*microparticles[i].vy
    momentumY += particle.m*particle.vy
    return momentumY




    ### update


    
    
def update(microparticles,particle,temperature,density,kinetic_energy):       
    for i in range(number_of_microparticles):    
        microparticles[i].x += deltaT*microparticles[i].vx
        microparticles[i].y += deltaT*microparticles[i].vy
        particle.vx,particle.vy = microparticles[i].check_event(particle)
    if(thermostat_flag==1):
        thermostat(microparticles,temperature,density,kinetic_energy)
    particle.x += deltaT*particle.vx
    particle.y += deltaT*particle.vy
    particle.check_event()




### classes
   
    
    

class Macroparticle:
    
    switch = 1
    
    x = 0.0 #location of each macroparticle
    y = 0.0 #location of each macroparticle
    vx = 0.0 #velocity of each macroparticle
    vy = 0.0 #velocity of each macroparticle
    fx = 0.0 #force of each macroparticle
    fy = 0.0 #force of each macroparticle
    m = 0.0 #mass
    r = 0.0 #radius
    left_temperature = 0.0
    right_temperature = 0.0    
 
    def check_event(self):
        if(self.x<0):
            self.x = L[0]
            self.switch = 0
        if(self.x>L[0]):
            self.x = 0
            self.switch = 0
        if(self.y<0):
            self.y = L[1]
            self.switch = 0
        if(self.y>L[1]):
            self.y = 0
            self.switch = 0
            
            
            
            
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
 
    def macroparticle_collision(self,macroparticle):
        if(surface_temperature_flag==1):
#            if(self.x-macroparticle.x<0):
#                temp = macroparticle.left_temperature
#            else:
#                temp = macroparticle.right_temperature
#            self.vx = np.random.normal(0,k*temp/self.m)
#            self.vy = np.random.normal(0,k*temp/self.m)        
#            macroparticle.vx = -1*self.m*self.vx/(macroparticle.m)
#            macroparticle.vy = -1*self.m*self.vy/(macroparticle.m)
            px = self.m*self.vx + macroparticle.m*macroparticle.vx
            py = self.m*self.vy + macroparticle.m*macroparticle.vy
            if(self.x-macroparticle.x<0):
                temp = macroparticle.left_temperature
            else:
                temp = macroparticle.right_temperature
            self.vx = np.random.normal(0,k*temp/self.m)
            self.vy = np.random.normal(0,k*temp/self.m)        
            macroparticle.vx =  (px-self.m*self.vx)/(macroparticle.m)
            macroparticle.vy =  (py-self.m*self.vy)/(macroparticle.m)            
        else:  
            teta = (self.vx-macroparticle.vx)*(self.x-macroparticle.x) + (self.vy-macroparticle.vy)*(self.y-macroparticle.y)
            teta = teta / ( (self.x-macroparticle.x)**2 + (self.y-macroparticle.y)**2 )
            self.vx = self.vx -2*macroparticle.m/(macroparticle.m+self.m)*teta*(self.x-macroparticle.x)
            self.vy = self.vy -2*macroparticle.m/(macroparticle.m+self.m)*teta*(self.y-macroparticle.y)
            macroparticle.vx = macroparticle.vx -2*self.m/(macroparticle.m+self.m)*teta*(macroparticle.x-self.x)
            macroparticle.vy = macroparticle.vy -2*self.m/(macroparticle.m+self.m)*teta*(macroparticle.y-self.y)
        return macroparticle.vx,macroparticle.vy
    
    def check_event(self,macroparticle):    
        VX = macroparticle.vx
        VY = macroparticle.vy
        if(self.x<0):
            self.x = L[0]
        if(self.x>L[0]):
            self.x = 0
        if(self.y<0):
            self.y = L[1]
        if(self.y>L[1]):
            self.y = 0                    
        if( ((self.x-macroparticle.x)**2+(self.y-macroparticle.y)**2) <= (macroparticle.r+0.05)**2):
            VX,VY = self.macroparticle_collision(macroparticle)
        return VX,VY



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
        
    ### mobile objects

particle = Macroparticle()            
microparticles = [Microparticle() for i in range(number_of_microparticles)]




### main


    ### initialization
    
temperature_initialization(temperature_gradient_flag,min_temperature,max_temperature)
particle.left_temperature = left_surface_temperature
particle.right_temperature = right_surface_temperature

        ### macroparticle

particle.x = L[0]/2
particle.y = L[1]/2
particle.vx = 0.0
particle.vy = 0.0
particle.m = macro_mass
particle.r = radius_of_the_macroparticle

        ###microparticles
        
            ###location
        
for i in range(number_of_microparticles): 
    microparticles[i].m = micro_mass
    switch = 0
    while (switch==0): # I do not want the microparticles to be in same place at first step
        microparticles[i].x = random()*L[0]
        microparticles[i].y = random()*L[1]
        switch = 1
        if( ((microparticles[i].x-particle.x)**2+(microparticles[i].y-particle.y)**2) <= (particle.r+0.05)**2 ):
            switch = 0
        else:
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
thermostat(microparticles,temperature,density,kinetic_energy)

fig = plt.figure()
for i in range(number_of_time_steps):
    if(particle.switch==1): # The program will be running until the particle will arive at sides
        if(i%show_step==0):  
            if(microparticles_show_flag==1):
                for i in range(number_of_microparticles):
                    color_code = Temperature[int(microparticles[i].region)]/max(Temperature)
                    microparticle_show = plt.Circle((microparticles[i].x,microparticles[i].y),0.05,color=(color_code, 0.0, 0.0))
                    fig = plt.gcf()
                    f = fig.gca()
                    f.add_artist(microparticle_show)
                particle_show = plt.Circle((particle.x,particle.y),particle.r,color=(0.0, 0.0, 0.5))
                fig = plt.gcf()
                f = fig.gca()
                f.add_artist(particle_show)
                plt.pause(0.0000001)
                fig.clf()
                #plt.plot(trackX,trackY,'b.')
                plt.xlim(0, L[0])
                plt.ylim(0, L[1])
                plt.gca().set_aspect('equal', adjustable='box')
        update(microparticles,particle,temperature,density,kinetic_energy)
#        if(i%10==0):
#            for c in range(number_of_microparticles):
#                phi = random()*2*np.pi
#                vxNew = microparticles[c].vx*np.cos(phi) - microparticles[c].vy*np.sin(phi)
#                vyNew = microparticles[c].vx*np.sin(phi) + microparticles[c].vy*np.cos(phi)
#                microparticles[c].vx = vxNew
#                microparticles[c].vy = vyNew
#        trackX.append(particle.x)
#        trackY.append(particle.y)  
fig.show()

#fig2 = plt.figure()
#plt.plot(trackX,trackY)
#fig2.show()







#fig = plt.figure()
#time = []
#sigma2 = []
#plt.xlabel('time')
#plt.ylabel('Mean squared displacement')
#
#for number_of_time_steps in range(100,5000,100):
# x2 = 0
# t = 0
#    ### initialization
#    
# temperature_initialization(temperature_gradient_flag,min_temperature,max_temperature)
# particle.left_temperature = left_surface_temperature
# particle.right_temperature = right_surface_temperature
#
#        ### macroparticle
#
# particle.x = L[0]/2
# particle.y = L[1]/2
# particle.vx = 0.0
# particle.vy = 0.0
# particle.m = macro_mass
# particle.r = radius_of_the_macroparticle
#
#        ###microparticles
#        
#            ###location
#        
# for i in range(number_of_microparticles): 
#    microparticles[i].m = micro_mass
#    switch = 0
#    while (switch==0): # I do not want the microparticles to be in same place at first step
#        microparticles[i].x = random()*L[0]
#        microparticles[i].y = random()*L[1]
#        switch = 1
#        if( ((microparticles[i].x-particle.x)**2+(microparticles[i].y-particle.y)**2) <= (particle.r+0.05)**2 ):
#            switch = 0
#        else:
#            for j in range(i):
#                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
#                    switch = 0        
#    for j in range(number_of_regions):
#        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
#            microparticles[i].region = j
#    
#        ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls
#
# for i in range(number_of_microparticles):
#    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
#    microparticles[i].vx = np.random.normal(0,Var)
#    microparticles[i].vy = np.random.normal(0,Var)
# thermostat(microparticles,temperature,density,kinetic_energy)
#
#
#
##calculate_temperature(temperature)
##plt.plot(location,Temperature,'g')
##plt.plot(location,temperature,'r')
# for i in range(number_of_time_steps):
#  if(particle.switch==1):
#    particle_show = plt.Circle((particle.x,particle.y),particle.r,color=(0.0, 0.0, 0.5))
#    fig = plt.gcf()
#    f = fig.gca()
#    f.add_artist(particle_show)
#    plt.pause(0.000001)
#    fig.clf()
#    plt.xlim(0, L[0])
#    plt.ylim(0, L[1])
#    plt.gca().set_aspect('equal', adjustable='box')
#    update(microparticles,particle,temperature,density,kinetic_energy)
##    if(i%10==0):
##        for c in range(number_of_microparticles):
##            phi = random()*2*np.pi
##            vxNew = microparticles[c].vx*np.cos(phi) - microparticles[c].vy*np.sin(phi)
##            vyNew = microparticles[c].vx*np.sin(phi) + microparticles[c].vy*np.cos(phi)
##            microparticles[c].vx = vxNew
##            microparticles[c].vy = vyNew
##    x2 += (particle.x-L[0]/2)**2
##    t +=1
##for i in range(49,N,50):
##    s = 0
##    for j in range(N-i-1):
##        s = s+(x[j+i]-x[j])**2
##    s = s/(N-i)
##    sigma2.append(s)
##    time.append(i*deltaT)
## x2 /= t    
## sigma2.append(x2)
## time.append(t*deltaT)
#
##plt.plot(time,sigma2,'b.')
##slope, intercept, r_value, p_value, std_err = stats.linregress(time, sigma2)
##print("slope",slope)
#plt.show()
                  




#fig = plt.figure()
#time = []
#sigma2 = []
#x = []
#plt.xlabel('time')
#plt.ylabel('Mean squared displacement')
#
#number_of_time_steps=1000
#
#    ### initialization
#    
#temperature_initialization(temperature_gradient_flag,min_temperature,max_temperature)
#particle.left_temperature = left_surface_temperature
#particle.right_temperature = right_surface_temperature
#
#        ### macroparticle
#
#particle.x = L[0]/2
#particle.y = L[1]/2
#particle.vx = 0.0
#particle.vy = 0.0
#particle.m = macro_mass
#particle.r = radius_of_the_macroparticle
#
#        ###microparticles
#        
#            ###location
#        
#for i in range(number_of_microparticles): 
#    microparticles[i].m = micro_mass
#    switch = 0
#    while (switch==0): # I do not want the microparticles to be in same place at first step
#        microparticles[i].x = random()*L[0]
#        microparticles[i].y = random()*L[1]
#        switch = 1
#        if( ((microparticles[i].x-particle.x)**2+(microparticles[i].y-particle.y)**2) <= (particle.r+0.05)**2 ):
#            switch = 0
#        else:
#            for j in range(i):
#                if( microparticles[j].x==microparticles[i].x and microparticles[j].y==microparticles[i].y ):
#                    switch = 0        
#    for j in range(number_of_regions):
#        if(microparticles[i].x<=location[j]+deltaX/2 and microparticles[i].x>location[j]-deltaX/2):
#            microparticles[i].region = j
#    
#        ### velocity initialzation based based on Boltzman distribution and linear temperature gradiant between two walls
#
#for i in range(number_of_microparticles):
#    Var = np.sqrt(k*Temperature[microparticles[i].region]/microparticles[i].m)
#    microparticles[i].vx = np.random.normal(0,Var)
#    microparticles[i].vy = np.random.normal(0,Var)
#thermostat(microparticles,temperature,density,kinetic_energy)
#
#
#
##calculate_temperature(temperature)
##plt.plot(location,Temperature,'g')
##plt.plot(location,temperature,'r')
#for i in range(number_of_time_steps):
#  if(particle.switch==1):
##    particle_show = plt.Circle((particle.x,particle.y),particle.r,color=(0.0, 0.0, 0.5))
##    fig = plt.gcf()
##    f = fig.gca()
##    f.add_artist(particle_show)
##    plt.pause(0.000001)
##    fig.clf()
##    plt.xlim(0, L[0])
##    plt.ylim(0, L[1])
##    plt.gca().set_aspect('equal', adjustable='box')
#    update(microparticles,particle,temperature,density,kinetic_energy)
#    if(i%1==0):
#        for c in range(number_of_microparticles):
#            phi = random()*2*np.pi
#            vxNew = microparticles[c].vx*np.cos(phi) - microparticles[c].vy*np.sin(phi)
#            vyNew = microparticles[c].vx*np.sin(phi) + microparticles[c].vy*np.cos(phi)
#            microparticles[c].vx = vxNew
#            microparticles[c].vy = vyNew
#    x.append(particle.x-L[0]/2)
#N = len(x)
#for i in range(49,N,50):
#    s = 0
#    for j in range(N-i-1):
#        s = s+(x[j+i]-x[j])**2
#    s = s/(N-i)
#    sigma2.append(s)
#    time.append(i*deltaT)
#
#plt.plot(time,sigma2,'b.')
#slope, intercept, r_value, p_value, std_err = stats.linregress(time, sigma2)
#print("slope",slope)
#plt.show()