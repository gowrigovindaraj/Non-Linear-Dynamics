import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
plt.close('all')
 
eps = 0.3 #setting a value of epsilon
 
np.random.seed(2) #generate most random numbers
plt.figure(1)
for iter in range(0,50):
 
    j_final = np.pi*(1.5*np.random.random()-0.5)  #action statement
    theta_final = 2*np.pi*np.random.random()         #angle statement
 
    orbit = np.int(200*(j_final+np.pi/2))  #initialising values
    rplot = np.zeros(shape=(orbit,))
    thetaplot = np.zeros(shape=(orbit,))
    x = np.zeros(shape=(orbit,))
    y = np.zeros(shape=(orbit,))    
    for iter2 in range(0,orbit):   #iterating
        j_new = j_final + eps*np.sin(theta_final)
        theta_new = np.mod(theta_final+j_new,2*np.pi)
         
        rplot[iter2] = j_new
        thetaplot[iter2] = np.mod(theta_new-np.pi,2*np.pi) - np.pi            
           
        j_final = j_new
        theta_final = theta_new
         
        x[iter2] = (j_new+np.pi+0.25)*np.cos(theta_new)
        y[iter2] = (j_new+np.pi+0.25)*np.sin(theta_new)
         
    plt.plot(x,y,'o',ms=1)  #making a plot for the various iterations
    plt.title("KAM Twist Map for ε = 0.3")
    plt.xlabel("J")
    plt.ylabel("θ")
plt.show() 
