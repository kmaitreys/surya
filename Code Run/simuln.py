import numpy as np
import matplotlib.pyplot as plt
import time
import os

plt.ion()

pi = np.pi

# plots final omega for r and theta
file_name = 'diffrot.dat'
fig1 = plt.figure()
ax1 = fig1.add_subplot(projection='polar')
ax1.grid(False)
#ax.set_yticks([0.55,0.7,0.85,1.0])
ax1.set_rorigin(0.25)
ax1.set_theta_zero_location('W', offset=90)
ax1.set_thetamin(0)
ax1.set_thetamax(180)

for i in range(0,10000):
    try :
        data = np.genfromtxt(file_name, usecols=np.arange(0,4))
        t = data[0][3]    # time
        print(t)
    except :
        time.sleep(0.0001)

     #azimuths = np.linspace(3.1415927, 0.0, 297)
    #zeniths = np.linspace(3.828, 6.96, 297)
    try :
        azimuths = np.array([data[i][0] for i in range(257)])
        zeniths = np.array([(data[i*257][1])/6.96 for i in range(257)])
        r, theta = np.meshgrid(zeniths, azimuths)
        values = np.zeros((257,257))
        for i in range(257):
            for j in range(257):
                values[i][j] = data[j*257+i][2]*10/(2*np.pi)
    except :
        time.sleep(0.0001)

    #print(r)
    #print(theta)
    # plotting is done here       
    C1 = ax1.contourf(theta, r, values,50,cmap=plt.cm.jet)
#    C2 = ax1.contour(theta,r, values,15,colors='black',linewidths=1)
    cb = plt.colorbar(C1,location='left',fraction=0.045, pad=0.0001)

    ax1.set_title(r'$\Omega/2\pi$ plot (t=%0.2f years)' %(t*3.17))

    fig1.canvas.draw()
    fig1.canvas.flush_events()

    cb.remove()
    print('done')