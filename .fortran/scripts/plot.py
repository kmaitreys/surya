import numpy as np
import matplotlib.pyplot as plt

pi = np.pi

# plots final omega for r and theta
file_name = 'diffrot.dat'

data = np.genfromtxt(file_name)

 #azimuths = np.linspace(3.1415927, 0.0, 297)
#zeniths = np.linspace(3.828, 6.96, 297)

azimuths = np.array([data[i][0] for i in range(257)])
zeniths = np.array([(data[i*257][1])/6.96 for i in range(257)])
r, theta = np.meshgrid(zeniths, azimuths)
values = np.zeros((257,257))
for i in range(257):
    for j in range(257):
        values[i][j] = data[j*257+i][2]*10/(2*np.pi)

#print(r)
#print(theta)
# plotting is done here    
fig1 = plt.figure()
ax1 = fig1.add_subplot(projection='polar')    
C1 = ax1.contourf(theta, r, values,100,cmap=plt.cm.plasma)
C2 = ax1.contour(theta,r, values,15,colors='black',linewidths=1)
ax1.grid(False)
ax1.set_title(r'$\Omega/2\pi$ plot')
#ax.set_yticks([0.55,0.7,0.85,1.0])
ax1.set_rorigin(0.25)
ax1.set_theta_zero_location('W', offset=90)
ax1.set_thetamin(0)
ax1.set_thetamax(180)
plt.colorbar(C1, ax=[ax1], location='left',pad=0.0001)

plt.savefig('omega.jpg',dpi=500, bbox_inches='tight')
plt.show()

################################################
################################################
'''
# plots omega at fixed r for different t and theta
file_name = 'omg_lat.dat'
data = np.genfromtxt(file_name)

num=700
jmax = 257
kmax = int(8/(0.0002*40))
print(kmax)

azimuth = np.array([data[j][0]*180/pi for j in range(jmax)])
time = np.array([(data[k*jmax][1]*3.169) for k in range(kmax-num,kmax)])

t, theta = np.meshgrid(time, azimuth)
vals = np.zeros((jmax,num))
for j in range(jmax):
    for k in range(kmax-num,kmax):
        vals[j][k-kmax] = data[k*jmax+j][2]*10/(2*np.pi)
    vals[j] = (vals[j]-np.mean(vals[j]))

#print(theta)
# plotting is done here        

fig2 = plt.figure(figsize=(20,10))
ax2 = fig2.add_subplot(111)
C1 = ax2.contourf(t, theta-90, vals,100,cmap=plt.cm.plasma)
C2 = ax2.contour(t,theta-90, vals,15,colors='black',linewidths=1)
ax2.grid(False)
ax2.set_title(r'$\Delta(\Omega/2\pi)$ plot')
#ax2.set_yticks([0.55,0.7,0.85,1.0])
#ax2.set_rorigin(0.25)
#ax2.set_theta_zero_location('W', offset=90)
#ax2.set_thetamin(0)
#ax2.set_thetamax(180)
#ax2.set_xlim(60,90)
ax2.set_xlabel('time (years)')
ax2.set_ylabel('latitude (deg)')
plt.colorbar(C1, ax=[ax2], location='right',pad=0.02)

plt.savefig('omega_time_lat.jpg',dpi=500, bbox_inches='tight')
plt.show()
'''