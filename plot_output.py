import numpy as np
import matplotlib.pyplot as plt

f = open("output_lc_np_2.dat",'r')

lines = f.readlines()
lines = lines[12:]
data = np.ones([len(lines),2])

for i in range(len(lines)):
    l = lines[i]
    data[i,:] = np.float128(l.split())

initial = 0
final = -1

x = data[initial:final,0]
y = data[initial:final,1]

#x = np.log10(data[:,0])
#y = np.log10(data[:,1])

#plt.scatter(x,y)
plt.loglog(x,y,label="-np 2")


f = open("output_lc_np_4.dat",'r')

lines = f.readlines()
lines = lines[12:]
data = np.ones([len(lines),2])

for i in range(len(lines)):
    l = lines[i]
    data[i,:] = np.float128(l.split())

initial = 0
final = -1

x = data[initial:final,0]
y = data[initial:final,1]

#x = np.log10(data[:,0])
#y = np.log10(data[:,1])

#plt.scatter(x,y)
plt.loglog(x,y,label = "-np 4")


f = open("output_lc_np_6.dat",'r')

lines = f.readlines()
lines = lines[12:]
data = np.ones([len(lines),2])

for i in range(len(lines)):
    l = lines[i]
    data[i,:] = np.float128(l.split())

initial = 0
final = -1

x = data[initial:final,0]
y = data[initial:final,1]

#x = np.log10(data[:,0])
#y = np.log10(data[:,1])

#plt.scatter(x,y)
plt.loglog(x,y,label = "-np 6")


f = open("output_lc_np_6_2.dat",'r')

lines = f.readlines()
lines = lines[12:]
data = np.ones([len(lines),2])

for i in range(len(lines)):
    l = lines[i]
    data[i,:] = np.float128(l.split())

initial = 0
final = -1

x = data[initial:final,0]
y = data[initial:final,1]

#x = np.log10(data[:,0])
#y = np.log10(data[:,1])

#plt.scatter(x,y)
plt.loglog(x,y,label = "-np 6_2")




plt.legend()
plt.xlabel("Time (days)")
plt.ylabel(r"Flux (mJy/$m^2$)")
plt.title("Afterglow Gaussian Jet")
plt.suptitle("On-axis JET code")
#plt.savefig("plots/afterglow_onaxis_gaussian_jet.png", dpi = 420)
plt.show()
