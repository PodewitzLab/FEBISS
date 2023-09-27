import math
import numpy as np
import timeit
import pandas as pd
import matplotlib.pyplot as plt
from numpy import array
start = timeit.default_timer() #according to: https://stackoverflow.com/a/16414144

Cl1=np.array([-1.678324450,-1.639160245,0.052476471])
C=np.array([-0.310918929,-0.497036348,0.251296752])-Cl1
Cl2=np.array([-0.833143382,0.890762121,1.258804266])-Cl1
Cl3=np.array([0.211305525,0.095815858,-1.357921146])-Cl1
Cl4=np.array([1.056486593,-1.335563125,1.051827416])-Cl1
Cl1=Cl1-Cl1

atom_list = [C,Cl1,Cl2,Cl3,Cl4] #TODO: Needs to be automated
label_list = ["C","Cl1","Cl2","Cl3","Cl4"] #TODO: Needs to be automated

print("Internal coordinates:")
for i in range(len(atom_list)):
    print("{0} = ".format(label_list[i]),"[",round(atom_list[i][0],4),round(atom_list[i][1],4),round(atom_list[i][2],4),"]")
print("")

def cos(ang):
    return np.cos(ang)

def sin(ang):
    return np.sin(ang)

def rotate(a,b,g,vec):
    m11 = cos(b)
    m21 = cos(a)*sin(b)
    m31 = sin(a)*sin(b)

    m12 = -cos(g)*sin(b)
    m22 = cos(a)*cos(b)*cos(g)-sin(a)*sin(g)
    m32 = cos(a)*sin(g)+cos(b)*cos(g)*sin(a)

    m13 = sin(b)*sin(g)
    m23 = -cos(g)*sin(a)-cos(a)*cos(b)*sin(g)
    m33 = cos(a)*cos(g)-cos(b)*sin(a)*sin(g)

    r = np.array([m11,m21,m31])*vec[0]
    s = np.array([m12,m22,m32])*vec[1]
    t = np.array([m13,m23,m33])*vec[2]
    rvec = r+s+t

    return rvec

def angle(vec1,vec2):
    norm1 = ((vec1**2).sum()**0.5)
    norm2 = ((vec2**2).sum()**0.5)
    vec1 = vec1/norm1 #normalization
    vec2 = vec2/norm2 #normalization
    angle = np.dot(vec1,vec2)/(norm1*norm2)

    return angle

def dist(vec1:np.array,vec2:np.array) -> int:
    diff = vec1-vec2
    if vec1.size == 3:
        dist = np.sqrt(diff[0]**2+diff[1]**2+diff[2]**2)
    elif vec1.size == 4:
        dist = np.sqrt(diff[0]**2+diff[1]**2+diff[2]**2+diff[3]**2)
    return dist

def quat(a,b,g): #from Euler: https://de.wikipedia.org/wiki/Quaternion#Bezug_zu_Eulerwinkeln
    w = cos(a/2)*cos(b/2)*cos(g/2)-sin(a)*cos(b)*sin(g)
    x = cos(a)*sin(b)*cos(g)+sin(a)*sin(b)*sin(g)
    y = -cos(a)*sin(b)*sin(g)+sin(a)*sin(b)*cos(g)
    z = sin(a)*cos(b)*cos(g)+cos(a)*cos(b)*sin(g)
    return np.array([w,x,y,z])


step = 0.05
threshold= 0.2
a = -np.pi
b = step #gimbal lock prevention
g = -np.pi
quat_list = []
count = 0
while a >= -np.pi and a < np.pi:
    while b > 0 and b < np.pi:
        while g >= -np.pi and g < np.pi:
            overlap1 = False
            overlap2 = False
            for i in atom_list:
                rvec = rotate(a,b,g,i)
                if overlap1 == False:
                    for j in atom_list:
                        if dist(rvec,j) <= 0.5:
                            overlap1 = True
                            break
                    if overlap1 == False:
                        break
                else:
                    for k in atom_list:
                        if dist(rvec,k) <= 0.5:
                            overlap2 = True
                            break
                    if overlap2 == False:
                        break
            if overlap1 == True and overlap2 == True:
                min_dist = 0
                for i in range(len(atom_list)):
                    rvec = rotate(a,b,g,atom_list[i])
                    min_dist += np.sqrt(min([
                        (i[0]-rvec[0])**2+(i[1]-rvec[1])**2+(i[2]-rvec[2])**2 for i in atom_list]))
                min_dist /= len(atom_list)
                if min_dist <= threshold and b != 0: #gimbal lock
                    quat_list.append(quat(a,b,g))
                    count += 1
            g += step
        g = -np.pi
        b += step
    b = step
    a += step
print("{0} representations found!\n".format(count))
dist_list = []
for i in range(len(quat_list)):
    for j in range(i,len(quat_list)):
        dist_list.append(dist(quat_list[i],quat_list[j]))
dist_list = np.array(dist_list)
print("max = ",max(dist_list))

df = pd.DataFrame(dist_list)
#df.hist()
#plt.show()
quat_list_filtered = []
for i in range(len(quat_list)):
    for j in range(i+1,len(quat_list)):
        if dist(quat_list[i],quat_list[j]) < 1.5:
            break
        quat_list_filtered.append(quat_list[i])

print("{0} quaternions remaining!".format(len(quat_list_filtered)))
stop = timeit.default_timer()
#print("{0} equivalent structures found in:".format(count), round(stop - start,3), "s")

# print("-----------(^_^)--------------\n"
#                           "Equivalent structure found. \n"
#                           "a = {0} rad\n"
#                           "b = {1} rad\n"
#                           "g = {2} rad\n\n"
#                           "The new coordinates are:".format(round(a,4), round(b,4), round(g,4)))
# for i in range(len(atom_list)):
#     print("{0} = ".format(label_list[i]),"[",round(rot_list[i][0],4),round(rot_list[i][1],4),round(rot_list[i][2],4),"]")
# print("-----------(^_^)--------------\n")