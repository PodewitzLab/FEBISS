import math
import numpy as np
from numpy import array
O_T = array([0,0,0])
H1_T = array([1,0,0])
H2_T = array([0,math.sqrt(0.5),math.sqrt(0.5)])

def norm(vec:np.ndarray) -> np.ndarray:
    vec = vec/((vec**2).sum()**0.5)
    return vec


def quat(H1, H2):
    w = 0
    x = 0
    y = 0
    z = 0
    X = norm(H1)
    print("X = {0}".format(X))
    Z = norm(np.cross(X, H2))
    print("Z = {0}".format(Z))
    Y = norm(np.cross(Z, X))
    print("Y= {0}".format(Y))

    m11 = X[0]
    m21 = X[1]
    m31 = X[2]

    m12 = Y[0]
    m22 = Y[1]
    m32 = Y[2]

    m13 = Z[0]
    m23 = Z[1]
    m33 = Z[2]

    trace = m11 + m22 + m33
    s = 0

    if trace > 0:
        s = 0.5 / math.sqrt(trace + 1)
        x = (m32 - m23) * s
        y = (m13 - m32) * s
        z = (m21 - m12) * s
    elif m11 > m22 and m11 > m33:
        s = 2 * math.sqrt(1 + m11 - m22 - m33)
        w = (m32 - m23) / s
        x = 0.25 * s
        y = (m12 + m21) / s
        z = (m13 + m31) / s
    elif m22 > m33:
        s = 2.0 * math.sqrt(1 + m22 - m11 - m33)
        w = (m13 - m31) / s
        x = (m12 + m21) / s
        y = 0.25 * s
        z = (m23 + m32) / s
    else:
        s = 2.0 * math.sqrt(1.0 + m33 - m11 - m22)
        w = (m21 - m12) / s
        x = (m13 + m31) / s
        y = (m23 + m32) / s
        z = 0.25 * s

    return np.array([w,x,y,z])

def multiply(q1,q2):
    w = q1[0]
    x = q1[1]
    y = q1[2]
    z = q1[3]
    w = w * q2[0] - x * q2[1] - y * q2[2] - z * q2[3]
    x = w * q2[1] + x * q2[0] + y * q2[3] - z * q2[2]
    y = w * q2[2] - x * q2[3] + y * q2[0] + z * q2[1]
    z = w * q2[3] + x * q2[2] - y * q2[1] + z * q2[0]
    return np.array([w,x,y,z])

def invert(q1):
    w = q1[0]
    x = -q1[1]
    y = -q1[2]
    z = -q1[3]
    return np.array([w,x,y,z])

def rotate(q1,vec):
    w = 0
    x = vec[0]
    y = vec[1]
    z = vec[2]
    vecQuat = np.array([w,x,y,z])
    vecQuat = multiply(q1,vecQuat)
    vecQuat = multiply(vecQuat,invert(q1))
    return np.array([vecQuat[1],vecQuat[2],vecQuat[3]])

with open("coords.txt", "r") as c:
    text = c.readlines()
    i = 0
    while i < len(text):
        #print(line.split())
        split = text[i].split()
        if split[-1] == "O":
            H1_line = text[i+1].split()
            H2_line = text[i+2].split()
            O = np.array([float(split[0]),float(split[1]),float(split[2])])
            H1 = np.subtract(np.array([float(H1_line[0]),float(H1_line[1]),float(H1_line[2])]),O)
            H2 = np.subtract(np.array([float(H2_line[0]),float(H2_line[1]),float(H2_line[2])]),O)
            O = np.subtract(O,O)
            print("O = {0}".format(O))
            print("H1 = {0}".format(H1))
            print("H2 = {0}".format(H2))
            q = quat(H1,H2)
            print("q = {0}".format(q))
            H1_new = norm(rotate(q,H1_T))
            print("H1_new = {0}".format(H1_new))
            print("---\n")
            i+=3
        else:
            continue
