import numpy as np

def normalize(v1: np.ndarray):
    return v1 / np.linalg.norm(v1)

def gigist_quat(X: np.ndarray, V2: np.ndarray):
    """
    way of calculating the quaternion representation of the orientation of the solvent according to
    https://github.com/maberl1/gigist/blob/9be781be4099f559a2413a553d3464f48ee32a57/Quaternion.h
    """

    X = normalize(X)
    Z = np.cross(X, V2)
    Z = normalize(Z)
    Y = np.cross(Z, X)
    Y = normalize(Y)

    m11 = X[0]
    m12 = Y[0]
    m13 = Z[0]

    m21 = X[1]
    m22 = Y[1]
    m23 = Z[1]

    m31 = X[2]
    m32 = Y[2]
    m33 = Z[2]

    trace = m11 + m22 + m33
    s = 0

    if trace > 0:
        s = 0.5 / np.sqrt(trace + 1)
        w = 0.25 / s
        x = (m32 - m23) * s
        y = (m13 - m31) * s
        z = (m21 - m12) * s

    elif (m11 > m22 and m11 > m33):
        s = 2 * np.sqrt(1.0 + m11 - m22 - m33)
        w = (m32 - m23) / s
        x = 0.25 * s
        y = (m12 + m21) / s
        z = (m13 + m31) / s

    elif m22 > m33:
        s = 2.0 * np.sqrt(1.0 + m22 - m11 - m33)
        w = (m13 - m31) / s
        x = (m12 + m21) / s
        y = 0.25 * s
        z = (m23 + m32) / s

    else:
        s = 2.0 * np.sqrt(1.0 + m33 - m11 - m22)
        w = (m21 - m12) / s
        x = (m13 + m31) / s
        y = (m23 + m32) / s
        z = 0.25 * s

    return w, x, y, z


def gist_quat(X: np.ndarray, V2: np.ndarray):
    """
    way of calculating the quaternion representation of the orientation of the solvent according to
    https://github.com/Amber-MD/cpptraj/blob/21def2c3a5b9c81c5a1fe68498c0ec19dba61574/src/Action_GIST.cpp
    """

    x_lab_ = np.array([1.0, 0.0, 0.0])
    y_lab_ = np.array([0.0, 1.0, 0.0])
    z_lab_ = np.array([0.0, 0.0, 1.0])

    X = normalize(X)
    V2 = normalize(V2)

    ar1 = np.cross(X, x_lab_)
    sar = ar1
    ar1 = normalize(ar1)

    dp1 = np.dot(x_lab_, X)
    theta = np.arccos(dp1)
    sign = np.dot(sar, X)

    if sign > 1e-14:
        theta = theta / 2.0
    else:
        theta = theta / -2.0

    w1 = np.cos(theta)
    sin_theta = np.sin(theta)
    x1 = ar1[0] * sin_theta
    y1 = ar1[1] * sin_theta
    z1 = ar1[2] * sin_theta
    w2 = w1
    x2 = x1
    y2 = y1
    z2 = z1

    temp_0 = ((w2 * w2 + x2 * x2) - (y2 * y2 + z2 * z2)) * X[0]
    temp_0 = (2 * (x2 * y2 + w2 * z2) * X[1]) + temp_0
    temp_0 = (2 * (x2 * z2 - w2 * y2) * X[2]) + temp_0

    temp_1 = 2 * (x2 * y2 - w2 * z2) * X[0]
    temp_1 = ((w2 * w2 - x2 * x2 + y2 * y2 - z2 * z2) * X[1]) + temp_1
    temp_1 = (2 * (y2 * z2 + w2 * x2) * X[2]) + temp_1

    temp_2 = 2 * (x2 * z2 + w2 * y2) * X[0]
    temp_2 = (2 * (y2 * z2 - w2 * x2) * X[1]) + temp_2
    temp_2 = ((w2 * w2 - x2 * x2 - y2 * y2 + z2 * z2) * X[2]) + temp_2

    temp = np.array([temp_0, temp_1, temp_2])

    X = temp

    temp2_0 = ((w2 * w2 + x2 * x2) - (y2 * y2 + z2 * z2)) * V2[0]
    temp2_0 = (2 * (x2 * y2 + w2 * z2) * V2[1]) + temp2_0
    temp2_0 = (2 * (x2 * z2 - w2 * y2) * V2[2]) + temp2_0

    temp2_1 = 2 * (x2 * y2 - w2 * z2) * V2[0]
    temp2_1 = ((w2 * w2 - x2 * x2 + y2 * y2 - z2 * z2) * V2[1]) + temp2_1
    temp2_1 = (2 * (y2 * z2 + w2 * x2) * V2[2]) + temp2_1

    temp2_2 = 2 * (x2 * z2 + w2 * y2) * V2[0]
    temp2_2 = (2 * (y2 * z2 - w2 * x2) * V2[1]) + temp2_2
    temp2_2 = ((w2 * w2 - x2 * x2 - y2 * y2 + z2 * z2) * V2[2]) + temp2_2

    temp2 = np.array([temp2_0, temp2_1, temp2_2])

    V2 = temp2

    ar2 = np.cross(temp, temp2)
    ar2 = normalize(ar2)
    dp2 = np.dot(ar2, z_lab_)
    theta = np.arccos(dp2)

    sar = np.cross(ar2, z_lab_)
    sign = np.dot(sar, temp)

    if sign < 0:
        theta /= 2.0
    else:
        theta /= -2.0

    w3 = np.cos(theta)
    sin_theta = np.sin(theta)
    x3 = x_lab_[0] * sin_theta

    w4 = w1 * w3 - x1 * x3
    x4 = w1 * x3 + x1 * w3
    y4 = y1 * w3 + z1 * x3
    z4 = -y1 * x3 + z1 * w3

    return w4, x4, y4, z4