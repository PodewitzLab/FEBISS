import numpy as np
from pymatgen.core import Molecule
import quaternion as quat
from .averageQuaternions import averageQuaternions as avg

def normalize(v1: np.ndarray):
    return v1 / np.linalg.norm(v1)

def inv(q:quat.quaternion):
    """
    Returns q^(-1)

    :param q: quaternion object
    :return: inverse of q
    """
    return quat.quaternion.inverse(q)

def norm(q:quat.quaternion):
    """
    Returns ||q||

    :param q: quaternion object
    :return: norm of q
    """
    return quat.quaternion.norm(q)

def unity(q:quat.quaternion):
    """
    Sets ||q|| -> 1

    :param q: quaternion object
    :return: normalized q
    """
    return q/norm(q)

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

def vec_to_quat(v) -> quat.quaternion:
    """
    :param v:
    :return:

    Takes 3D-coordinates v and returns quaternion [0,v]
    """
    help_list = [0]
    help_list.extend(v)
    return quat.from_float_array(np.array(help_list))

def vec_part(q:quat.quaternion):
    help_list = quat.as_float_array(q)
    vec = list(help_list[1:4])
    return vec

def calc_quats(mol: Molecule, i1: int, i2: int):
    """
    Calculates characteristic quaternion from positions of center of mass and of atoms i1 and i2 the way the quaternion
    is calculated in cpptraj GIST.

    :param mol: Molecule object
    :param i1: index number of atom 1
    :param i2: index number of atom 2
    :return: quaternion object
    """
    at1 = (mol.cart_coords[i1] - mol.center_of_mass) / np.linalg.norm(
        mol.cart_coords[i1] - mol.center_of_mass)
    at2 = (mol.cart_coords[i2] - mol.center_of_mass) / np.linalg.norm(
        mol.cart_coords[i2] - mol.center_of_mass)
    q = quat.from_float_array(np.array(gist_quat(at1, at2)))
    return q


def distance(q1:quat.quaternion, q2:quat.quaternion):
    """
    Calculates the distance between two quaternions using the formula theta = 2 * arccos (|q1*q2|) where |q1*q2| denotes
    the absolute value of the dot product between two quaternions.

    This metric is used also for S(orient) calculation in GIST:
    Ramsey, S., et al. Journal of Computational Chemistry 2016, 37 (21), 2029–2037. https://doi.org/10.1002/jcc.24417.

    For formula see for example:
    p. 273 in Hanson, A. J. Visualizing Quaternions; Morgan Kaufmann Publishers Inc.: San Francisco, CA, USA, 2006.

    :param q1: quaternion object
    :param q2: quaternion object
    :return: theta (scalar)
    """

    q1_array = quat.as_float_array(q1)
    q2_array = quat.as_float_array(q2)
    scalar_product = np.abs(q1_array[0]*q2_array[0]+
                            q1_array[1]*q2_array[1]+
                            q1_array[2]*q2_array[2]+
                            q1_array[3]*q2_array[3])
    theta = 2*np.arccos(scalar_product)
    return theta


def new_coord_gen(mol:Molecule, q:quat.quaternion,t: np.ndarray = np.array([0,0,0])) -> list[list[float]]:
    """
    Calculates the new coordinates after rotation of mol by q using the formula:
    [0,x'] = q * [0,x] * q^(-1)

    NEW LM20231124: Created possibility for translation of molecules after rotation with vector t. For the addition of
    translation vector the vec_part first has to be transformed into np.ndarray. For usage with structure_randomizer.write(),
    the np.ndarray is transformed back to a list (probably not needed, though).

    :param mol: Molecule object
    :param q: quaternion object
    :param t: translation vector as np.ndarray
    :return: new coords as list of arrays
    """
    q = unity(q)
    new_coords = []
    for m in range(len(mol.cart_coords)):
        q_coord = unity(vec_to_quat(mol.cart_coords[m]))
        q_coord_trans = unity(q * q_coord * inv(q))
        new_coords.append(list(np.array(vec_part(q_coord_trans))+t))
    return new_coords

def create_Q_matrix(list_of_q:list[quat.quaternion]):
    """
    Creates the Q matrix that is needed for averageQuaternion.py.
    Q is a N x 4 matrix having the quaternions that shall be averaged as row vectors.

    :param list_of_q: Contains quaternion objects
    :return: Q matrix
    """
    q_arrays = []
    for q in list_of_q:
        q_arr = quat.as_float_array(q)
        q_arrays.append(q_arr)
    Q = np.matrix(q_arrays)
    return Q

def average_quats(Q):
    return avg(Q)
