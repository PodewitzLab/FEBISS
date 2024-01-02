import os

import quaternion

from utilities.structure_randomizer import write_xyz as write
from pymatgen.symmetry import analyzer as ana
from pymatgen.core import Molecule
import numpy as np
import quaternion as quat
import quat_calc
from utilities.averageQuaternions import averageQuaternions as avg
from typing import Union


def setup(ori_path, return_path, write_out=False, load=False, i1 = 3, i2 = 4) -> (dict[quat.quaternion], Molecule):
    """
    Takes a molecular xyz-file (ori_path) centers it (puts COM at origin) and analyzes the symmetry elements using the
    pymatgen.symmetry package. Only the rotation matrices are then converted into quaternions which are stored in
    quat_dict. quat_dict is a dictionary with keys from 1 to N(rotations) since index 0 should be left for the original
    file.

    :param ori_path: Path to reference xyz-file.
    :param return_path: Used only if write_out == True. Path of created file
    :param write_out: Default: False. If true, equivalent structures will be written out to return_path.
    :param load: Default: False. If true, the written out molecules get read in again as separate Molecule objects.
    :param i1: Index of first rigid atom.
    :param i2: Index of second rigid atom.
    :return: quat_dict, cref, eq_struc_dict, char_quat_dict
    """

    ref = Molecule.from_file(ori_path)
    cref = ref.get_centered_molecule()
    cref_pga = ana.PointGroupAnalyzer(cref)
    symmops_list = cref_pga.get_symmetry_operations()

    num = 1
    quat_dict = {}
    eq_struc_dict = {}
    char_quat_dict = {}
    for symm in symmops_list:  # filtering adopted from analyzer.get_rotational_symmetry_number() (line 1284)
        coord_list = []
        rot = symm.rotation_matrix
        if np.abs(np.linalg.det(rot) - 1) < 1e-4:
            quat_dict[num] = quat.from_rotation_matrix(rot)
            for j in range(len(ref.cart_coords)):
                coord_list.append(symm.apply_rotation_only(ref.cart_coords[j]))
            if write_out == True:
                written_path = write(ori_path, return_path, coord_list, idx=num)

                if load:
                    eq_struc_dict[num] = Molecule.from_file(written_path)
                    char_quat_dict[num] = create_gist_quat(eq_struc_dict[num],i1,i2)

            num += 1


    return quat_dict, cref, eq_struc_dict, char_quat_dict

def inv(q:quaternion.quaternion):
    """
    Returns q^(-1)

    :param q: quaternion object
    :return: inverse of q
    """
    return quaternion.quaternion.inverse(q)

def norm(q:quaternion.quaternion):
    """
    Returns ||q||

    :param q: quaternion object
    :return: norm of q
    """
    return quaternion.quaternion.norm(q)

def unity(q:quaternion.quaternion):
    """
    Sets ||q|| -> 1

    :param q: quaternion object
    :return: normalized q
    """
    return q/norm(q)

def create_gist_quat(mol: Molecule, i1 : int, i2 : int):
    """
    Calculates characteristic quaternion from positions of center of mass and of atoms i1 and i2 the way the quaternion
    is calculated in cpptraj GIST.

    :param mol: Molecule object
    :param i1: index number of atom 1
    :param i2: index number of atom 2
    :return: quaternion object
    """
    at1 = (mol.cart_coords[i1] - mol.center_of_mass) / np.linalg.norm(mol.cart_coords[i1] - mol.center_of_mass)  # H1
    at2 = (mol.cart_coords[i2] - mol.center_of_mass) / np.linalg.norm(mol.cart_coords[i2] - mol.center_of_mass)  # H2
    q = quat.from_float_array(np.array(quat_calc.gist_quat(at1, at2)))
    return q


def new_coord_gen(mol:Molecule, q:quaternion.quaternion):
    """
    Calculates the new coordinates after rotation of mol by q using the formula:
    [0,x'] = q * [0,x] * q^(-1)

    :param mol: Molecule object
    :param q: quaternion object
    :return: new coords as list of arrays
    """
    q = unity(q)
    new_coords = []
    for m in range(len(mol.cart_coords)):
        q_coord = unity(quat_calc.vec_to_quat(mol.cart_coords[m]))
        q_coord_trans = unity(q * q_coord * inv(q))
        new_coords.append(quat_calc.vec_part(q_coord_trans))
    return new_coords

def distance(q1:quaternion.quaternion, q2:quaternion.quaternion):
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

    q1_array = quaternion.as_float_array(q1)
    q2_array = quaternion.as_float_array(q2)
    scalar_product = np.abs(q1_array[0]*q2_array[0]+
                            q1_array[1]*q2_array[1]+
                            q1_array[2]*q2_array[2]+
                            q1_array[3]*q2_array[3])
    theta = 2*np.arccos(scalar_product)
    return theta

""" average_q is currently not in use. 
an implementation based on the power method to find the principal eigenvector will follow"""
# def average_q(list_of_q: list[quaternion.quaternion]):
#     """
#     This power method implementation to find the principal eigenvector of the matrix M (see:
#     Markley, F. L.; Cheng, Y.; Crassidis, J. L.; Oshman, Y. Averaging Quaternions.
#     Journal of Guidance, Control, and Dynamics 2007, 30 (4), 1193–1197. https://doi.org/10.2514/1.28949.)
#     was adopted from https://pythonnumericalmethods.berkeley.edu/notebooks/chapter15.02-The-Power-Method.html (accessed
#     06 Nov 2023).
#
#     :param list_of_q: List of quaternions to be averaged
#     :return: Average quaternion
#     """
#     import numpy as np
#
#     def normalize(x):
#         fac = abs(x).max()
#         x_n = x / x.max()
#         return fac, x_n
#
#     x = np.array([1,0,0,0])
#
#     #https://stackoverflow.com/questions/17428621/python-differentiating-between-row-and-column-vectors (6 Nov 2023),
#     #https://stackoverflow.com/questions/34179413/numpy-column-and-row-vectors (6 Nov 2023),
#     M = np.zeros((4,4)) #https://datagy.io/numpy-zeros/ (6 Nov 2023)
#     for q in list_of_q:
#         q = np.reshape(quat.as_float_array(q),(4,1))
#         M = M + np.dot(q,np.transpose(q))
#
#     for i in range(50):
#         x = np.dot(M, x)
#         lambda_1, x = normalize(x)
#
#
#     print('Eigenvalue:', lambda_1)
#     print('Eigenvector:', x)
#     eigenValues, eigenVectors = np.linalg.eig(M)
#     eigenVectors = eigenVectors[:,eigenValues.argsort()[::-1]]
#     print(np.real(eigenVectors[:,0]))

def create_Q_matrix(list_of_q:list[quaternion.quaternion]):
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


"""
--------------------------------------------------------------------------------------------
------------------------------------TESTING------------------------------------------------
--------------------------------------------------------------------------------------------
"""

"""
Just applying the quaternions which transform the equivalent structures into the reference on the arbitrary
structure obviously does not lead to equivalent structures.
"""
def test_1(ori_path, test_path, return_path):
    from scipy.spatial.transform import Rotation as R

    stp = setup(ori_path,return_path)
    test = Molecule.from_file(test_path)
    for i in range(len(stp[0].keys())):
        q = quat.as_float_array(stp[0][i+1])
        write(ori_path,
              return_path,
              R.from_quat(q).apply(test.cart_coords),
              i+1,
              test.labels)

"""
Therefore, we first have to find the quaternion that transforms the characteristic quaternion from the arbitrarily
oriented molecule into the characteristic quaternion of the reference.
q2 = q1*qt
q1_inv * q2 = qt
"""

def test_2(ori_path,return_path):
    stp = setup(ori_path,return_path)
    q1 = stp[0][1]
    q2 = stp[0][2]
    q1_inv = quaternion.quaternion.inverse(stp[0][1])
    qt = q1_inv * q2
    q2_test = q1 * qt
    print("q2: {0}".format(q2))
    print("q2_test: {0}".format(q2_test))


def test_3(ori_path,test_path,i1=3,i2=4):
    ori = Molecule.from_file(ori_path)
    test = Molecule.from_file(test_path)
    q_ori = create_gist_quat(ori,i1,i2)
    q_test = create_gist_quat(test,i1,i2)
    stp = setup(ori_path,test_path)
    print(stp[0][1])
    print(stp[0][2])

    q_ori_inv = quat.quaternion.inverse(q_ori)
    qt = q_ori_inv * q_test
    q2_test = q_ori*qt
    print("q_test: {0}".format(q_test))
    print("q2_test: {0}".format(q2_test))

"""
How to apply the rotation quaternions to get equivalent structures. 
Apparently, running scipy.spatial.transform.Rotation.from_quat(qt).apply(coords) does not calculate R(x) = q*[0,x]*q' 
(see e.g. Karney, C. Journal of Molecular Graphics and Modelling 25 (2007) 595–604 (p.596)). Maybe this also explains
why calling R.from_quat(quaternion(1,0,0,0)).as_matrix gives [[1 0 0], [0 -1 0], [0 0 -1]] and not the identity matrix.
Keeping the quaternions at unit length is mandatory.
"""
def test_4(ori_path,return_path,i1=3,i2=4,write_out=False):
    from scipy.spatial.transform import Rotation as R
    stp = setup(ori_path,return_path)
    quat_dict = stp[0]
    num = 1

    #this does not work
    for q in list(quat_dict.values()):
        new_coords = R.from_quat(quaternion.as_float_array(q)).apply(stp[1].cart_coords)
        if write_out:
            write(ori_path,return_path+'_wrong.xyz',new_coords,num,stp[1].labels)
        print(num)
        print(new_coords,end='\n\n')
        num += 1

    num = 1
    #this works
    for q in list(quat_dict.values()):
        new_coords = new_coord_gen(stp[1],q)
        if write_out:
            write(ori_path,return_path+'_true.xyz',new_coords,num,stp[1].labels)
        print(num)
        print(np.matrix(new_coords),end='\n\n')
        num += 1

"""
Here it is shown how to transform the original structure into an equivalent structure by using the characteristic
quaternion.
"""

def relation_between_quats_of_eq_structs(ori_path,test_path,return_path,write_out = False):
    ori = Molecule.from_file(ori_path) #the original structure
    test = Molecule.from_file(test_path) #equivalent structure
    q_ori = create_gist_quat(ori,3,4)
    q_test = create_gist_quat(test,3,4)

    qt = q_test * inv(q_ori)
    new_coords = new_coord_gen(ori,qt)

    if write_out:
        write(ori_path,return_path,new_coords,0,ori.labels)

    else:
        print(new_coords)

"""
Here it is finally shown how to get the equivalent quaternion structure for an arbitrary orientation
"""

def test_with_randoms(ori_path,test_path,return_path,i1=3,i2=4,write_out=False):
    stp = setup(ori_path,return_path)
    ori = Molecule.from_file(ori_path)  # the original structure
    test = Molecule.from_file(test_path)  # random structure
    ref_quats = stp[0]  # holds quaternions that transform ori into equivalent structures

    for rot in range(len(ref_quats)):

        q_ori = create_gist_quat(ori,i1,i2)
        q_test = create_gist_quat(test,i1,i2)

        qt = q_test * inv(q_ori) * ref_quats[rot+1] #+1 as always because ref_quats is a dict with keys 1 - [number of rots]

        new_coords = new_coord_gen(ori, qt)
        if write_out:
            write(ori_path, return_path, new_coords, rot, ori.labels)
            print('rot: {0}'.format(rot))
            print("equivalent quaternion by quaternion multiplication:\n {0}".format(qt * q_ori))
            print("equivalent quaternion by calculating the characteristic quaternion from the created file:\n {0}".
            format(create_gist_quat(Molecule.from_file(return_path.format(rot)),i1,i2)), end="\n\n")

        else:
            print('rot: {0}'.format(rot))
            print("equivalent quaternion:\n {0}".format(qt * q_ori))
            print("new_coords: {0}".format(np.matrix(new_coords)), end="\n\n")


def distance_test(ori_path,ref_path_list: Union[list[str],str],i1 = 3, i2 = 4):

    ori = Molecule.from_file(ori_path)
    q_ori = create_gist_quat(ori,i1,i2)
    print("q_ori: {0}".format(q_ori))

    if type(ref_path_list) == str:
        ref_path_list = [ref_path_list]

    for i in ref_path_list:
        ref = Molecule.from_file(i)
        q_ref = create_gist_quat(ref, i1, i2)
        print("q_ref: {0}".format(q_ref))
        print("theta: {0}".format(distance(q_ori, q_ref)))


"""
This is a test of the actual clean up done in FEBISS. A random test orientation of a solvent gets compared to the 
characteristic quaternions of the equivalent structures. This is done by measuring the quaternion distance to each
equivalent structure. Two cases can then apply:
1) The test structure is closest to the original structure -> nothing happens
2) The test structure is closest to an equivalent structure -> the test structure gets rotated by the same rotation that
is needed to convert the equivalent structure to the original structure, i.e. the inverse of the rotation quaternion
needed to rotate the original structure to the equivalent structure.
"""

def clean_up_test(ori_path,test_path,refdir_path, i1 = 3, i2 = 4, abb = 'ACN'):
    ori = Molecule.from_file(ori_path)
    ori_char_q = create_gist_quat(ori,i1,i2)
    test = Molecule.from_file(test_path)
    test_char_q = create_gist_quat(test,i1,i2)
    path_template = refdir_path + "/{0}".format(abb) + "_{0}.xyz"  # last placeholder is for enumeration of files used in write()
    stp = setup(ori_path,path_template,write_out=True,load=True, i1=i1, i2=i2)
    key_list = list(stp[3].keys())
    distance_list = []

    for j in range(len(key_list)):
        distance_list.append(distance(test_char_q, stp[3][j+1]))
    print("rotation quaternions: {0}".format(stp[0]))
    print("distance_list before: {0}".format(distance_list))

    min_idx = distance_list.index(min(distance_list))
    print("min_idx = {0}".format(min_idx))

    print("test_char_q before: {0}".format(test_char_q))
    test_char_q = test_char_q * inv(ori_char_q) * inv(stp[0][min_idx+1]) * ori_char_q
    print("test_char_q after: {0}".format(test_char_q))

    distance_list = []
    for j in range(len(key_list)):
        distance_list.append(distance(test_char_q, stp[3][j + 1]))
    print("distance_list after: {0}".format(distance_list))


def average_test(ori_path,theta:float,i1=3,i2=4, return_path = None):
    """
    Rotates a given molecule at ori_path for theta degrees around the z-axis. If a return path is given it also creates
    xyz-files of the rotated molecules.

    :param ori_path: Path of the original molecule
    :param theta: rotation angle
    :param i1: rigid atom 1 for characteristic quaternion calculation
    :param i2: rigid atom 2 for characteristic quaternion calculation
    :param return_path: Is None by default. If given, prints out xyz-files of rotated molecules.
    :return: Prints out and xyz-files to return path if return path is given.
    """
    ori = Molecule.from_file(ori_path)
    q_ori = create_gist_quat(ori, i1, i2)

    ori.rotate_sites(theta=theta)
    q_1 = create_gist_quat(ori, i1, i2)
    if return_path is not None:
        write(ori_path,return_path,ori.cart_coords,1)

    ori.rotate_sites(theta=theta)
    q_2 = create_gist_quat(ori, i1, i2)
    if return_path is not None:
        write(ori_path,return_path,ori.cart_coords,2)

    print("q_ori: {0}".format(q_ori))
    print("q_1: {0}".format(q_1))
    print("q_2: {0}".format(q_2))

    Q = create_Q_matrix([q_ori,q_1,q_2])
    print(quat.from_float_array(avg(Q)))