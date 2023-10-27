import numpy as np
import quaternion as quat #installed on 06 Sept 2023 in febiss_pymatgen env via "conda install -c conda-forge quaternion" according to https://quaternion.readthedocs.io/en/latest/
import rmsd
from pymatgen.symmetry import analyzer as ana
from pymatgen.core import Molecule

class Util:
    @staticmethod
    def normalize(v1: np.ndarray):
        return v1/np.linalg.norm(v1)

class Representation:

    def __init__(self, mol: Molecule):
        self.mol = mol
        self.com = mol.center_of_mass
        self.red = self.choose_substructure(mol) #contains the indices of the chosen rigid atoms as well as the normalized relative coordinate vectors wrt the com
        self.rigid_atom_1_idx = self.red[0] #https://stackoverflow.com/questions/12646326/calling-a-class-function-inside-of-init
        self.rigid_atom_2_idx = self.red[1] #https://stackoverflow.com/questions/12646326/calling-a-class-function-inside-of-init
        self.quaternion = self.change_coord_sys()[1]


    """deprecated"""
    def choose_substructure(self, reference):
        for atom1 in range(len(self.mol.cart_coords)-1):
            for atom2 in range(atom1+1, len(self.mol.cart_coords)):
                at1_v = (self.mol.cart_coords[atom1] - self.com)/np.linalg.norm(self.mol.cart_coords[atom1] - self.com)
                at2_v = (self.mol.cart_coords[atom2] - self.com)/np.linalg.norm(self.mol.cart_coords[atom2] - self.com)
                if np.linalg.norm(np.cross(at1_v, at2_v)) < 0.1: #this corresponds to around 5.739° and 174.261°
                    if atom1 == len(self.mol.cart_coords)-1 and atom2 == len(self.mol.cart_coords):
                        print("Ohoh, the passed solvent is linear.")
                        raise NotImplementedError
                    continue
                else:
                    return atom1, atom2, at1_v, at2_v

    def choose_substructure_gist(self):
        pass

    def change_coord_sys(self):
        x_axis = self.red[2]
        y_axis = Util.normalize(np.cross(x_axis,self.red[3]))
        z_axis = Util.normalize(np.cross(x_axis,y_axis))
        transform_matrix = np.c_[x_axis,y_axis,z_axis] #https://stackoverflow.com/questions/27513246/combining-vectors-as-column-matrix-in-numpy
        return transform_matrix, quat.from_rotation_matrix(transform_matrix)


    def find_equivalent_rotations(self):
        pass


def main():
    name_list = ["ACN","ACT","BNZ","CCl4","CHX","CL3","CL4","DCM","DMF","DMS","ETL","HEX","MTL","NH3","OCT","PYR","THF","TOL"]
    for i in name_list:
        solvent = Representation(Molecule.from_file("/home/lum/Desktop/temp/{0}.xyz".format(i)))
        solvent.choose_substructure()
        com = solvent.com
        solv_pga = ana.PointGroupAnalyzer(solvent.mol)
        #print(solv_pga.get_pointgroup())
        #print(solv_pga.get_equivalent_atoms()["sym_ops"])
        #print(solv_pga.get_symmetry_operations())

        coords = {}
        coords_ori = {}
        for j in range(len(solvent.cart_coords)):
            coords[j] = tuple([solvent.labels[j],solvent.cart_coords[j]])
            coords_ori[j] = tuple([solvent.labels[j], solvent.cart_coords[j]-com])
            #print("{1}{0} = {2}".format(i,coords[i][0],coords[i][1]))

        #print("\nCOM was set to (0,0,0)! The new coordinates are:")
        #for i in range(len(solvent.cart_coords)):
            #print("{1}{0} = {2}".format(i, coords_ori[i][0], coords_ori[i][1]))
        #print()

        #print(solv_pga.get_symmetry_operations())
        symmops_list = solv_pga.get_symmetry_operations()
        num = 1
        coord_dict = {}
        quat_dict = {}
        with open("/home/lum/Desktop/temp/{0}_rot.txt".format(i),'w') as rotfile:
            rotfile.write("ROTMATRICES START\n")
            for symm in symmops_list: #filtering adopted from analyzer.get_rotational_symmetry_number() (line 1284)
                coord_list = []
                rot = symm.rotation_matrix
                if np.abs(np.linalg.det(rot) - 1) < 1e-4:
                    rotfile.write("rotation {0}:\n".format(num))
                    rotfile.write("{0}\n".format(rot))
                    quat_dict[num] = quat.from_rotation_matrix(rot)
                    for j in range(len(solvent.cart_coords)):
                        coord_list.append("{0}{1}_{3} = {2}\n".format(coords_ori[j][0], j, symm,symm.apply_rotation_only(coords_ori[j][1]),num))
                    coord_dict[num] = coord_list
                    num += 1
            rotfile.write("ROTMATRICES STOP. TOTAL NUMBERS: {0}\n\n".format(num-1))
            rotfile.write("COORDINATES START\n")
            keylist = coord_dict.keys()
            for k in keylist:
                rotfile.write("rotation {0}:\n".format(k))
                for l in coord_dict[k]:
                    rotfile.write("{0}".format(l))
            rotfile.write("COORDINATES STOP\n\n")
            rotfile.write("QUATERNIONS START\n")
            keylist = quat_dict.keys()
            for m in keylist:
                rotfile.write("rotation {0}: ".format(m))
                rotfile.write("{0}\n".format((quat_dict[m])))
            rotfile.write("QUATERNIONS STOP\n\n")
            rotfile.write("")
        print("{0} end".format(i))