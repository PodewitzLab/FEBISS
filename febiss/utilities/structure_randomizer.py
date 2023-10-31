from pymatgen.core import Molecule
from scipy.spatial.transform import Rotation as R

def structure_randomizer(num,ori_path,return_path):
    """
    :param num: number of randomized structures to be created
    :param ori_path:
    :param return_path:

    Takes an xyz-file from ori_path and writes a [num]ber of xyz-files with randomized structures into return_path.
    """

    euler = R.random(num).as_euler('zxy',degrees=True)
    r = [R.from_euler('zxz',euler[i]) for i in range(len(euler))]

    mol = Molecule.from_file(ori_path)
    labels = mol.labels
    cmol = mol.get_centered_molecule()

    for i in range(num):
        xyzlines = []
        text = open(ori_path, 'r').readlines()
        xyzlines.extend([text[0], text[1]])
        new_coords = r[i].apply(cmol.cart_coords)

        for k in range(len(new_coords)):
            line = "{0:<2} {1:>12}{2:>12}{3:>12}\n".format(labels[k], round(new_coords[k][0],5), round(new_coords[k][1],5), round(new_coords[k][2],5))
            xyzlines.append(line)
        #print(xyzlines,end='\n\n')
        with open(return_path.format(i),'w') as newfile:
            newfile.writelines(xyzlines)

structure_randomizer(50,'/home/lum/FEBISS/febiss/solvents/ACN.xyz','/home/lum/FEBISS/febiss/solvents/ACN/ACN.xyz')