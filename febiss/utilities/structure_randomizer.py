from pymatgen.core import Molecule
from scipy.spatial.transform import Rotation as R

def structure_randomizer(num,ori_path,return_path):
    """
    :param num: number of randomized structures to be created
    :param ori_path: path of the original xyz file
    :param return_path: path of the files to be created. accepts one placeholder {0} for an index that goes from 0 to num

    Takes an xyz-file from ori_path and writes a [num]ber of xyz-files with randomized structures into return_path.
    """

    r = R.random(num)

    mol = Molecule.from_file(ori_path)
    labels = mol.labels
    cmol = mol.get_centered_molecule()

    for i in range(num):
        new_coords = r[i].apply(cmol.cart_coords)
        write_xyz(ori_path,return_path,new_coords,i,labels)


def write_xyz(ori_path,return_path,new_coords, idx=0, labels = None):
    """
    :param ori_path: ori_path: path of the original xyz file that acts as a template
    :param return_path: path of the files to be created. accepts one placeholder {0} for an index that goes from 0 to num
    :param new_coords: newly created coordinates after rotation, translation, ...
    :param idx: Optional. Used for numbering the created xyz-files
    :param labels: Optional. If labels of atoms are already known they can be passed as list. Otherwise this function
           will determine them from the xyz file using the pymatgen.core package "Molecule".
    :return path of written file


    Takes an xyz-file (ori_path) as template to create an xyz-file (return_path) with the new_coords. Return_path
    accepts a placeholder {0} for numbering with idx (default: 0).
    """

    if labels is None:
        from pymatgen.core import Molecule
        mol = Molecule.from_file(ori_path)
        labels = mol.labels

    xyzlines = []
    text = open(ori_path, 'r').readlines()
    xyzlines.extend([text[0], text[1]])
    for k in range(len(new_coords)):
        line = "{0:<2} {1:>12}{2:>12}{3:>12}\n".format(labels[k], round(new_coords[k][0], 5),
                                                       round(new_coords[k][1], 5), round(new_coords[k][2], 5))
        xyzlines.append(line)
    with open(return_path.format(idx),'w') as newfile:
        newfile.writelines(xyzlines)

    return return_path.format(idx)