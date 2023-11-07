start = "@<TRIPOS>ATOM"
stop = "@<TRIPOS>BOND"
name_list = ["ACN","ACT","BNZ","CCl4","CHX","CL3","CL4","DCM","DMF","DMS","ETL","HEX","MTL","NH3","OCT","PYR","THF","TOL"]
path = input("Full path of folder which contains solvent mol2 files: ")
path_mol2 = path+"/{0}.mol2"
path_xyz = path+"/{0}.xyz"
for i in name_list:
    with open(path_mol2.format(i),'r') as mol2file:
        lines = mol2file.readlines()
        start_idx = 0
        stop_idx = 0
        for j in range(len(lines)):
            if lines[j].strip()==start:
                start_idx = j
            elif lines[j].strip()==stop:
                stop_idx = j
                break

    xyz_lines = []
    xyz_lines.append("{0}\n".format(lines[stop_idx-1].split()[0]))
    xyz_lines.append("{0}\n".format(lines[1].strip()))
    for k in range(start_idx+1,stop_idx):
        parts = lines[k].split()
        parts[1] = parts[1].strip("0123456789")
        line = "{0:<2}{1:>12}{2:>12}{3:>12}\n".format(parts[1],parts[2],parts[3],parts[4])
        xyz_lines.append(line)

    with open(path_xyz.format(i),'w') as xyzfile:
        xyzfile.writelines(xyz_lines)

    print("{0} end".format(i))
