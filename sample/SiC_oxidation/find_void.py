import os
import glob
import re
import collections
import read_cfg
import find_o2space
import atom_network

AtomStr = collections.namedtuple("AtomStr", "x, vid, q, name")

dumpfiles = glob.glob(os.getcwd()+"/cfg/run.*.cfg")
dumpfiles_num = [(re.search("[0-9]+", re.search("\.[0-9]+\.cfg", x).group()).group(), x) for x in dumpfiles]
dumpfiles_num.sort(key = lambda x:int(x[0]))

timestep = dumpfiles_num[-1][0]
dumpfilename = dumpfiles_num[-1][1]

atomArr, blxyz = read_cfg.readCFG(dumpfilename)
atomIdDict = {x.vid:x for x in atomArr}

atom_net = atom_network.findAtomNetwork(atomArr, blxyz, False)

O2_dict = {"O":2}
O_dict = {"O":1}

co_co2_count = 0
delete_atoms_CO = []
delete_atoms_vac = []

atomMainSet = set()

atom_unsearched = set([x.vid for x in atomArr])
while len(atom_unsearched) > 0:
    atom_found = [atom_unsearched.pop()]
    atom_connected = set(atom_found)
    while len(atom_found) > 0:
        last = atom_found.pop()
        for atom2 in atom_net.get(last, []):
            if not(atom2[0] in atom_connected):
                atom_connected.add(atom2[0])
                atom_unsearched.remove(atom2[0])
                atom_found.append(atom2[0])
    atomNameDict = {}
    xm = [0.0]*3
    for atom in atom_connected:
        atomNameDict[atomIdDict[atom].name] = atomNameDict.get(atomIdDict[atom].name, 0)+1
        xm[0] += atomIdDict[atom].x[0]/len(atom_connected)
        xm[1] += atomIdDict[atom].x[1]/len(atom_connected)
        xm[2] += atomIdDict[atom].x[2]/len(atom_connected)

    if len(atom_connected) > 20:
        atomMainSet = atomMainSet.union(set(atom_connected))
    if len(atom_connected) < 7:
        if not "Si" in atomNameDict and "C" in atomNameDict and "O" in atomNameDict:
            co_co2_count += 1
            delete_atoms_CO.append((atom_connected, atomNameDict))
        elif atomNameDict == O2_dict or atomNameDict == O_dict:
            delete_atoms_vac.append((atom_connected, atomNameDict))

print("CO/CO2: ", co_co2_count)
delete_atoms = []
nCtotal = 0
for mol, nameDict in delete_atoms_CO:
    delete_atoms.extend(list(mol))
    nCtotal += nameDict["C"]
    nC = str(nameDict["C"])
    nO = str(nameDict["O"])
    if nC == "1":
        nC = ""
    if nO == "1":
        nO = ""
    print("Delete C-O mol: C"+nC+"O"+nO)
print("Delete C total: "+str(nCtotal))
delete_vac_count = 0
for mol, nameDict in delete_atoms_vac:
    delete_atoms.extend(list(mol))
    delete_vac_count += len(mol)
print("Delete O/O2 (vac): ", delete_vac_count)

atomArr_deleted = [x for x in atomArr if (x.vid in atomMainSet)]
atomArr_deleted_O2 = sorted([x for x in atomArr_deleted if x.name == "O"], key=lambda atom:atom.x[2])
zMin = atomArr_deleted_O2[5].x[2]    + 1.0
zMax = atomArr_deleted_O2[-10].x[2]- 6.0
zMergin = 1.0

inner_nrho_single = 1.0/blxyz[0]/blxyz[1]/(zMax-zMin+2.0*zMergin)
inner_atoms = int(0.001/inner_nrho_single)
print("Insert O2 range: ["+str(zMin)+", "+str(zMax)+"]")
print("Insert O atoms: "+str(inner_atoms*2))

Oatoms = find_o2space.find_o2space(atomArr, blxyz, zMin, zMax, inner_atoms, 10000)

if len(Oatoms) > 0:
    for i in range(len(Oatoms)//2):
        print("O["+str(i)+"]: ", Oatoms[2*i].x, Oatoms[2*i+1].x)
else:
    print("No O2 insertion.")

fwts = open("cfg/in.o2add."+str(timestep), "w")
if len(delete_atoms) > 0:
    fwts.write("group gdelete subtract all all\n")
    for atom in delete_atoms:
        fwts.write("group gdelete id "+str(atom)+"\n")
    fwts.write("delete_atoms group gdelete\n")
    fwts.write("group gdelete delete\n")
for atom in Oatoms:
    fwts.write("create_atoms 3 single "+str(atom.x[0])+" "+str(atom.x[1])+" "+str(atom.x[2])+" remap yes units box\n")
else:
    fwts.write("")
fwts.close()

