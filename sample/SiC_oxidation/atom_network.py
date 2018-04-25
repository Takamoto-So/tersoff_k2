import bisect
import math
import itertools
import collections

AtomStr = collections.namedtuple("AtomStr", "x, vid, q, name")

Atom_Radius_Raw = {"C": 0.72, "Si": 1.07, "H":0.35, "O":0.6}
Atom_Radius = {}
Atom_Radius_Max = 0.0
for key, value in Atom_Radius_Raw.items():
    radius = value*(0.5+0.5**0.5)
    Atom_Radius[key] = radius
    Atom_Radius_Max = max(Atom_Radius_Max, radius)

KdNode = collections.namedtuple("KdNode", "divide_value, divide_key, left, right, data")

def makeKdTree(vertices, region, indexes=None):
    if indexes == None:
        indexes = list(range(len(vertices)))
    if region == None:
        region = [(float("-inf"), float("inf"))]*len(vertices[0].x)

    if len(indexes) == 0:
        return None

    divide_key = 0
    region_max = 0
    for i in range(len(region)):
        if region_max < region[i][1]-region[i][0]:
            region_max = region[i][1]-region[i][0]
            divide_key = i
    if len(indexes) < 32:#region[divide_key][1]-region[divide_key][0] < 5.0:
        return KdNode(-1, -1, None, None, indexes)
    indexes.sort(key=lambda i:vertices[i].x[divide_key])
    index_mid = len(indexes)//2
    divide_value = vertices[indexes[index_mid]].x[divide_key]

    region_left = region[:]
    region_right = region[:]
    region_left[divide_key] = (region[divide_key][0], divide_value)
    region_right[divide_key] = (divide_value, region[divide_key][1])
    kd_left = makeKdTree(vertices, region_left, indexes[:index_mid])
    kd_right = makeKdTree(vertices, region_right, indexes[index_mid+1:])

    ans = KdNode(divide_value, divide_key, kd_left, kd_right, [indexes[index_mid]])
    return ans

def printKdTree(kdTree, depth=0):
    print(" "*depth, kdTree.vertex.x, kdTree.divide_key)
    if kdTree.left != None:
        printKdTree(kdTree.left, depth+1)
    if kdTree.right != None:
        printKdTree(kdTree.right, depth+1)

def findVertices(kdTree, vertices, region):
    ans = []
    for kd_datum in kdTree.data:
        for v_x, r_x in zip(vertices[kd_datum].x, region):
            if not(r_x[0] <= v_x and v_x <= r_x[1]):
                break
        else:
            ans.append(kd_datum)
    if kdTree.left != None and region[kdTree.divide_key][0] < kdTree.divide_value:
        ans.extend(findVertices(kdTree.left, vertices, region))
    if kdTree.right != None and region[kdTree.divide_key][1] > kdTree.divide_value:
        ans.extend(findVertices(kdTree.right, vertices, region))
    return ans

def findAtomNetwork(atomArr, blxyz, zflag):
    atomArr_9x = []
    for atom in atomArr:
        if zflag:
            for xlat, ylat, zlat in itertools.product(range(-1,2), repeat=3):
                atomArr_9x.append(AtomStr((atom.x[0]+xlat*blxyz[0]+ylat*blxyz[3]+zlat*blxyz[4], atom.x[1]+ylat*blxyz[1]+zlat*blxyz[5], atom.x[2]+zlat*blxyz[2]), atom.vid, atom.q, atom.name))
        else:
            zlat = 0
            for xlat, ylat in itertools.product(range(-1,2), repeat=2):
                atomArr_9x.append(AtomStr((atom.x[0]+xlat*blxyz[0]+ylat*blxyz[3]+zlat*blxyz[4], atom.x[1]+ylat*blxyz[1]+zlat*blxyz[5], atom.x[2]+zlat*blxyz[2]), atom.vid, atom.q, atom.name))
    region_global = []
    for i in range(3):
        xmin = float("-inf")
        xmax = float("inf")
        #for atom in atomArr_9x:
            #xmin = min(xmin, atom.x[i])
            #xmax = max(xmax, atom.x[i])
        region_global.append((xmin, xmax))
    kdTree = makeKdTree(atomArr_9x, region_global)
    atom_network = {}
    for atom in atomArr:
        len_max = 2.0*Atom_Radius_Max
        region = []
        for i in range(len(atom.x)):
            if i == 0:
                region.append((atom.x[i], atom.x[i]+len_max))
            else:
                region.append((atom.x[i]-len_max, atom.x[i]+len_max))
        neigh_rect = findVertices(kdTree, atomArr_9x, region)
        for atom2_index in neigh_rect:
            atom2 = atomArr_9x[atom2_index]
            if atom2.vid == atom.vid:
                continue
            if atom2.x[0] == atom.x[0]:
                if atom2.x[1] > atom.x[1]:
                    continue
                if atom2.x[1] == atom.x[1]:
                    if atom2.x[2] > atom.x[2]:
                        continue
            con_max = 0.0
            r1 = Atom_Radius.get(atom.name, 0.0)
            r2 = Atom_Radius.get(atom2.name, 0.0)
            con_max = r1+r2
            dz = atom2.x[2]-atom.x[2]
            if dz<-con_max or dz>con_max:
                continue
            dy = atom2.x[1]-atom.x[1]
            if dy<-con_max or dy>con_max:
                continue
            dx = atom2.x[0]-atom.x[0]
            if dx<-con_max or dx>con_max:
                continue
            iaLength = math.sqrt(dx*dx+dy*dy+dz*dz)

            if iaLength < con_max:
                atom_network.setdefault(atom.vid, []).append((atom2.vid, dx, dy, dz))
                atom_network.setdefault(atom2.vid, []).append((atom.vid, -dx, -dy, -dz))
    return atom_network

