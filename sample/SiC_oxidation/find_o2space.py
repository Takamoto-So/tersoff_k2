import math
import random
import collections

AtomStr = collections.namedtuple("AtomStr", "x, vid, q, name")

def sphere_random():
    while True:
        X, Y = random.random() * 2 - 1, random.random() * 2 - 1
        S = X * X + Y * Y
        if S <= 1:
            break
    Z = 2 * S - 1
    S = math.sqrt((1 - Z * Z) / S)
    return (X * S, Y * S, Z)

def fillVoidArr(voidArr, atom, divide, blxyz, zMin, zRange):
    r = 2.0
    x1 = int((atom.x[0]-r)/blxyz[0]*divide)
    x2 = int((atom.x[0]+r)/blxyz[0]*divide)+1
    y1 = int((atom.x[1]-r)/blxyz[1]*divide)
    y2 = int((atom.x[1]+r)/blxyz[1]*divide)+1
    z1 = max(0, int((atom.x[2]-r-zMin)/zRange*divide))
    z2 = min(divide, int((atom.x[2]+r-zMin)/zRange*divide)+1)
    for i in range(x1, x2):
        i = i % divide
        for j in range(y1, y2):
            j = j % divide
            for k in range(z1, z2):
                #k = k % divide
                x = (1.0*i+0.5)*blxyz[0]/float(divide)
                y = (1.0*j+0.5)*blxyz[1]/float(divide)
                z = (1.0*k+0.5)*zRange/float(divide)+zMin
                dx = (x-atom.x[0]+0.5*blxyz[0])%blxyz[0]-0.5*blxyz[0]
                dy = (y-atom.x[1]+0.5*blxyz[1])%blxyz[1]-0.5*blxyz[1]
                dz = (z-atom.x[2]+0.5*blxyz[2])%blxyz[2]-0.5*blxyz[2]
                dr = math.sqrt(dx*dx+dy*dy+dz*dz)
                if dr < r:
                    voidArr[i][j][k] = 1

def find_o2_single(voidArr, divide, blxyz, zMin, zRange, tryMax):
    foundO2pos = False
    for O2tryCount in range(tryMax):
        while True:
            O1x, O1y, O1z = random.uniform(0, blxyz[0]), random.uniform(0, blxyz[1]), random.uniform(zMin, zMin+zRange)
            O1xi = int(O1x/blxyz[0]*divide)%divide
            O1yi = int(O1y/blxyz[1]*divide)%divide
            O1zi = max(min(divide-1, int((O1z-zMin)/zRange*divide)), 0)
            if voidArr[O1xi][O1yi][O1zi] == 0:
                break
        O2x, O2y, O2z = sphere_random()
        O2x *= 1.2
        O2y *= 1.2
        O2z *= 1.2
        O2x += O1x
        O2y += O1y
        O2z += O1z
        O2xi = int(O2x/blxyz[0]*divide)%divide
        O2yi = int(O2y/blxyz[1]*divide)%divide
        O2zi = max(min(divide-1, int((O2z-zMin)/zRange*divide)), 0)
        if voidArr[O2xi][O2yi][O2zi] == 0:
            foundO2pos = True
            break
    return(foundO2pos, [O1x, O1y, O1z], [O2x, O2y, O2z])


def find_o2space(atomArr, blxyz, zMin, zMax, nO2, tryMax):
    zRange = zMax-zMin

    divide = 100
    voidArr = [[[0 for i in range(divide)] for j in range(divide)] for k in range(divide)]
    for i, atom in enumerate(atomArr):
        fillVoidArr(voidArr, atom, divide, blxyz, zMin, zRange)

    ans = []
    for i in range(nO2):
        Oarr = []
        foundO2pos, O1, O2 = find_o2_single(voidArr, divide, blxyz, zMin, zRange, tryMax)
        if foundO2pos == False:
            break
        else:
            Oarr.append(AtomStr(O1, -1, 0.0, "O"))
            Oarr.append(AtomStr(O2, -1, 0.0, "O"))
            for atom in Oarr:
                fillVoidArr(voidArr, atom, divide, blxyz, zMin, zRange)
            ans.extend(Oarr)

    return ans
