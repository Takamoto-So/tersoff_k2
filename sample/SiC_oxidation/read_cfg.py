import functools
import collections

AtomStr = collections.namedtuple("AtomStr", "x, vid, q, name")

def readCFG(filename):

    cfg_words = [".NO_VELOCITY.", ".DUMMY."]

    fr = open(filename, "r")
    
    isBoundData = False
    isAtomData = False
    cfgData = {}
    atomArr = []
    colCount = 0
    cfgSetting = True
    blxyz = []
    atomName = ""
    for line in fr.readlines():
        if functools.reduce(lambda x,y: x or y, map(line.startswith, cfg_words)):
            continue
        if cfgSetting:
            eqind = line.find("=")
            if eqind > -1:
                cfgData[line[:eqind-1]] = line[eqind+1:-1]
                continue
            else:
                cfgSetting = False
                blxyz = [float(x.split()[0]) for x in [cfgData["H0(1,1)"], cfgData["H0(2,2)"], cfgData["H0(3,3)"], cfgData["H0(2,1)"], cfgData["H0(3,1)"], cfgData["H0(3,2)"]]]
    
        words = line.split()
        if colCount == 0:
            None
        elif colCount == 1:
            atomName = words[0]
        else:
            atomArr.append((int(words[3]), float(words[0]), float(words[1]), float(words[2]), float(words[6]), atomName))
        colCount = (colCount+1)%3
    atomArr.sort(key=lambda x:x[0])

    posToV = lambda x:AtomStr((x[1]*blxyz[0]+x[2]*blxyz[3]+x[3]*blxyz[4], x[2]*blxyz[1]+x[3]*blxyz[5], x[3]*blxyz[2]), x[0], x[4], x[5])
    atomArr_v = [posToV(x) for x in atomArr]
    return atomArr_v, blxyz

