import numpy as np
import pandas as pd
import argparse
from Bio.PDB import *
import math
def alphaAngle(poprReszta, reszta):
    v1 = poprReszta['O3\''].get_vector()
    v2 = reszta['P'].get_vector()
    v3 = reszta['O5\''].get_vector()
    v4 = reszta['C5\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def betaAngle(reszta):
    v1 = reszta['P'].get_vector()
    v2 = reszta['O5\''].get_vector()
    v3 = reszta['C5\''].get_vector()
    v4 = reszta['C4\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def gammaAngle(reszta):
    v1 = reszta['O5\''].get_vector()
    v2 = reszta['C5\''].get_vector()
    v3 = reszta['C4\''].get_vector()
    v4 = reszta['C3\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def deltaAngle(reszta):
    v1 = reszta['C5\''].get_vector()
    v2 = reszta['C4\''].get_vector()
    v3 = reszta['C3\''].get_vector()
    v4 = reszta['O3\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def epsilonAngle(reszta, nastReszta):
    v1 = reszta['C4\''].get_vector()
    v2 = reszta['C3\''].get_vector()
    v3 = reszta['O3\''].get_vector()
    v4 = nastReszta['P'].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def zetaAngle(reszta, nastReszta):
    v1 = reszta['C3\''].get_vector()
    v2 = reszta['O3\''].get_vector()
    v3 = nastReszta['P'].get_vector()
    v4 = nastReszta['O5\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def chiAngle(reszta):
    if(reszta.get_resname() in ['A', 'G']):
        v3 = reszta['N9'].get_vector()
        v4 = reszta['C4'].get_vector()
    else:
        v3 = reszta['N1'].get_vector()
        v4 = reszta['C2'].get_vector()
    v1 = reszta['O4\''].get_vector()
    v2 = reszta['C1\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def v0Angle(reszta):
    v1 = reszta['C4\''].get_vector()
    v2 = reszta['O4\''].get_vector()
    v3 = reszta['C1\''].get_vector()
    v4 = reszta['C2\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def v1Angle(reszta):
    v1 = reszta['O4\''].get_vector()
    v2 = reszta['C1\''].get_vector()
    v3 = reszta['C2\''].get_vector()
    v4 = reszta['C3\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def v2Angle(reszta):
    v1 = reszta['C1\''].get_vector()
    v2 = reszta['C2\''].get_vector()
    v3 = reszta['C3\''].get_vector()
    v4 = reszta['C4\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def v3Angle(reszta):
    v1 = reszta['C2\''].get_vector()
    v2 = reszta['C3\''].get_vector()
    v3 = reszta['C4\''].get_vector()
    v4 = reszta['O4\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def v4Angle(reszta):
    v1 = reszta['C3\''].get_vector()
    v2 = reszta['C4\''].get_vector()
    v3 = reszta['O4\''].get_vector()
    v4 = reszta['C1\''].get_vector()
    return math.degrees(calc_dihedral(v1, v2, v3, v4))

def getPDBFile(sciezka):
    czyPlikJestPoprawny = True
    try:
        parser = PDBParser()
        PDBFile = parser.get_structure("Struct", sciezka)
    except FileNotFoundError:
        czyPlikJestPoprawny = False  
        PDBFile = ""
    return czyPlikJestPoprawny, PDBFile

def getTorsionAngleMatrix(PDBFile):
    reszty= [reszta for reszta in PDBFile.get_residues()]
    resztyBezH = []
    for reszta in reszty:
        if(reszta.get_resname() in ['A', 'C', 'U', 'G']):
            resztyBezH.append(reszta)

    katy = np.array([])
    katyZReszty = []
    for id, reszta in enumerate(resztyBezH):
        if id != 0:
            katyZReszty.append(alphaAngle(resztyBezH[id-1],reszta))
        else:
            katyZReszty.append(None)
        katyZReszty.append(betaAngle(reszta))
        katyZReszty.append(gammaAngle(reszta))
        katyZReszty.append(deltaAngle(reszta))
        if id < (len(resztyBezH) - 1):
            katyZReszty.append(epsilonAngle(reszta, resztyBezH[id+1]))
            katyZReszty.append(zetaAngle(reszta, resztyBezH[id+1]))
        else:
            katyZReszty.append(None)
            katyZReszty.append(None)
        katyZReszty.append(chiAngle(reszta))
        if id == 0:
            katy = np.array(katyZReszty)
        else:
            katy = np.vstack([katy,katyZReszty])
        katyZReszty=[]
    return katy

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input PDB file with RNA structure.", required=True)
    parser.add_argument("-o", "--output", help="CSV output file.", required=True)
    args = parser.parse_args()
    czyPlikJestPoprawny, PDBFile = getPDBFile(args.input)
    if czyPlikJestPoprawny:
        katy = getTorsionAngleMatrix(PDBFile)
        DF = pd.DataFrame(katy)
        if not args.output.endswith(".csv"):
            args.output+=".csv"
        DF.to_csv(args.output, header=False, index=False)
