import glob
import os
from rdkit import Chem
from rdkit.Chem import MolToMolFile, MolFromMolFile, AllChem, MolToMolBlock
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMolecule
from rdkit.Chem.rdMolTransforms import SetDihedralDeg, GetDihedralDeg
class MolToGjfFile:
    def __init__(self, mol, fileName, directory):
        self.mol = mol

        self.cutBlock()
        self.gjfFromBlock(fileName = fileName)
        self.writeGjf(directory, fileName)

    def cutBlock(self):
        self.shortBlock = []

        molBlock = MolToMolBlock(self.mol).splitlines()[4:]

        for i in range(len(molBlock) - 1, 0, -1):
            if len(molBlock[i].split()) < 5:
                self.shortBlock = molBlock[:i]

    def gjfFromBlock(self, fileName):
        self.gjfArray = []

        self.gjfArray.append(f"%nprocshare = 8\n"  # 12 for G7, 8 for G6
                             f"%mem = 16\n"  # 24 for G7, 16 for G6
                             f"%chk={fileName}.chk\n"
                             f"#p opt=vtight b3lyp/6-311++g(d,p) nosymm int=superfinegrid scf=verytight\n\n"
                             f"[No Title]\n\n"
                             f"0 1\n")

        for string in self.shortBlock:
            string = string.split()
            self.gjfArray.append(f"{string[3]} {string[0]} {string[1]} {string[2]}\n")

        self.gjfArray.append(f"\n\n\n")

    def writeGjf(self, directory, fileName):
        with open(rf"{directory}\{fileName}.gjf", "a") as newFile:
            for string in self.gjfArray:
                newFile.write(string)
            pass
def readSdf(directory):

    sdfFilePath = glob.glob(fr"{directory}\*.sdf")[0]
    mols = Chem.SDMolSupplier(sdfFilePath)

    return mols

def findAngles(mol):
    rotatable = []
    dihAnglesList = []

    def fingRotatebleBond():
        rotatableRingBonds = Chem.MolFromSmarts("[!$(*#*)&X3&R]-!@[!$(*#*)&X3&R]")
        rotatableMatches = mol.GetSubstructMatches(rotatableRingBonds)

        for bond in rotatableMatches:
            if bond[0] < bond[1]:
                rotatable.append((bond[0], bond[1]))
            else:
                rotatable.append((bond[1], bond[0]))

    def getTorsion():
        for bond in rotatable:
            fourAtomBond = []

            for atom in mol.GetAtomWithIdx(bond[0]).GetNeighbors():
                idx = atom.GetIdx()

                if idx != bond[1]:
                    fourAtomBond.append(idx)
                    fourAtomBond.append(bond[0])
                    break
            for atom in mol.GetAtomWithIdx(bond[1]).GetNeighbors():
                idx = atom.GetIdx()

                if idx != bond[0]:
                    fourAtomBond.append(bond[1])
                    fourAtomBond.append(idx)
                    break

            dihAnglesList.append(fourAtomBond)

    fingRotatebleBond()
    getTorsion()
    return dihAnglesList

def processMol(mol, dihAnglesList):
    def findAnglesDeg(conf, dihAnglesList):
        originalAngles = []

        for atoms in dihAnglesList:
            originalAngle = GetDihedralDeg(conf, atoms[0], atoms[1], atoms[2], atoms[3])
            originalAngles.append(originalAngle)

        return originalAngles

    def setRightAngles(conf, originalAngles):
        for atoms, angle in zip(dihAnglesList, originalAngles):
            SetDihedralDeg(conf, atoms[0], atoms[1], atoms[2], atoms[3], angle)

    molH = Chem.AddHs(mol)

    dihAnglesList = findAngles(mol)

    originalAngles = findAnglesDeg(molH.GetConformer(), dihAnglesList)

    AllChem.EmbedMolecule(molH)
    setRightAngles(molH.GetConformer(), originalAngles)

    MMFFOptimizeMolecule(molH, maxIters=500, mmffVariant="MMFF94s")

    return mol

def doMagic(directory):
    mols = readSdf(fr"{directory}")

    i = 0
    if not os.path.exists(fr"{directory}\mol"):
        os.makedirs(fr"{directory}\mol")

    for mol in mols:
        MolToMolFile(mol, fr"{directory}\mol\{i}.mol")
        i += 1

    fileNames = next(os.walk(fr"{directory}\mol"), (None, None, []))[2]

    if not os.path.exists(fr"{directory}\gjf"):
        os.makedirs(fr"{directory}\gjf")

    for fileName in fileNames:
        mol = MolFromMolFile(fr"{directory}\mol\{fileName}")

        dihAnglesList = findAngles(mol)
        mol = processMol(mol, dihAnglesList)
        MolToMolFile(mol, fr"{directory}\gjf\{i}.mol")

        MolToGjfFile(mol, fileName[:-4], fr"{directory}\gjf")

#directory = r"C:\Users\Leo\Desktop\test"
#mols = readSdf(directory)
#
#
#i = 0
#for mol in mols:
#    dihAnglesList = findAngles(mol)
#    mol = processMol(mol, dihAnglesList)
#    MolToMolFile(mol, fr"{directory}\gjf\{i}.mol")
#    i += 1

#directory = r"C:\Users\Leo\Desktop\test"
#mols = readSdf(fr"{directory}")
#
#i = 0
#if not os.path.exists(fr"{directory}\mol"):
#    os.makedirs(fr"{directory}\mol")
#
#for mol in mols:
#    MolToMolFile(mol, fr"{directory}\mol\{i}.mol")
#    i += 1
#
#fileNames = next(os.walk(fr"{directory}\mol"), (None, None, []))[2]
#
#if not os.path.exists(fr"{directory}\gjf"):
#    os.makedirs(fr"{directory}\gjf")
#
#for fileName in fileNames:
#    mol = MolFromMolFile(fr"{directory}\mol\{fileName}")
#
#    dihAnglesList = findAngles(mol)
#    mol = processMol(mol, dihAnglesList)
#    MolToMolFile(mol, fr"{directory}\gjf\{i}.mol")
#
#    MolToGjfFile(mol, fileName[:-4], fr"{directory}\gjf")

doMagic(r"C:\Users\Leo\Desktop\test")