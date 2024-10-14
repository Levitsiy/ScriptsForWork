import os
import numpy as np
from rdkit import Chem

from ScriptsForWork.get_out_data.xyz2mol import xyz2mol

#directory = fr"{os.getcwd()}"
directory = fr"C:\Users\Leo\Desktop\Работа НИОХ\сalcilations\annulation library\3 rings\08.10.24 готовые расчеты"

fileNames = next(os.walk(directory), (None, None, []))[2]
newFile = open(fr"{directory}\moleculesInfo.csv", "a")

for fileName in fileNames:
    if fileName[-3:] != "out":
        continue

    iOLineNumberStart = 0  # iO = Input orientation
    iOLineNumberFinish = 0  # iO = Input orientation

    HOMOLineNumber = 0
    LUMOLineNumber = 0

    periodicTable = {"1":"H",
                     "6":"C",
                     "7":"N",
                     "8":"O",
                     "9":"F",
                     "15":"P",
                     "16":"S"}

    with open(fr"{directory}\{fileName}", "r") as file:
        linesArray = file.readlines()
        pass

    if linesArray[len(linesArray)-1].find("Normal termination of Gaussian") == -1:
        os.rename(fr"{directory}\{fileName}", fr"{directory}\{fileName[:-4]}_fail.out")
        print(fr"{fileName} fail termination")
        continue
    else: print(fr"{fileName} normal termination")

    for lineNumber in range(len(linesArray) - 1, 0, -1):

        if "Distance matrix (angstroms):" in linesArray[lineNumber]:
            iOLineNumberFinish = lineNumber - 1

        if "Input orientation:" in linesArray[lineNumber]:
            iOLineNumberStart = lineNumber + 5
            break

    iOArrayCoordinates = np.zeros((iOLineNumberFinish - iOLineNumberStart, 3))  # iO = Input orientation
    iOArrayAtomTypes = []  # iO = Input orientation
    y = 0

    for i in range(iOLineNumberStart, iOLineNumberFinish): #Заполняет массив с данными о геомертии
        string = linesArray[i].split()
        iOArrayAtomTypes.append(int(string[1]))

        iOArrayCoordinates[y][0] = string[3]
        iOArrayCoordinates[y][1] = string[4]
        iOArrayCoordinates[y][2] = string[5]

        y += 1


    for lineNumber in range(len(linesArray) - 1, 0, -1):

        if " Alpha virt. eigenvalues" in linesArray[lineNumber]:
            LUMOLineNumber = lineNumber

        if "Alpha  occ. eigenvalues" in linesArray[lineNumber]:
            HOMOLineNumber = lineNumber
            break

    HOLO_energy = float(linesArray[HOMOLineNumber].split()[-1]) * 27.2114
    LUMO_energy = float(linesArray[LUMOLineNumber].split()[4]) * 27.2114

    mol = xyz2mol(atoms = iOArrayAtomTypes, coordinates = iOArrayCoordinates)[0]
    smiles = Chem.MolToSmiles(mol)

    #stringToNewFile = f"{fileName[0:-4]} {smiles} {HOLO_energy} {LUMO_energy}\n"
    stringToNewFile = f"{fileName[0:-4]} {smiles} {LUMO_energy - HOLO_energy}\n"
    newFile.write(stringToNewFile)

newFile.close()

