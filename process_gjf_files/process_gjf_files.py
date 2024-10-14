import os
from os import walk

#directory = fr"{os.getcwd()}"
directory = fr"C:\Users\Leo\Desktop\Работа НИОХ\сalcilations\annulation library\3 rings\08.10.24 аннелирование на расчеты"

fileNames = next(walk(directory), (None, None, []))[2]

for fileName in fileNames:
    with open(fr"{directory}\{fileName}", "r") as file:
        linesArray = file.readlines()
        shortArray = linesArray[5:]
        finalArray = []
        pass

    os.rename(fr"{directory}\{fileName}", fr"{directory}\{fileName[:-4]}_raw.gjf") #Помечает рабочий файл флагом _raw

    for s in shortArray:
        if s != "\n" and s[0] != "L":
            split = s.split()
            newS = fr"{split[0]} {split[2]} {split[3]} {split[4]}"
            finalArray.append(newS)

    for i in range(4, 0, -1):
        finalArray.insert(0, linesArray[i].rstrip())

    with open(fr"{directory}\{fileName[:-4]}_processed.gjf", "a") as newFile: #Создает новый файл с измененными удаленным слолбцом нулей и помечает флагом _processed
        newFile.write(f"%nprocshare = 12\n"  # 12 for G7, 8 for G6
                      f"%mem = 24Gb\n" # 24 for G7, 16 for G6
                      f"%chk={fileName[:-4]}.chk\n"
                      f"#p opt=vtight b3lyp/6-311++g(d,p) nosymm int=superfinegrid scf=verytight\n")
        for s in finalArray:
            newS = f"{s}\n"
            newFile.write(newS)
        newFile.write(f"\n\n")
        pass

