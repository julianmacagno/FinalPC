cFile = open("cFile.txt", "r")
pyFile = open("pyFile.txt", "r")

cList = cFile.readlines()
pyList = pyFile.readlines()

counter = 0

for i in range (0, len(cList)) :
    string =  str(pyList[i])
    if(string.count("e") != 0):
        string = string[:7] + string[string.index("e"):]
    

    # print string
    res = float(cList[i]) - float(string)
    if(res != 0) :
        print "pyList: " + string.replace("\n", '')
        print "cList: " + str(cList[i]).replace("\n", '')
        print "Res: " + str(res) + "\n"
        counter += 1

print "\n\n"
print "Cantidad de resultados distintos: " + str(counter)
cFile.close()
pyFile.close()