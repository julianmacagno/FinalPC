cFile = open("cFile.txt", "r")
pyFile = open("pyFile.txt", "r")
cList = cFile.readlines()
pyList = pyFile.readlines()
counter = 0
precision = 5

for i in range (0, len(pyList)) :
    py_string =  str(pyList[i])
    if(py_string.count("e") != 0):
        py_string = py_string[:py_string.index(".") + precision] + py_string[py_string.index("e"):]
    else:
        py_string = py_string[:py_string.index(".") + precision]

    c_string =  str(cList[i])
    if(c_string.count("e") != 0):
        length = len(c_string[c_string.index("."):c_string.index("e")])
        if(length>precision):
            length = precision
        c_string = c_string[:c_string.index(".") + length] + c_string[c_string.index("e"):]
    else:
        if(c_string.count(".") != 0):
            c_string = c_string[:c_string.index(".") + precision]

    res = abs(float(c_string) - float(py_string))
    if(res != 0) :
        print "pyList:   " + py_string.replace("\n", '')
        print "cList:    " + str(c_string).replace("\n", '')
        print "Res:      " + str(res) + "\n"
        counter += 1

print "\n\n"
print "Cantidad de resultados distintos: " + str(counter)
cFile.close()
pyFile.close()