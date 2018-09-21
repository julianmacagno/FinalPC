import subprocess
lista = [] 
prom = 0.0

for _ in range (0,300):
    p = subprocess.Popen("./bin/finalPC", stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    lista.append(float(output))
    prom += float(output)

prom = prom / len(lista)

print lista
print "\nPromedio de 100 iteraciones en secuencial: " + str(prom) + " milisegundos\n\n"

prom = 0.0
lista = []

for _ in range (0,300):
    p = subprocess.Popen("./bin/finalPC-OpenMP", stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    lista.append(float(output))
    prom += float(output)

prom = prom / len(lista)

print lista
print "\nPromedio de 100 iteraciones en paralelo con OpenMP: " + str(prom) + " milisegundos"