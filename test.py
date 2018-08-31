import numpy as np

t = np.linspace(1, 10, 10)
tau = t

T, Tau = np.meshgrid(t, tau)

print "\n\n"
print t
print "\n\n"
print tau
print "\n\n"
print T
print "\n\n"
print Tau
print "\n\n"
print T.size
print "\n\n"




for i in range(0, T.size):
    print T[i]