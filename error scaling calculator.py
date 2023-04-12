import numpy as np
k=3
e=1e-3
for m in range(7):
    t = m+k
    f= 1 - 8*e*pow(2,k)*pow(t,2)
    print(f)
