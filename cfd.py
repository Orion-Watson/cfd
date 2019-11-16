import numpy as np


def make_superindexed(L): #Only use L as multiple of 8
    if((L%8) != 0):
        raise Exception("L not a nultiple of 8!")
    n = (L+1)**2
    A = np.zeros((n, n))
    b = np.zeros(n)
    for j in range(0, L+1): #y-ccordinate j/L
        for l in range(0, L+1): #x-coordinate l/L
            i = (j*(L+1))+l
            if(l == 0): #Upstream
                A[i][i] = 1
                b[i] = 1
            elif(l == L): #Downstream
                A[i][i] = 1
                b[i] = 3
            elif(j == L): #Side opposite plate
                A[i][i] = 1
                b[i] = 2
            elif(j == 0): #Plate side
                A[i][i] = 1
                b[i] = 4
            elif((l == L/4) and (j < L/4)): #Front of plate
                A[i][i] = 1
                b[i] = 7
            elif((l == 3*L/8) and (j < L/4)): #Back of plate
                A[i][i] = 1
                b[i] = 5
            elif((j == L/4) and (l >= L/4) and (l <= 3*L/8)): #Top of plate
                A[i][i] = 1
                b[i] = 6
            elif((j < L/4) and (l > L/4) and (l < 3*L/8)): #Inside plate
                A[i][i] = 1
                b[i] = 8
            else:
                A[i][i] = -1
                b[i] = 0

    return A, b


"""Takes a superindexed vector and returns a normally indexed matrix"""
def desuperindex(super_index, L):
    regular_index = np.zeros((L+1, L+1))
    for j in range(0, L+1):
        for l in range(0, L+1):
            i = (j*(L+1))+l
            regular_index[l][j] = super_index[i]
    return regular_index

L = 16
a, b = make_superindexed(L)
print(desuperindex(b, L))
