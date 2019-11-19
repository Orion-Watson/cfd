import numpy as np

class Flow:
    def __init__(X,Y,L):
        self.L = L
        #length of the surface in x direction
        self.X = X
        #length of the surface in y direction
        self.Y = Y
        self.W = np.zeros((L+1,L+1))
        self.Psi = np.zeros((L+1,L+1))
        #indicates wheather Psi or W should be updated
        self.updatePsi = True
        self.solutionFound = False

    def initialize(self):
        self.freeFlowInit()

    def freeFlowInit(self):
        self.W = np.zeros((L+1,L+1))
        self.Psi = np.zeros((L+1,L+1))

    def solve(self):
        #initialize Psi matrix and W matrix
        self.initialize()
        while not self.solutionFound:
            if self.updatePsi:
                self.updatePsi()
            else:
                self.updateW()
            self.updatePsi = !self.updatePsi

    def updatePsi(self):
        self.updatePsiInternal()
        self.updatePsiBoundary()

    def updateW(self):
        self.updateWInternal()
        self.updateWBoundary()

    def updatePsiInternal(self):
        print()

    def updatePsiBoundary(self):
        print()

    def updateWInternal(self):
        for x in range(1,L+1):
            for y in range(1,L+1):
                psiStencil = self.PsiStencil(x,y)
                if w != False:
                    self.W[x,y] = psiStencil

    def updateWBoundary(self):
        print()

    def PsiStencil(self,x,y):
        if x != 0 and y != 0 and x !=  (self.L +1) and y != (self.L +1):
            Psi = self.Psi[x,y]
            stencil = -4*Psi + self.Psi[x+1,y] + self.Psi[x-1,y] + self.Psi[x,y+1]+ self.Psi[x,y-1]
            return stencil
        else:
            return False

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


    """I don't think we need any super indexing anymore """
    # """Takes a superindexed vector and returns a normally indexed matrix"""
    # def desuperindex(super_index, L):
    #     regular_index = np.zeros((L+1, L+1))
    #     for j in range(0, L+1):
    #         for l in range(0, L+1):
    #             i = (j*(L+1))+l
    #             regular_index[l][j] = super_index[i]
    #     return regular_index

L = 16
a, b = make_superindexed(L)
print("a: ", a)
print("b: ", b)
print(desuperindex(b, L))
