import numpy as np

class Flow:
    def __init__(X,Y,PlateY,L):
        self.L = L
        #length of the surface in x direction
        self.X = X
        #length of the surface in y direction
        self.PlateXFront = .25
        self.PlateXBack = .375
        self.PlateYTop = .25
        self.Y = Y
        self.W = np.zeros((L+1,L+1))
        self.Psi = np.zeros((L+1,L+1))
        #indicates wheather Psi or W should be updated
        self.updatePsi = True
        self.solutionFound = False

    def initialize(self):
        self.freeFlowInit()

    def freeFlowInit(self):
        #initialize all vorticities as 0
        self.W = np.zeros((L+1,L+1))
        """TODO: need to figure out how we should initialize PSI """
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
        for l in range(1,L+1):
            for j in range(1,L+1):
                if self.inPlate(l,j):
                    psiStencil = self.PsiStencil(l,j)
                    if w != False:
                        self.W[l,j] = psiStencil
                else:
                    self.W[l,j] = 0

    def updateWBoundary(self):
        print()

    def PsiStencil(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            Psi = self.Psi[l,j]
            stencil = -4*Psi + self.Psi[l+1,j] + self.Psi[l-1,j] + self.Psi[l,j+1]+ self.Psi[l,j-1]
            return stencil
        else:
            return False

    def inPlate(self,l,j):
        normalizedDistanceFromFront = self.PlateXFront / self.X
        normalizedDistanceOfBack = (self.X - self.PlateXBack)/ self.X
        if normalizedDistanceFromFront > (l / self.L) and normalizedDistanceOfBack > (l / self.L):
            normalizedDistanceOfTop = (self.Y - self.PlateYTop)/ self.Y
            if (j / self.J) < normalizedDistanceOfTop:
                return True
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
