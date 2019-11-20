import numpy as np

class Flow:
    def __init__(X,Y,PlateY,L):
        self.L = L
        #length of the surface in x direction
        self.X = X
        #length of the surface in y direction
        self.Y = Y
        self.V0 = 1
        self.PlateXFront = .25
        self.PlateXBack = .375
        self.PlateYTop = .25
        self.plateFront = self.PlateXFront/self.X
        self.plateBack = self.PlateXBack/self.X
        self.plateTop = self.PlateYTop/self.Y
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
        self.Psi = np.zeros((L+1,L+1))
        for j in range(0, self.L+1):
            freeflow = self.V0*self.X*j/(self.L+1)
            for l in range(0, self.L+1):
                if(inPlate(l,j)):
                    self.Psi[l][j] = 0
                else:
                    self.Psi[l][j] = freeflow


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
        for j in range(0, self.L+1): #y-ccordinate j/L
            for l in range(0, self.L+1): #x-coordinate l/L
                if(l == 0): #Upstream
                    self.psi[l][j] = self.psi[l+1][j]
                elif(l == L): #Downstream
                    self.psi[l][j] = self.psi[l-1][j]
                elif(j == L): #Side opposite plate
                    self.psi[l][j] = self.psi[l][j-1]*self.V0/self.h
                elif(j == 0): #Plate side
                    self.psi[l][j] = 0
                elif((l/(L+1) == self.plateFront) and (j/(L+1) < self.plateTop)): #Front of plate
                    self.psi[l][j] = -1 #TODO
                elif((l/(L+1) == self.plateBack) and (j/(L+1) < self.plateTop)): #Back of plate
                    self.psi[l][j] = -1 #TODO
                elif((j/(L+1) == self.plateTop) and (l/(L+1) >= self.plateFront) and (l/(L+1) <= self.plateBack)): #Top of plate
                    self.psi[l][j] = -1 #TODO
                elif(self.inplate(l,j)): #Inside plate
                    self.psi[l][j] = 0


    def updateWInternal(self):
        for l in range(1,L+1):
            for j in range(1,L+1):
                if not self.inPlate(l,j):
                    psiStencil = self.PsiStencil(l,j)
                    if w != False:
                        self.W[l,j] = psiStencil
                else:
                    self.W[l,j] = 0

    def updateWBoundary(self):
        for j in range(0, self.L+1): #y-ccordinate j/L
            for l in range(0, self.L+1): #x-coordinate l/L
                if(l == 0): #Upstream
                    self.w[l][j] = 0
                elif(l == L): #Downstream
                    self.w[l][j] = self.w[l-1][j]
                elif(j == L): #Side opposite plate
                    self.w[l][j] = 0
                elif(j == 0): #Plate side
                    self.w[l][j] = 0
                elif((l/(L+1) == self.plateFront) and (j/(L+1) < self.plateTop)): #Front of plate
                    self.w[l][j] = -2*(self.h*self.h)*self.psi[l][j-1]
                elif((l/(L+1) == self.plateBack) and (j/(L+1) < self.plateTop)): #Back of plate
                    self.w[l][j] = -2/(self.h*self.h)*self.psi[l][j+1]
                elif((j/(L+1) == self.plateTop) and (l/(L+1) >= self.plateFront) and (l/(L+1) <= self.plateBack)): #Top of plate
                    self.w[l][j] = -2/(self.h*self.h)*self.psi[l+1][j]
                elif(self.inplate(l,j): #Inside plate
                    self.w[l][j] = 0

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
