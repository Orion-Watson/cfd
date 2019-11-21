import numpy as np
import matplotlib.pyplot as plt

class Flow:
    def __init__(self, X, Y, L):
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
        self.h = self.X/self.Y
        #indicates wheather Psi or W should be updated
        self.updatePsi = True
        self.solutionFound = False

    def inPlate(self, l, j):
        front = l/(self.L+1) > self.plateFront
        back = l/(self.L+1) < self.plateBack
        top = j/(self.L+1) < self.plateTop
        return (front and back and top)

    def initialize(self):
        self.freeFlowInit()

    def freeFlowInit(self):
        #initialize all vorticities as 0
        self.W = np.zeros((self.L+1,self.L+1))
        self.Psi = np.zeros((self.L+1,self.L+1))
        for j in range(0, self.L+1):
            freeflow = self.V0*self.X*j/(self.L+1)
            for l in range(0, self.L+1):
                if(self.inPlate(l, j)):
                    self.Psi[l][j] = 0 #Inside Plate
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
            self.updatePsi = not self.updatePsi

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
                if(j == self.L): #Side opposite plate
                    self.Psi[l][j] = self.Psi[l][j-1]*self.V0/self.h
                elif(l == 0): #Upstream
                    self.Psi[l][j] = self.Psi[l+1][j]
                elif(l == self.L): #Downstream
                    self.Psi[l][j] = self.Psi[l-1][j]
                elif(j == 0): #Plate side
                    self.Psi[l][j] = 0
                elif((l/(self.L+1) == self.plateFront) and (j/(self.L+1) < self.plateTop)): #Front of plate
                    self.Psi[l][j] = -1 #TODO
                elif((l/(self.L+1) == self.plateBack) and (j/(self.L+1) < self.plateTop)): #Back of plate
                    self.Psi[l][j] = -1 #TODO
                elif((j/(self.L+1) == self.plateTop) and (l/(self.L+1) >= self.plateFront) and (l/(self.L+1) <= self.plateBack)): #Top of plate
                    self.Psi[l][j] = -1 #TODO
                elif(self.inPlate(l,j)): #Inside plate
                    self.Psi[l][j] = 0


    def updateWInternal(self):
        for l in range(1,self.L+1):
            for j in range(1,self.L+1):
                if not self.inPlate(l,j):
                    psiStencil = self.PsiStencil(l,j)
                    if w != False:
                        self.W[l,j] = psiStencil
                else:
                    self.W[l,j] = 0

    def updateWBoundary(self):
        for j in range(0, self.L+1): #y-ccordinate j/L
            for l in range(0, self.L+1): #x-coordinate l/L
                if(j == self.L): #Side opposite plate
                    self.W[l][j] = 0
                elif(l == 0): #Upstream
                    self.W[l][j] = 0
                elif(l == self.L): #Downstream
                    self.W[l][j] = self.W[l-1][j]
                elif(j == 0): #Plate side
                    self.W[l][j] = 0
                elif((l/(self.L+1) == self.plateFront) and (j/(self.L+1) < self.plateTop)): #Front of plate
                    self.W[l][j] = -2*(self.h*self.h)*self.psi[l][j-1]
                elif((l/(self.L+1) == self.plateBack) and (j/(self.L+1) < self.plateTop)): #Back of plate
                    self.W[l][j] = -2/(self.h*self.h)*self.psi[l][j+1]
                elif((j/(self.L+1) == self.plateTop) and (l/(self.L+1) >= self.plateFront) and (l/(self.L+1) <= self.plateBack)): #Top of plate
                    self.W[l][j] = -2/(self.h*self.h)*self.psi[l+1][j]
                elif(self.inPlate(l,j)): #Inside plate
                    self.W[l][j] = 0

    def PsiStencil(self,l,j):
        if l != 0 and j != 0 and l !=  (self.self.L +1) and j != (self.self.L +1):
            Psi = self.Psi[l,j]
            stencil = -4*Psi + self.Psi[l+1,j] + self.Psi[l-1,j] + self.Psi[l,j+1]+ self.Psi[l,j-1]
            return stencil
        else:
            return False

    def PsiStencilDy(self,l,j):
        if l != 0 and j != 0 and l !=  (self.self.L +1) and j != (self.self.L +1):
            h = j/(L+1) * self.Y
            stencil = (self.Psi[l,j+1] - self.Psi[l,j-1])/h
            return stencil
        else:
            return False

    def PsiStencilDx(self,l,j):
        if l != 0 and j != 0 and l !=  (self.self.L +1) and j != (self.self.L +1):
            h = l/(L+1) * self.X
            stencil = (self.Psi[l+1,j] - self.Psi[l-1,j]) / h
            return stencil
        else:
            return False

    def getXYCoords(self,A):
        cordsX = np.zeros((len(A),len(A)))
        cordsY = np.zeros((len(A),len(A)))
        for l in range(len(A)):
            for j in range(len(A)):
                x = (l/self.L)*self.X
                y = (j/self.L)*self.Y
                cordsX[l,j] = x
                cordsY[l,j] = y

        return cordsX, cordsY

    def flowGraph(self):
        conlines = np.linspace(0, 1, 20)
        X,Y = self.getXYCoords(self.Psi)
        plt.contour(X,Y,self.Psi, conlines, cmap = 'plasma')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title(r"$\Psi$")
        plt.show()


free = Flow(1, 1, 32)
free.initialize()
free.updateWBoundary()
free.updatePsiBoundary()
free.flowGraph()
