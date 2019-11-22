import numpy as np
import matplotlib.pyplot as plt

class Flow:
    def __init__(self, X, Y, L):
        self.L = L #highest index in matrix
        self.X = X #length of surface in x direction
        self.Y = Y #length of the surface in y direction
        self.V0 = 1 #Velocity of incoming fluid
        #Plate boundaries
        self.PlateXFront = .25
        self.PlateXBack = .375
        self.PlateYTop = .25
        #Plate boundaries as fractions
        self.plateFront = self.PlateXFront/self.X
        self.plateBack = self.PlateXBack/self.X
        self.plateTop = self.PlateYTop/self.Y
        self.W = np.zeros((L+1,L+1)) #Vorticity matrix
        self.Psi = np.zeros((L+1,L+1)) #Stream function matrix
        self.hX = X/(L+1) #step sizes
        self.hY = Y/(L+1)
        self.updatePsiMatrix = True #Whether w has been updates since last psi update
        self.solutionFound = False #True when within tolerance
        self.maxIterations = 10 #Max iterations before giving up

    """Returns true if (l,j) is inside the plate"""
    def inPlate(self, l, j):
        front = l/(self.L+1) > self.plateFront
        back = l/(self.L+1) < self.plateBack
        top = j/(self.L+1) < self.plateTop
        return (front and back and top)

    """Sets reasonable initial values"""
    def initialize(self):
        self.freeFlowInit()

    """Sets all values to free flow conditions"""
    def freeFlowInit(self):
        #initialize all vorticities as 0
        self.W = np.zeros((self.L+1,self.L+1))

        #set psi to free flow conditions except for 0 in plate
        self.Psi = np.zeros((self.L+1,self.L+1))
        for j in range(0, self.L+1):
            freeflow = self.V0*self.X*j/(self.L+1)
            for l in range(0, self.L+1):
                if(self.inPlate(l, j)):
                    self.Psi[l][j] = 0 #Inside Plate
                else:
                    self.Psi[l][j] = freeflow

    """Alternate between updating psi matrix and w matrix"""
    def solve(self):
        #initialize Psi matrix and W matrix
        self.initialize()
        print("init: ")
        print("W: ")
        print(self.W)
        print("Psi: ")
        print(self.Psi)
        i = 0
        while not self.solutionFound and i < self.maxIterations:
            if self.updatePsiMatrix:
                self.updatePsi()
                print("updatePsi: ")
                print("W: ")
                print(self.W)
                print("Psi: ")
                print(self.Psi)
            else:
                self.updateW()
                print("updateW: ")
                print("W: ")
                print(self.W)
                print("Psi: ")
                print(self.Psi)
            self.updatePsiMatrix = not self.updatePsiMatrix #Change which we update next time
            i +=1

    def updatePsi(self):
        self.updatePsiInternal()
        self.updatePsiBoundary()

    def updateW(self):
        self.updateWInternal()
        self.updateWBoundary()

    """TODO: fix this, gives crazy noise instead of real values"""
    def updatePsiInternal(self):
        for j in range(0, self.L+1):
            for l in range(0, self.L+1):
                if l != 0 and j != 0 and l != (self.L) and j != (self.L):
                    if not self.inPlate(l,j):
                        Psi = self.PsiStencilLaPlaz(l,j) + self.W[j][l] #TODO: fix this
                        self.Psi[j,l] = Psi

    """Imposes boundary conditions, call after updating psi internal"""
    def updatePsiBoundary(self):
        for j in range(0, self.L+1): #y-ccordinate j/L
            for l in range(0, self.L+1): #x-coordinate l/L
                if not self.inPlate(l,j):
                    if(j == self.L): #Side opposite plate
                        self.Psi[l][j] = self.V0*self.Y
                    elif(l == 0): #Upstream
                        self.Psi[l][j] = self.Psi[l+1][j]
                    elif(l == self.L): #Downstream
                        self.Psi[l][j] = self.Psi[l-1][j]
                    elif(j == 0): #Plate side
                        self.Psi[l][j] = 0
                    elif((l/(self.L+1) == self.plateFront) and (j/(self.L+1) < self.plateTop)): #Front of plate
                        self.Psi[l][j] = 0
                    elif((l/(self.L+1) == self.plateBack) and (j/(self.L+1) < self.plateTop)): #Back of plate
                        self.Psi[l][j] = 0
                    elif((j/(self.L+1) == self.plateTop) and (l/(self.L+1) >= self.plateFront) and (l/(self.L+1) <= self.plateBack)): #Top of plate
                        self.Psi[l][j] = 0
                    elif(self.inPlate(l,j)): #Inside plate
                        self.Psi[l][j] = 0
                else:
                    self.W[l,j] = 0

    def updateWInternal(self):
        for l in range(1,self.L+1):
            for j in range(1,self.L+1):
                if not self.inPlate(l,j):
                    if l != 0 and j != 0 and l != (self.L) and j != (self.L):
                        PsiStencilLaPlaz = self.PsiStencilLaPlaz(l,j)
                        self.W[l,j] = PsiStencilLaPlaz
                else:
                    self.W[l,j] = 0

    def updateWBoundary(self):
        for j in range(0, self.L+1): #y-ccordinate j/L
            for l in range(0, self.L+1): #x-coordinate l/L
                if not self.inPlate(l,j):
                    if(j == self.L): #Side opposite plate
                        self.W[l][j] = 0
                    elif(l == 0): #Upstream
                        self.W[l][j] = 0
                    elif(l == self.L): #Downstream
                        self.W[l][j] = self.W[l-1][j]
                    elif(j == 0): #Plate side
                        self.W[l][j] = 0
                    elif((l/(self.L+1) == self.plateFront) and (j/(self.L+1) < self.plateTop)): #Front of plate
                        self.W[l][j] = (-2/(self.hY*self.hY))*self.Psi[l][j-1]
                    elif((l/(self.L+1) == self.plateBack) and (j/(self.L+1) < self.plateTop)): #Back of plate
                        self.W[l][j] = (-2/(self.hY*self.hY))*self.Psi[l][j+1]
                    elif((j/(self.L+1) == self.plateTop) and (l/(self.L+1) >= self.plateFront) and (l/(self.L+1) <= self.plateBack)): #Top of plate
                        self.W[l][j] = (-2/(self.hX*self.hX))*self.Psi[l+1][j]
                    elif(self.inPlate(l,j)): #Inside plate
                        self.W[l][j] = 0
                else:
                    self.W[l,j] = 0

    """Stencil for Laplacian of psi"""
    def PsiStencilLaPlaz(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            Psi = self.Psi[l,j]
            stencil = -Psi + (1/4)*(self.Psi[l+1,j] + self.Psi[l-1,j] + self.Psi[l,j+1]+ self.Psi[l,j-1])
            return stencil
        else:
            return False

    def PsiStencilLaPlazDy(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            h = j/(L+1) * self.Y
            stencil = (self.Psi[l,j+1] - self.Psi[l,j-1])/h
            return stencil

    def PsiStencilLaPlazDx(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            h = l/(L+1) * self.X
            stencil = (self.Psi[l+1,j] - self.Psi[l-1,j]) / h
            return stencil
        else:
            return False

    """Compute X and Y coordinates"""
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

    """Graph contour lines of psi"""
    def flowGraph(self):
        print("graph: ")
        print("W: ")
        print(self.W)
        print("Psi: ")
        print(self.Psi)
        #conlines = np.linspace(0, 1, 20)
        X,Y = self.getXYCoords(self.Psi)
        plt.contour(X,Y,self.Psi, cmap = 'plasma')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title(r"$\Psi$")
        plt.show()


free = Flow(1, 1, 7)
free.initialize()
#free.updateWBoundary()
#free.updatePsiBoundary()

free.solve()
free.flowGraph()
