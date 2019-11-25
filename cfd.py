import numpy as np
import matplotlib.pyplot as plt

class Flow:
    def __init__(self, X, Y, L):
        self.L = L #highest index in matrix
        self.X = X #length of surface in x direction
        self.Y = Y #length of the surface in y direction
        self.V0 = 1 #Velocity of incoming fluid
        self.w = .3 #over relaxation factor
        self.viscosity = .8 #viscosity of the fluid
        #Plate boundaries
        self.PlateXFront = .25
        self.PlateXBack = .375
        self.PlateYTop = .25
        #Plate boundaries as fractions
        self.plateFront = self.PlateXFront/self.X
        self.plateBack = self.PlateXBack/self.X
        self.plateTop = self.PlateYTop/self.Y
        #Vorticity matrix
        self.W = np.zeros((L+1,L+1))
        #Stream function matrix
        self.Psi = np.zeros((L+1,L+1))
        #step sizes
        self.hX = X/(L+1)
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
            freeflow = self.V0*self.hX*j
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

    """ Updates Psi matrix by calling functions to update Psi internal points
    and Psi boundary points"""
    def updatePsi(self):
        self.updatePsiInternal()
        self.updatePsiBoundary()

    """ Updates W matrix by calling functions to update Psi internal points
    and Psi boundary points"""
    def updateW(self):
        self.updateWInternal()
        self.updateWBoundary()

    """Updates the internal values of the Psi Matrix based on previous Psi
    values and updated w values"""
    def updatePsiInternal(self):
        for j in range(0, self.L+1):
            for l in range(0, self.L+1):
                #Do not update boundary Psi values
                if l != 0 and j != 0 and l != (self.L) and j != (self.L):
                    #do not update Psi values in the plate
                    if not self.inPlate(l,j):
                        W = self.W[l,j]
                        PsiSquareStencil = self.PsiSquareStencil(l,j)
                        Psi = ((1-self.w) * self.Psi[l,j]) + (self.w * PsiSquareStencil) - W
                        self.Psi[j,l] = Psi
                    else:
                        #set values inside plate to 0
                        self.Psi[j,l] = 0

    """Imposes boundary conditions on the Psi matrix"""
    def updatePsiBoundary(self):
        for j in range(0, self.L+1): #y-ccordinate j/L
            for l in range(0, self.L+1): #x-coordinate l/L
                if not self.inPlate(l,j):
                    #Side opposite plate
                    if(j == self.L):
                        self.Psi[l][j] = self.V0*self.Y
                    #Upstream
                    elif(l == 0):
                        self.Psi[l][j] = self.Psi[l+1][j]
                    #Downstream
                    elif(l == self.L):
                        self.Psi[l][j] = self.Psi[l-1][j]
                    #Plate side
                    elif(j == 0):
                        self.Psi[l][j] = 0
                    #Front of plate
                    elif((l/(self.L+1) == self.plateFront) and (j/(self.L+1) < self.plateTop)):
                        self.Psi[l][j] = 0
                    #Back of plate
                    elif((l/(self.L+1) == self.plateBack) and (j/(self.L+1) < self.plateTop)):
                        self.Psi[l][j] = 0
                    #Top of plate
                    elif((j/(self.L+1) == self.plateTop) and (l/(self.L+1) >= self.plateFront) and (l/(self.L+1) <= self.plateBack)):
                        self.Psi[l][j] = 0
                    #Inside plate
                    elif(self.inPlate(l,j)):
                        self.Psi[l][j] = 0

    """Updates the internal values of the w Matrix based on previous w
    values and updated Psi values"""
    def updateWInternal(self):
        for l in range(1,self.L+1):
            for j in range(1,self.L+1):
                #only update internal points
                if l != 0 and j != 0 and l != (self.L) and j != (self.L):
                    #make sure index is not in plate
                    if not self.inPlate(l,j):
                        partial = 1/(self.viscosity) * (self.PsiStencilDy(l,j) * self.WStencilDx(l,j) -
                        (self.PsiStencilDx(l,j) * self.WStencilDy(l,j)))
                        squareStencil = self.WSquareStencil(l,j)
                        W = ((1-self.w)*self.W[l,j]) + (self.w)*squareStencil + partial
                        self.W[l,j] = W
                else:
                    #set w in plate to 0
                    self.W[l,j] = 0

    """Imposes boundary conditions on the W matrix"""
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
            PsiSquareStencil = self.PsiSquareStencil(l,j)
            stencil = -Psi + PsiSquareStencil
            return stencil

    def PsiSquareStencil(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (1/4)*(self.Psi[l+1,j] + self.Psi[l-1,j] + self.Psi[l,j+1]+ self.Psi[l,j-1])
            return stencil

    def WSquareStencil(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (1/4)*(self.W[l+1,j] + self.W[l-1,j] + self.W[l,j+1]+ self.W[l,j-1])
            return stencil

    def PsiStencilDy(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (self.Psi[l,j+1] - self.Psi[l,j-1])/(2*self.hY)
            return stencil

    def PsiStencilDx(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (self.Psi[l+1,j] - self.Psi[l-1,j]) /(2*self.hX )
            return stencil
        else:
            return False

    def WStencilDy(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (self.W[l,j+1] - self.W[l,j-1])/(2*self.hY)
            return stencil

    def WStencilDx(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (self.W[l+1,j] - self.W[l-1,j])/(2*self.hX)
            return stencil

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
