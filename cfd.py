import numpy as np
import matplotlib.pyplot as plt

class Flow:
    def __init__(self, X, Y, L):
        self.L = L #highest index in matrix
        self.X = X #length of surface in x direction
        self.Y = Y #length of the surface in y direction
        self.V0 = 1 #Velocity of incoming fluid
        self.w = 1.5 #over relaxation factor
        self.viscosity = .1 #viscosity of the fluid
        #Plate boundaries
        self.PlateXFront = .2 * X
        self.PlateXBack = .5 * X
        self.PlateYTop = .2 * Y
        #Vorticity matrix
        self.W = np.zeros((L+1,L+1))
        #Stream function matrix
        self.Psi = np.zeros((L+1,L+1))
        #step sizes
        self.hX = X/(L)
        self.hY = Y/(L)
        self.updatePsiMatrix = True #Whether w has been updates since last psi update
        self.solutionFound = False #True when within tolerance
        self.maxIterations = 7 #Max iterations before giving up

    """Returns true if (l,j) is inside the plate"""
    def inPlate(self, l, j):
        front = l*self.hX > self.PlateXFront
        back = l*self.hX < self.PlateXBack
        top = j*self.hY < self.PlateYTop
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
            freeflow = self.V0*self.hY*j
            for l in range(0, self.L+1):
                if(self.inPlate(l, j)):
                    self.Psi[l][j] = 0 #Inside Plate
                else:
                    self.Psi[l][j] = freeflow

    """Alternate between updating psi matrix and w matrix"""
    def solve(self):
        #initialize Psi matrix and W matrix
        self.initialize()
        self.flowGraph()
        i = 0
        while not self.solutionFound and i < self.maxIterations:
            if self.updatePsiMatrix:
                print("updatePsi")
                self.updatePsi()
                print("Psi: ", self.Psi)
                self.flowGraph()
            else:
                print("updateW")
                self.updateW()
                print("W: ", self.W)
            self.updatePsiMatrix = not self.updatePsiMatrix #Change which we update next time
            i +=1

    """ Updates Psi matrix by calling functions to update Psi internal points
    and Psi boundary points"""
    def updatePsi(self):
        for l in range(0, self.L+1):
            for j in range(0, self.L+1):
                if not self.inPlate(l,j):
                    #Side opposite plate
                    if(j == self.L):
                        self.Psi[l][j] = self.V0
                    #Upstream
                    elif(l == 0):
                        self.Psi[l][j] = self.V0*self.hY*j
                    #Downstream
                    elif(l == self.L):
                        self.Psi[l][j] = self.Psi[l-1][j]
                    #Plate side
                    elif(j == 0):
                        self.Psi[l][j] = 0
                    #Front of plate
                    elif((l*self.hX == self.PlateXFront) and (j*self.hY < self.PlateYTop)):
                        self.Psi[l][j] = 0
                    #Back of plate
                    elif((l*self.hX == self.PlateXBack) and (j*self.hY < self.PlateYTop)):
                        self.Psi[l][j] = 0
                    #Top of plate
                    elif((j*self.hY) == self.PlateYTop) and (l*self.hX >= self.PlateXFront) and (l*self.hX <= self.PlateXBack):
                        self.Psi[l][j] = 0
                    #Interal Point Not on Boundary
                    else:
                        W = self.W[l,j]
                        PsiSquareStencil = self.PsiSquareStencil(l,j)
                        Psi = ((1-self.w) * self.Psi[l,j]) + (self.w * PsiSquareStencil) + (W *self.hX*self.hY)
                        self.Psi[l,j] = Psi
                else:
                    #set values inside plate to 0
                    self.Psi[l,j] = 0

    """ Updates W matrix by calling functions to update Psi internal points
    and Psi boundary points"""
    def updateW(self):
        for l in range(0, self.L+1): #y-ccordinate j/L
            for j in range(0, self.L+1): #x-coordinate l/L
                if not self.inPlate(l,j):
                    #Side opposite plate
                    if(j == self.L):
                        self.W[l][j] = 0
                     #Upstream
                    elif(l == 0):
                        self.W[l][j] = 0
                     #Downstream
                    elif(l == self.L):
                        """TODO"""
                        self.W[l][j] = self.W[l-1][j]
                    #Plate side
                    elif(j == 0):
                        self.W[l][j] = 0
                    #Front of plate
                    elif((l*self.hX == self.PlateXFront) and (j*self.hY) < self.PlateYTop):
                        self.W[l][j] = (-2/(self.hY*self.hY))*self.Psi[l-1][j]
                        print("---")
                        print(l, ", ", j)
                        print("front of plate: ", self.W[l][j])
                    #Back of plate
                    elif((l*self.hX == self.PlateXBack) and (j*self.hY) < self.PlateYTop):
                        self.W[l][j] = (-2/(self.hY*self.hY))*self.Psi[l+1][j]
                        print("---")
                        print(l, ", ", j)
                        print("back of plate: ", self.W[l][j])
                    #Top of plate
                    elif((j*self.hY) == self.PlateYTop) and (l*self.hX >= self.PlateXFront) and (l*self.hX) <= self.PlateXBack:
                        self.W[l][j] = (-2/(self.hX*self.hX))*self.Psi[l][j+1]
                        print("---")
                        print(l, ", ", j)
                        print("top of plate: ", self.W[l][j])
                    else:
                        partial = self.w/(4*self.viscosity) * (self.PsiStencilDy(l,j) * self.WStencilDx(l,j) -
                        (self.PsiStencilDx(l,j) * self.WStencilDy(l,j)))
                        squareStencil = self.WSquareStencil(l,j)
                        W = ((1-self.w)*self.W[l,j]) + (self.w)*squareStencil - partial
                        self.W[l,j] = W
                        print("internal point: ", self.W[l][j])
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
        conlines = np.linspace(0, 1, 30)
        X,Y = self.getXYCoords(self.Psi)
        plt.contour(X,Y,self.Psi, conlines, cmap = 'plasma')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.title(r"$\Psi$")
        plt.show()


free = Flow(1, 1, 50)
free.initialize()
#free.updateWBoundary()
#free.updatePsiBoundary()

free.solve()
free.flowGraph()
