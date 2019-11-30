import numpy as np
import matplotlib.pyplot as plt
import math

class Flow:
    def __init__(self, X, Y, L):
        self.L = L #highest index in matrix
        self.X = X #length of surface in x direction
        self.Y = Y #length of the surface in y direction
        self.V0 = 1 #Velocity of incoming fluid
        self.w = 1.5 #over relaxation factor
        self.viscosity = 0.9 #viscosity of the fluid
        #Plate boundaries
        self.PlateXFront = (.3) * (L)
        self.PlateXBack = (.5) * (L)
        self.PlateYTop = (.3) * (L)
        #Vorticity matrix
        self.W = np.zeros((L+1,L+1))
        #Stream function matrix
        self.Psi = np.zeros((L+1,L+1))
        #step sizes
        self.hX = X/(L)
        self.hY = Y/(L)
        self.updatePsiMatrix = True #Whether w has been updates since last psi update
        self.solutionFound = False #True when within tolerance
        self.maxIterations = 200 #Max iterations before giving up
        self.residualNorms = []

    """Returns true if (l,j) is inside the plate"""
    def inPlate(self, l, j):
        front = l > self.PlateXFront
        back = l < self.PlateXBack
        top = j < self.PlateYTop
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
        i = 0
        while not self.solutionFound and i < self.maxIterations:
            if self.updatePsiMatrix:
                self.updatePsi()
            else:
                self.updateW()
            self.updatePsiMatrix = not self.updatePsiMatrix #Change which we update next time
            self.residualNorms += [self.getResidualNorm()]
            i +=1

    """ Updates Psi matrix by calling functions to update Psi internal points
    and Psi boundary points"""
    def updatePsi(self):
        for l in range(0, self.L+1):
            for j in range(0, self.L+1):
                if not self.inPlate(l,j):
                    #Side opposite plate
                    if(j == self.L):
                        self.Psi[l][j] = self.V0*self.Y
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
                    elif((l == self.PlateXFront) and (j < self.PlateYTop)):
                        self.Psi[l][j] = 0
                    #Back of plate
                    elif((l == self.PlateXBack) and (j < self.PlateYTop)):
                        self.Psi[l][j] = 0
                    #Top of plate
                    elif((j) == self.PlateYTop) and (l>= self.PlateXFront) and (l <= self.PlateXBack):
                        self.Psi[l][j] = 0
                    #Interal Point Not on Boundary
                    else:
                        W = self.W[l,j]
                        PsiSquareStencil = self.PsiSquareStencil(l,j)
                        Psi = ((1-self.w) * self.Psi[l,j]) + (self.w * PsiSquareStencil) + (((self.w*self.hX*self.hY)/4)*W)
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
                        self.W[l][j] = self.W[l-1][j]
                    #Plate side
                    elif(j == 0):
                        self.W[l][j] = 0
                    #Front of plate
                    elif(l == self.PlateXFront) and (j < self.PlateYTop):
                        self.W[l][j] = (-2/(self.hX*self.hX))*self.Psi[l-1][j]
                    #Back of plate
                    elif((l == self.PlateXBack) and (j) < self.PlateYTop):
                        self.W[l][j] = (-2/(self.hX*self.hX))*self.Psi[l+1][j]
                    #Top of plate
                    elif((j) == self.PlateYTop) and (l >= self.PlateXFront) and (l) <= self.PlateXBack:
                        self.W[l][j] = (-2/(self.hY*self.hY))*self.Psi[l][j+1]
                    else:
                        psiywx = self.PsiStencilDy(l,j) * self.WStencilDx(l,j)
                        psixwy = self.PsiStencilDx(l,j) * self.WStencilDy(l,j)
                        partial = ((self.w*self.hX*self.hY)/(self.viscosity*4))*(psiywx-psixwy)
                        avg = (self.w)*self.WSquareStencil(l,j)
                        factor = (1-self.w)*self.W[l,j]
                        final = factor+avg-partial
                        self.W[l,j] = final
                else:
                    self.W[l,j] = 0

    """Stencil for Laplacian of psi"""
    def PsiStencilLaPlaz(self,l,j):
        Psi = self.Psi[l,j]
        PsiSquareStencil = 0
        constant = 4
        if l != 0 and j != 0 and l !=  (self.L) and j != (self.L):
            if l >= self.PlateXFront and l <= self.PlateXBack and  j <= self.PlateYTop:
                return -self.W[l,j]
            else:
                PsiSquareStencil = (1/4)*(self.Psi[l+1,j] + self.Psi[l-1,j] + self.Psi[l,j+1]+ self.Psi[l,j-1])
                stencil = constant*(PsiSquareStencil - Psi)/(self.hX*self.hY)
                return stencil
        else:
            return -self.W[l,j]

    def PsiSquareStencil(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L) and j != (self.L):
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
            stencil = (self.Psi[l+1,j] - self.Psi[l-1,j]) /(2*self.hX)
            return stencil

    def WStencilDy(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (self.W[l,j+1] - self.W[l,j-1])/(2*self.hY)
            return stencil

    def WStencilDx(self,l,j):
        if l != 0 and j != 0 and l !=  (self.L +1) and j != (self.L +1):
            stencil = (self.W[l+1,j] - self.W[l-1,j])/(2*self.hX)
            return stencil

    def getResidualMatrix(self):
        residualMatrix = np.zeros((self.L+1,self.L+1))
        for j in range(0, self.L+1):
            for l in range(0, self.L+1):
                residual = self.PsiStencilLaPlaz(l,j) + self.W[l,j]
                residualMatrix[l,j] = residual
        return residualMatrix

    def getResidualNorm(self):
        residualMatrix = self.getResidualMatrix()
        norm = 0
        for j in range(0, self.L+1):
            for l in range(0, self.L+1):
                nij = (residualMatrix[l,j]**2)*(self.hX*self.hY)
                norm += nij
        norm = math.sqrt(norm)
        return norm

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

    def graph(self):
        self.flowGraph()
        self.vorticityGraph()
        self.residualGraph()
        self.normGraph()

    """Graph contour lines of psi"""
    def flowGraph(self):
        #graph stream function
        plt.cla()
        conlines = np.linspace(0, 1, 50)
        X,Y = self.getXYCoords(self.Psi)
        plt.contour(X,Y,self.Psi, conlines, cmap = 'plasma')
        title = r"$\Psi$ (" + str(self.maxIterations) +  " Iterations, N = " + str(self.L) + r" $\mu$=" + str(self.viscosity) +  ")"
        plt.title(title)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

    def vorticityGraph(self):
        #graph vorticity
        plt.cla()
        plt.imshow(self.W.T, cmap='viridis',origin = 'lower')
        title = r"$\omega$ (" + str(self.maxIterations) +  " Iterations, N = " + str(self.L) + r" $\mu$=" + str(self.viscosity) +  ")"
        plt.title(title)
        x_positions = np.linspace(0,len(self.W),5)
        x_labels = np.linspace(0,self.X,5)
        x_labels = [round(x,2) for x in x_labels]
        plt.xticks(x_positions, x_labels)
        y_positions = x_positions
        y_labels = np.linspace(0,self.Y,5)
        y_labels = [round(y,2) for y in y_labels]
        plt.yticks(y_positions, y_labels)
        plt.colorbar()
        plt.show()

    def residualGraph(self):
        #graph risidual
        plt.cla()
        residualMatrix = self.getResidualMatrix()
        x_positions = np.linspace(0,len(self.W),5)
        x_labels = np.linspace(0,self.X,5)
        x_labels = [round(x,2) for x in x_labels]
        y_positions = x_positions
        y_labels = np.linspace(0,self.Y,5)
        y_labels = [round(y,2) for y in y_labels]
        plt.imshow(residualMatrix.T, cmap='viridis',origin = 'lower')
        title = r"$r_{\Psi}$ (" + str(self.maxIterations) +  " Iterations, N = " + str(self.L) + r" $\mu$=" + str(self.viscosity) +  ")"
        plt.title(title)
        plt.xticks(x_positions, x_labels)
        plt.yticks(y_positions, y_labels)
        plt.colorbar()
        plt.show()

    def normGraph(self):
        #graph risidual norm
        plt.cla()
        plt.plot(self.residualNorms,label="norm")
        title = r"$R_{\Psi}$ (" + str(self.maxIterations) +  " Iterations, N = " + str(self.L) + r" $\mu$=" + str(self.viscosity) +  " w=" + str(self.w) + " )"
        plt.title(title)
        plt.xlabel("Iteration")
        plt.ylabel(r"$R_{\Psi}$")
        plt.show()

free = Flow(1, 1, 30)
free.solve()
free.graph()
