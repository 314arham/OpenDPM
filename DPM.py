import math as m
import numpy as np
import scipy as sp

def Minimum(array):
    """
    Return minimum value and argmin for a given array.
    """                                                                         #returns [argmin, min]
    return [[i for i in range(len(array)) if array[i]==min(array)][0], min(array)]

def Rotate(vector, theta): 
    """
    Counter-Clockwised rotatation of xy-Cartesian coordinate system through the angle theta
    """                                                                 #rotation of coordinate system
    [x, y] = vector
    return np.array([x*m.cos(theta) - y*m.sin(theta), x*m.sin(theta) + y*m.cos(theta)])

def Abs(a, b):
    """
    Return absolute value of a vector (a, b)
    """                                                                              #distance between 2 points
    return m.sqrt((a.x-b.x)**2 + (a.y-b.y)**2)

class geometry:
    """
    Geometry object
    ver 1.1
    Events are occured inside geometry. length, width and center should be specified.
    """
    def __repr__(self):
        S = str(
                        'class          '   +   'geometry'                  +   '\n'    +
                        'name           '   +       self.name               +   '\n'    +
                        'length         '   +   str(self.length)            +   '\n'    +
                        'width          '   +   str(self.width)             +   '\n'    +
                        'center         '   +   str(self.center)            +   '\n'
                )
        return S
    def __init__(self, Length, Width, Center = [0, 0], Name = 'geometry01'):
        self.length, self.width, self.center, self.name = Length, Width, Center, Name
    def __eq__(self, other):
        return self.length == other.length and self.width == other.width and self.center == other.center
    def area(self):
        return self.length*self.width

class particle:
    def __repr__(self):
        S = str(
                        'class          '   +   'particle'                  +   '\n'    +
                        'type           '   +   'hard sphere'               +   '\n'    +
                        'radius         '   +   str(self.r)                 +   '\n'    +
                        'COG            '   +   str([self.x, self.y])       +   '\n'    +
                        'velocity       '   +   str([self.u, self.v])       +   '\n'    +
                        'volume         '   +   str(self.vol)               +   '\n'    +
                        'mass           '   +   str(self.M)                 +   '\n'    +
                        'omega          '   +   str(self.Omega)             +   '\n'    +
                        'MOI            '   +   str(self.momentOfInertia)   +   '\n'    +
                        'density        '   +   str(self.rho)               +   '\n'
                )
        return S
    def __init__(self, X, Y, Radius=0.1, velocityVector=[0, 0], Omega=0, density=2600, GeometryObj = geometry(10, 10)):
        self.r, self.x, self.y, self.eps_p, self.eps_w, self.Omega = Radius, X + GeometryObj.center[0], Y + GeometryObj.center[1], 0.8, 0.8, Omega
        [self.u, self.v] = velocityVector
        self.rho = density
        self.GeometryObj = GeometryObj
        self.__props()
    def __add__(self, other):
        if self.GeometryObj == other.GeometryObj:
            p = particle(
                            self.x, 
                            self.y, 
                            self.r + other.r, 
                            [self.u + other.u, self.v + other.v], 
                            self.Omega + other.Omega, 
                            self.rho + other.rho, 
                            self.GeometryObj
                        )
        return p
    def __neg__(self):
        self.u *= -1
        self.v *= -1
        return self
    def __mul__(self, num):
        self.r *= num
        return self
    def __props(self):
        self.vol = sp.pi*(4/3)*self.r**3
        self.M = self.rho*self.vol
        self.momentOfInertia = 0.4*self.M*self.r**2
        return True
    def move(self, dt):
        self.x += self.u*dt
        self.y += self.v*dt
    def gravity(self, dt):
        self.v += -9.81*dt
    def distanceFromWall(self):
        d = []
        U = [-self.u, self.u, -self.v, self.v]
        D = [self.x, -self.x+self.GeometryObj.length, self.y, -self.y+self.GeometryObj.width]
        D = np.array(D)-np.array([self.r]*4)
        for i in range(len(U)):
            if U[i] != 0 and D[i]/U[i]>0:
                d += [D[i]/U[i]]
            else:
                d += [np.inf]
        return Minimum(d)
    def collisionDynamics(self, p):
        if type(p) == int:
            if p==0 or p==1:
                theta = np.deg2rad(90)
            elif p==2 or p==3:
                theta = 0
            rotated_velocity_self = Rotate([self.u, self.v], theta)
            S0 = (rotated_velocity_self[0] + self.Omega*self.r)
            B1 = 1/self.M + self.r**2/self.momentOfInertia
            rotated_velocity_self[0] += -S0/(B1*self.M)
            C0 = rotated_velocity_self[1]
            B2 = 1/self.M
            rotated_velocity_self[1] += -(1 + self.eps_p)*C0/(B2*self.M)
            self.Omega += -S0*self.r/(B1*self.momentOfInertia)
            [self.u, self.v] = Rotate(rotated_velocity_self, -theta)
            return 'wall'
        else:
            theta = np.arctan(np.true_divide(self.x - p.x, self.y - p.y))
            rotated_velocity_self = Rotate([self.u, self.v], theta)
            rotated_velocity_p = Rotate([p.u, p.v], theta)
            S0 = (rotated_velocity_self[0] + self.Omega*self.r) - (rotated_velocity_p[0] - p.Omega*p.r)
            B1 = 1/self.M + self.r**2/self.momentOfInertia + 1/p.M + p.r**2/p.momentOfInertia
            rotated_velocity_self[0] += -S0/(B1*self.M)
            rotated_velocity_p[0] += S0/(B1*p.M)
            C0 = rotated_velocity_self[1] - rotated_velocity_p[1]
            B2 = 1/self.M + 1/p.M
            rotated_velocity_self[1] += -(1 + self.eps_p)*C0/(B2*self.M)
            rotated_velocity_p[1] += (1 + self.eps_p)*C0/(B2*p.M)
            self.Omega += -S0*self.r/(B1*self.momentOfInertia)
            p.Omega += -S0*p.r/(B1*p.momentOfInertia)   
            [self.u, self.v] = Rotate(rotated_velocity_self, -theta)
            [p.u, p.v] = Rotate(rotated_velocity_p, -theta)
        return 'particle'

class solidPhase:
    def __init__(self, Gravity = True, GeometryObj = geometry(10, 10), name='Project01'):
        self.save = False
        self.g_const = -9.81
        self.GeometryObj = GeometryObj
        self.projectName = name
        self.isInGravityField = Gravity
    def arrayLattice(self, Mapping=[3, 3, 1], Rlims=[0.6, 0.2] , Rho = 2600, velocitySeeds=[2, -2], OmegaSeed=0):
        """
        Constructs a n*m grid of particles, where n = Mapping[0] and m = Mapping[1],
        Radii are randomly generated between Rlims[0] and Rlims[0] + Rlims[1]
        Velocities are randomly generated between velocitySeed[0] and 
        velocitySeeds[1]. Omega = sp.rand()*OmegaSeed
        """
        R = sum(Rlims)
        self.p = []
        self.p = np.array(
                            [
                                particle(
                                            X = i, 
                                            Y = j, 
                                            Radius = sp.rand()*Rlims[0]+Rlims[1], 
                                            velocityVector = np.array([sp.rand()*velocitySeeds[0]-velocitySeeds[0]/2, 
                                                                    sp.rand()*velocitySeeds[1]-velocitySeeds[1]/2]), 
                                            Omega = 0, 
                                            density = Rho, 
                                            GeometryObj = self.GeometryObj
                                        ) 
                                for i in np.linspace(R, self.GeometryObj.length - R, Mapping[0]) 
                                for j in np.linspace(R, self.GeometryObj.width - R, Mapping[1])
                            ]
                        )
        self.numberOfParticles = len(self.p)
    def custom(self, particlesList):
        """
        User-defined particles are constructed to a solid phase.
        
        >>> a = particle(4, 3, 0.1)
        >>> q = particle(3.85, 5, 0.2)
        >>> phase = solidPhase()
        >>> phase.custom([particle1, particle2])
        """
        self.p = []
        self.p = np.array(particlesList)
        self.numberOfParticles = len(self.p)
        return True
    def saveAs(self, projectName = 'Project01', saveTo = 'Documents/CSV/recent/', dataStepSize = 1):
        """
        Data will not be saved unless .saveAs() method be called.
        Note: step size of output *.csv files could be changed by dataStepSize value. 
        """
        self.save = True
        self.projectName, self.saveTo, self.dataStepSize = projectName, saveTo, dataStepSize
    def collisionListGenerator(self):
        """
        Collision list is consist of nearest particle-particle or particle-wall collision times.
        Each element of it, is a list that gives information of a certain event.
        [particle A, particle B, time] or [particle A, int wallNumber, time]
        """
        collisionTimesList = []
        for i in range(self.numberOfParticles):
            for j in range(i+1, self.numberOfParticles):
                if Abs(self.p[i], self.p[j]) <= 5*max(self.p[i].r, self.p[j].r):
                    R12 = [self.p[i].x - self.p[j].x, self.p[i].y - self.p[j].y]
                    V12 = [self.p[i].u - self.p[j].u, self.p[i].v - self.p[j].v]
                    Delta = np.dot(R12, V12)**2 - (V12[0]**2 + V12[1]**2)*((R12[0]**2 + R12[1]**2) - (self.p[i].r + self.p[j].r)**2)
                    A = -np.dot(R12, V12)
                    B = np.dot(V12, V12)
                    if A>0 and B!=0 and Delta>=0 and (A-m.sqrt(Delta))/B>0:
                        collisionTimesList += [[self.p[i], self.p[j], np.true_divide(A-m.sqrt(abs(Delta)), B)]]
                    else:
                        collisionTimesList += [[i, j, np.inf]]
                else:
                    collisionTimesList += [[i, j, np.inf]]
        wallCollisionList = [[i, i.distanceFromWall()[0], i.distanceFromWall()[1]] for i in self.p]
        return collisionTimesList + wallCollisionList
    def run(self, TIME=5, DeltaT=0.01):
        Time = 0.0
        DT = DeltaT
        timeStep = 0
        ROUNDOFF_ERROR_CORECTION = 0.99
        while Time <= TIME:
            acctim = 0
            collisionList = self.collisionListGenerator()                                           ### set up collision list
            collisionTimes = [t[2] for t in collisionList]                                          #Extrcting collision times from collisionList
            minimumCollisionData = Minimum(collisionTimes)                                          #Find minimum collision time in collision list. minimumCollisionData = [number of element in array, element itself]
            [primary, secondary, tab] = collisionList[minimumCollisionData[0]]                      ### locate minimum collision time "tab"
            acctim += ROUNDOFF_ERROR_CORECTION*tab                                                  ### increment acctim by tab
            while acctim < DT:
                for particles in self.p: particles.move(ROUNDOFF_ERROR_CORECTION*tab)
                primary.collisionDynamics(secondary)                                                ### collision dynamics
                collisionList = self.collisionListGenerator()                                       ### reset list
                collisionTimes = [t[2] for t in collisionList]
                minimumCollisionData = Minimum(collisionTimes)
                [primary, secondary, tab] = collisionList[minimumCollisionData[0]]                  ### locate minimum collision time
                acctim += ROUNDOFF_ERROR_CORECTION*tab                                              ### increment acctim by tab
            for particles in self.p: particles.move(DT - (acctim - ROUNDOFF_ERROR_CORECTION*tab))
            if self.isInGravityField == True: 
                for particles in self.p: particles.gravity(DT)
            Time += DT
            if self.save == True and timeStep % self.dataStepSize == 0:                             ### Write the output data in *.csv files, by a given timeStep
                file = open(
                                self.saveTo + self.projectName      +       '.csv.'     +
                                str(timeStep//self.dataStepSize),           mode='w'
                            )
                file.write('h,x,z,d\n')                                                             ### Order of variables in file Height, x, z=0, diameter, time. z is given so visualization in Paraview will have a better quality
                for i in range(self.numberOfParticles):
                    file.write(
                                    str(round(self.p[i].y, 3))      +       ','         +
                                    str(round(self.p[i].x,3))       +       ','         +
                                                '0'                 +       ','         +
                                    str(round(2*self.p[i].r,3))     +       '\n'
                                )
            timeStep+=1
