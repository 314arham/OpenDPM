import math as m
import numpy as np
import scipy as sp

#DPM Hard sphere approach event driven                             #Based on B.P.B Hoomans et al. 1995

def Minimum(a):
    return [[i for i in range(len(a)) if a[i]==min(a)][0], min(a)]

def Rotate(vector, theta):
    [x, y] = vector
    return np.array([x*m.cos(theta) - y*m.sin(theta), x*m.sin(theta) + y*m.cos(theta)])

class particle:
    def __init__(self, X, Y, Radius=0.1, velocityVector=[0, 0], Omega=0, Mass=0.1, Cube=[10, 10, 10]):
        self.r, self.x, self.y = Radius, X, Y
        [self.u, self.v] = velocityVector
        self.velocityVector = velocityVector
        self.M = Mass
        self.I = 0.4*self.M*self.r**2
        self.eps_w = 0.8
        self.eps_p = 0.8
        self.Omega = Omega
        self.cube = Cube
    def volume(self):
        return sp.pi*(4/3)*self.r**3
    def move(self, dt):
        self.x += self.u*dt
        self.y += self.v*dt
    def dist(self):
        d = []
        U = [-self.u, self.u, -self.v, self.v]
        D = [self.x, -self.x+self.cube[0], self.y, -self.y+self.cube[1]]
        D = np.array(D)-np.array([self.r]*4)                                                            #Radius term
        for i in range(len(U)):
            if U[i] != 0 and D[i]/U[i]>0:
                d += [D[i]/U[i]]
            else:
                d += [np.inf]
        return Minimum(d)
    def collisionDynamics(self, p):           #19:01 Tue
        #print(type(p))
        if type(p) == int:
            if p==0 or p==1:
                self.u = -self.u*self.eps_w
            elif p==2 or p==3:
                self.v = -self.v*self.eps_w
            return 'wall'
        else:
            eps_p = 0.8                                                                                     #eps_p
            theta = np.arctan(np.true_divide(self.x - p.x, self.y - p.y))
            rotated_velocity_self = Rotate([self.u, self.v], theta)
            rotated_velocity_p = Rotate([p.u, p.v], theta)
            S0 = (rotated_velocity_self[0] + self.Omega*self.r) - (rotated_velocity_p[0] - p.Omega*p.r)
            B1 = 1/self.M + self.r**2/self.I + 1/p.M + p.r**2/p.I
            rotated_velocity_self[0] += -S0/(B1*self.M)
            rotated_velocity_p[0] += S0/(B1*p.M)                                                            #Vx = u changed here; due to Hoomans et al. 2D
            C0 = rotated_velocity_self[1] - rotated_velocity_p[1]
            B2 = 1/self.M + 1/p.M
            rotated_velocity_self[1] += -(1 + eps_p)*C0/(B2*self.M)                                         #Vy = v changed here; due to Hoomans et al. 2D for a better resoloution in Paraview, it's better to simulate particles in XZ plane. for this, the gravity() and collision_p must be changed (e.g: self.v -> self.w)
            rotated_velocity_p[1] += (1 + eps_p)*C0/(B2*p.M)
            self.Omega += -S0*self.r/(B1*self.I)
            p.Omega += -S0*p.r/(B1*p.I)   
            [self.u, self.v] = Rotate(rotated_velocity_self, -theta)
            [p.u, p.v] = Rotate(rotated_velocity_p, -theta)
        return 'particle'

class solidPhase:
    def __init__(self):
        self.g_const = -9.81
    def sameSizeStationary(self, Radius, Cube=[10, 10, 10]):                   #initializes a stable square grid of same size particles
        self.p = []
        [Length, Width, Height] = Cube
        for i in np.arange(Radius, Length, 2*Radius):
            for j in np.arange(Radius, Width, 2*Radius):
                self.p.append(particle(i, j, Radius))
        self.p = np.array(self.p)
        self.numberOfParticles = len(self.p)
    def arrayLattice(self, Mapping=[3, 3, 1], velocitySeeds=[2, -2], OmegaSeed=0, Rlims=[0.6, 0.2], Cube=[10, 10, 10]):
        R = sum(Rlims)
        self.p = []
        self.p = np.array([particle(X = i, 
                                    Y = j, 
                                    Radius = sp.rand()*Rlims[0]+Rlims[1], 
                                    velocityVector = np.array([sp.rand()*velocitySeeds[0]-velocitySeeds[0]/2, 
                                                               sp.rand()*velocitySeeds[1]-velocitySeeds[1]/2]), 
                                    Omega = 0, 
                                    Mass = 1, 
                                    Cube = Cube) 
                                    for i in np.linspace(R, Cube[0]-R, Mapping[0]) 
                                    for j in np.linspace(R, Cube[1]-R, Mapping[1])])
        self.numberOfParticles = len(self.p)
    def custom(self, particlesList):                                           #a solid phase with custom given particles will be constructed
        self.p = []
        self.p = np.array(particlesList)
        self.numberOfParticles = len(self.p)
    def gravity(self, dt):
        for particle in self.p:
            particle.v += self.g_const*dt
    def collisionListGenerator(self):
        collisionTimesList = []
        for i in range(self.numberOfParticles):
            for j in range(i+1, self.numberOfParticles):
                R12 = [self.p[i].x - self.p[j].x, self.p[i].y - self.p[j].y]
                V12 = [self.p[i].u - self.p[j].u, self.p[i].v - self.p[j].v]
                Delta = np.dot(R12, V12)**2 - (V12[0]**2 + V12[1]**2)*((R12[0]**2 + R12[1]**2) - (self.p[i].r + self.p[j].r)**2)
                A = -np.dot(R12, V12)
                B = np.dot(V12, V12)
                #if Delta<0: print(Delta)
                if A>0 and B!=0 and Delta>=0 and (A-m.sqrt(Delta))/B>0:
                    collisionTimesList += [[self.p[i], self.p[j], np.true_divide(A-m.sqrt(abs(Delta)), B)]]
                else:
                    collisionTimesList += [[i, j, np.inf]]
        wallCollisionList = [[i, i.dist()[0], i.dist()[1]] for i in self.p]
        return collisionTimesList + wallCollisionList
    def run(self, TIME=5):
        n = self.numberOfParticles
        Time = 0.0
        DT = 0.01
        timeStep = 0
        dataStepSize = 1
        ROUNDOFF_ERROR_CORECTION = 0.99
        while Time<=TIME:
            acctim = 0
            collisionList = self.collisionListGenerator()                                    ### set up collision list
            collisionTimes = [t[2] for t in collisionList]                                   #Extrcting collision times from collisionList
            minimumCollisionData = Minimum(collisionTimes)                                   #Find minimum collision time in collision list. minimumCollisionData = [number of element in array, element itself]
            [primary, secondary, tab] = collisionList[minimumCollisionData[0]]               ### locate minimum collision time "tab"
            acctim += ROUNDOFF_ERROR_CORECTION*tab                                           ### increment acctim by tab
            while acctim < DT:
                for particles in self.p: particles.move(ROUNDOFF_ERROR_CORECTION*tab)
                x = primary.collisionDynamics(secondary)                                     ### collision dynamics
                collisionList = self.collisionListGenerator()                                ### reset list
                collisionTimes = [t[2] for t in collisionList]
                minimumCollisionData = Minimum(collisionTimes)
                [primary, secondary, tab] = collisionList[minimumCollisionData[0]]           ### locate minimum collision time
                acctim += ROUNDOFF_ERROR_CORECTION*tab                                       ### increment acctim by tab
            for particles in self.p: particles.move(DT - (acctim - ROUNDOFF_ERROR_CORECTION*tab))
            self.gravity(DT)
            Time += DT
            if timeStep % dataStepSize == 0:                                                 ### Write the output data in *.csv files, by a given timeStep
                file = open('F:/CSV/N/P.csv.'+
                            str(timeStep//dataStepSize), mode='w')                           # Change the address #
                file.write('h,x,z,d,Time\n')                                                 ### Order of variables in file Height, x, z=0, diameter, time. z is given so visualization in Paraview will have a better quality
                for i in range(n):
                    file.write(str(round(self.p[i].y, 3))+','+
                               str(round(self.p[i].x,3))+',0,'+
                               str(round(2*self.p[i].r,3))+','+
                               str(Time)+'\n')
            timeStep+=1

F = solidPhase()
p = particle(1, 9, 0.1, [0, -10], 0)
q = particle(1.15, 8, 0.1, [0, 0], 0)
F.custom([p, q])
#F.arrayLattice()
F.run(1)
