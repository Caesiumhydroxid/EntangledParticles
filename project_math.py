from scipy import constants
import numpy as np
from scipy import sparse
from scipy import integrate

def generateGaussianWavePacket(x, x0, a, k0):
    return (1/(2* np.pi * a**2))**(1/4) * np.exp(-(x-x0)**2/(4*a**2)) * np.exp(1j*k0*x)

def generateMagneticVectorPotential(zs,B0,zOn, fadelength=1e-8):
    index = 0
    Bs = np.zeros(len(zs))
    for i,z in enumerate(zs):
        if(z<zOn):
            Bs[i] = 0
        elif(z>= zOn and z<zOn+fadelength):
            Bs[i] = B0
        else:
            Bs[i] = 0
    As = integrate.cumtrapz(Bs,zs,initial=0)
    return As,Bs

    

def mapZSToLin(z,s):
    return int(z*(2)+s)

def zsFromLinear(n):
    return (n//2,int(n%(2)))

N = 201

def mapZ12StoLin(z1,z2,s,Nz=N):
    return int(z2*(Nz*2)+2*z1+s)

def z12SFromLinear(n,Nz=N):
    return ((n//2)%Nz,n//2//Nz,int(n%2))

print(z12SFromLinear(6,3))

# Coordinate System 
#   y  
#   | /x
#   |/
#   ----->z

# The homogenous magnetic field is oriented in y-direction
# Hence, we have a vector Potential in x Direction A(z) = H(z-z0)*B₀(z-z0) ex => B = H(z-z0)*B₀ * ey)  

sigma1 = [[0,1],[1,0]]
sigma2 = [[0,-1j],[1j,0]]
sigma3 = [[1,0],[0,-1]]

# Since the B-Field is only directed in the y direction, we only actually need sigma2

g = 2
def hamiltonian(zs, delta, As):
    N = len(zs)
    H = sparse.lil_matrix((N*2,N*2), dtype=complex)
    q = constants.elementary_charge #electrons

    for i,z in enumerate(zs):
        if(i==0):
            B = (As[i+1]-As[i])/(delta)
        elif(i == N-1):
            B = (As[i]-As[i-1])/(delta)
        else:
            B = As[i-1]-As[i+1]/(2*delta)
        
        m = constants.m_e
        # Spinor ⬆
        if(i>0):   
              H[mapZSToLin(i,0),mapZSToLin(i-1,0)] = -(constants.hbar)**2/(delta)**2 
        #if(i>0):
        #    H[mapZSToLin(i,0),mapZSToLin(i-1,0)] = - 2*(1j * constants.hbar)/ (2*delta) * q * As[i]
        H[mapZSToLin(i,0),mapZSToLin(i,0)] = 2*(constants.hbar)**2/(delta)**2 + q**2 * As[i]**2
        #if(i+1<N):
        #    H[mapZSToLin(i,0),mapZSToLin(i+1,0)] =  2*(1j * constants.hbar)/ (2*delta) * q * As[i]
        if(i+1<N):
            H[mapZSToLin(i,0),mapZSToLin(i+1,0)] = -(constants.hbar)**2/(delta)**2

        # Spin, magnetic field coupling
        H[mapZSToLin(i,0),mapZSToLin(i,0)] -= g*(q * constants.hbar) * B * sigma1[0][0]/2
        H[mapZSToLin(i,0),mapZSToLin(i,1)] -= g*(q * constants.hbar) * B * sigma1[0][1]/2

        # Spinor ⬇

        if(i>0):   
            H[mapZSToLin(i,1),mapZSToLin(i-1,1)] = -(constants.hbar)**2/(delta)**2 
        #if(i>0):
        #    H[mapZSToLin(i,0),mapZSToLin(i-1,0)] = - 2*(1j * constants.hbar)/ (2*delta) * q * As[i]
        H[mapZSToLin(i,1),mapZSToLin(i,1)] = 2*(constants.hbar)**2/(delta)**2 + q**2 * As[i]**2
        #if(i+1<N):
        #    H[mapZSToLin(i,0),mapZSToLin(i+1,0)] =  2*(1j * constants.hbar)/ (2*delta) * q * As[i]
        if(i+1<N):
            H[mapZSToLin(i,1),mapZSToLin(i+1,1)] = -(constants.hbar)**2/(delta)**2

        # Spin, magnetic field coupling
        H[mapZSToLin(i,1),mapZSToLin(i,0)] -= g*(q * constants.hbar) * B * sigma1[1][0]/2
        H[mapZSToLin(i,1),mapZSToLin(i,1)] -= g*(q * constants.hbar) * B * sigma1[1][1]/2
    
    return (1/(2*constants.m_e))*H.tocsc()


g = 2
def hamiltonian2Particles(zs, delta, As):
    N = len(zs)
    H = sparse.lil_matrix((N*N*2,N*N*2), dtype=complex)
    q = constants.elementary_charge
    for z1 in range(len(zs)):
        for z2 in range(len(zs)):

            m = constants.m_e
            if(z1-1 < 0):
                Bz1 = As[z1+1]-As[z1]/ delta
            elif(z1+1 >= N): 
                Bz1 = As[z1]-As[z1-1]/ delta
            else:
                Bz1 = As[z1+1]-As[z1-1]/ (2*delta)

            if(z2-1 < 0):
                Bz2 = As[z2+1]-As[z2]/ delta
            elif(z2+1 >= N): 
                Bz2 = As[z2]-As[z2-1]/ delta
            else:
                Bz2 = As[z2+1]-As[z2-1]/ (2*delta)

            # Spinor ⬆
            if(z1>0):
                H[mapZ12StoLin(z1,z2,0),mapZ12StoLin(z1-1,z2,0)] = -(constants.hbar)**2/(delta**2) + (1j * constants.hbar)/ (delta) * q * As[z1]
            if(z2>0):
                H[mapZ12StoLin(z1,z2,0),mapZ12StoLin(z1,z2-1,0)] =-(constants.hbar)**2/(delta**2) + (1j * constants.hbar)/ (delta) * q * As[z2]
            H[mapZ12StoLin(z1,z2,0),mapZ12StoLin(z1,z2,0)] = 4*(constants.hbar)**2/(delta**2) + q**2 * As[z1]**2 + q**2 * As[z2]**2
            if(z1+1<N):
                H[mapZ12StoLin(z1,z2,0),mapZ12StoLin(z1+1,z2,0)] = -(constants.hbar)**2/(delta**2) - (1j * constants.hbar)/ delta * q * As[z1]
            if(z2+1<N):
                H[mapZ12StoLin(z1,z2,0),mapZ12StoLin(z1,z2+1,0)] = -(constants.hbar)**2/(delta**2) - (1j * constants.hbar)/ delta * q * As[z2]

            H[mapZ12StoLin(z1,z2,0),mapZ12StoLin(z1,z2,0)] -= g*(q * constants.hbar) * (Bz1 ) * sigma2[0][0]/2
            H[mapZ12StoLin(z1,z2,0),mapZ12StoLin(z1,z2,1)] -= g*(q * constants.hbar) * (Bz1 ) * sigma2[0][1]/2

            # Spinor ⬇
            if(z1>0):
                H[mapZ12StoLin(z1,z2,1),mapZ12StoLin(z1-1,z2,1)] = -(constants.hbar)**2/(delta**2) + (1j * constants.hbar)/ (delta) * q * As[z1]
            if(z2>0):
                H[mapZ12StoLin(z1,z2,1),mapZ12StoLin(z1,z2-1,1)] = -(constants.hbar)**2/(delta**2) + (1j * constants.hbar)/ (delta) * q * As[z2]
            H[mapZ12StoLin(z1,z2,1),mapZ12StoLin(z1,z2,1)] = 4*(constants.hbar)**2/(delta**2) + q**2 * As[z1]**2 + q**2 * As[z2]**2
            if(z1+1<N):
                H[mapZ12StoLin(z1,z2,1),mapZ12StoLin(z1+1,z2,1)] = -(constants.hbar)**2/(delta**2) - (1j * constants.hbar)/ delta * q * As[z1]
            if(z2+1<N):
                H[mapZ12StoLin(z1,z2,1),mapZ12StoLin(z1,z2+1,1)] = -(constants.hbar)**2/(delta**2) - (1j * constants.hbar)/ delta * q * As[z2]

            H[mapZ12StoLin(z1,z2,1),mapZ12StoLin(z1,z2,0)] -= g*(q * constants.hbar) * (Bz1 ) * sigma2[1][0]/2
            H[mapZ12StoLin(z1,z2,1),mapZ12StoLin(z1,z2,1)] -= g*(q * constants.hbar) * (Bz1 ) * sigma2[1][1]/2

    return (1/(2 * constants.m_e))*H.tocsc()