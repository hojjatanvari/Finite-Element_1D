#HW 1
#MEMA 646
#Hojjat Anvari
#Prof. Strouboulis
#Texas A&M University
#-----------------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
###-------------- Question 1: -u'' = x     |     u(0) = 0 , u(1) = 0 --------------###

#----- Exact Solution -----#
x = np.linspace(0,1,100)
u_exact = 1/6 * (x - x**3)
L = 1
#----- Weak Formulation -----#
# Global K Matrix
def K_global(K_g, n):
    for row in range (1,n):
        if (row == 1):
            K_g[0][0] = Ke[1][1]
        else:
            K_g[row -2][row - 2] += Ke[0][0]
            K_g[row-1][row-1] = Ke[1][1]
            K_g[row -2][row-1] = Ke[0][1]
            K_g[row-1][row-2] = Ke[1][0]
        if(row == n-1):
            K_g[row-1][row-1] += Ke[0][0]

# Global F vector
def F_global(F_g, n):
    for row in range (1,n):
        if (row == 1):
            Xa = 0
            Xb = h
            fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
            F_g[0][0] = fe[1][0]
        else:
            Xa += h
            Xb += h
            fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
            F_g[row-2][0] += fe[0][0]
            F_g[row-1][0] += fe[1][0] 
        if(row == n-1):
            Xa += h
            Xb += h
            fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])            
            F_g[row-1][0] += fe[0][0]
#-----------------------------------------------------------------------------------------------------------------------------------------#
print("--------------------------------------- Question 1: -u'' = x     |     u(0) = 0 , u(1) = 0 ---------------------------------------")
#-----------------------------------------------||Uex||^2 -------------------------------------------------#
print("\n")
def d_uex_2(xx):
    return (1/6 - (xx - xx**3))**2
uex_E_2, _ = spi.quad(d_uex_2, 0, 1)
print("\n||Uex||^2:")
print(uex_E_2)
# ------------------------------------------------ n = 2 ------------------------------------------------ #
n = 2
h = L/n
# Local element stiffness and load vectors
Xa = 0
Xb = h
Ke = np.array([[1/h, -1/h], [-1/h, 1/h]])
fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
#--------------------------------------
K_g2 = np.zeros((n-1,n-1))
F_g2 = np.zeros((n-1,1))
K_global(K_g2, n)
F_global(F_g2, n)
U2 = np.linalg.inv(K_g2) @ F_g2
U2 = np.insert(U2, 0, 0, axis = 0) # Boundary condition
U2 = np.insert(U2, n, 0, axis = 0)
F_g2 = np.insert(F_g2, 0, 0, axis = 0)
F_g2 = np.insert(F_g2, n, 0, axis = 0)
Uh2 = 1/2 * np.transpose(U2) @ F_g2
error_2_2 = abs(uex_E_2 - Uh2)
print("\n")
print("-------------------------- Number of elements: 2 --------------------------")
print("Global stiffness matrix(K_g):")
print(K_g2)
print("\nGlobal load vector(F_g):")
print(F_g2)
print("\nThe finite element approximation (uh) coefficients:")
print(U2)
print("\nEnergy in the finite element solution (1/2*uT*F):")
print(Uh2)
print("\n||eE||^2:")
print(error_2_2)
# ------------------------------------------------ n = 4 ------------------------------------------------ #
n = 4
h = L/n
# Local element stiffness and load vectors
Xa = 0
Xb = h
Ke = np.array([[1/h, -1/h], [-1/h, 1/h]])
fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
#-----------------------------------
K_g4 = np.zeros((n-1,n-1))
F_g4 = np.zeros((n-1,1))
K_global(K_g4, n)
F_global(F_g4, n)
U4 = np.linalg.inv(K_g4) @ F_g4
U4 = np.insert(U4, 0, 0, axis = 0) # Boundary condition
U4 = np.insert(U4, n, 0, axis = 0)
F_g4 = np.insert(F_g4, 0, 0, axis = 0)
F_g4 = np.insert(F_g4, n, 0, axis = 0)
Uh4 = 1/2 * np.transpose(U4) @ F_g4
error_2_4 = abs(uex_E_2 - Uh4)
print("\n")
print("-------------------------- Number of elements: 4 --------------------------")
print("Global stiffness matrix(K_g):")
print(K_g4)
print("\nGlobal load vector(F_g):")
print(F_g4)
print("\nThe finite element approximation (uh) coefficients:")
print(U4)
print("\nEnergy in the finite element solution (1/2*uT*F):")
print(Uh4)
print("\n||eE||^2:")
print(error_2_4)
# ------------------------------------------------ n = 6 ------------------------------------------------ #
n = 6
h = L/n
# Local element stiffness and load vectors
Xa = 0
Xb = h
Ke = np.array([[1/h, -1/h], [-1/h, 1/h]])
fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
#------------------------------
K_g6 = np.zeros((n-1,n-1))
F_g6 = np.zeros((n-1,1))
K_global(K_g6, n)
F_global(F_g6, n)
U6 = np.linalg.inv(K_g6) @ F_g6
U6 = np.insert(U6, 0, 0, axis = 0) # Boundary condition
U6 = np.insert(U6, n, 0, axis = 0)
F_g6 = np.insert(F_g6, 0, 0, axis = 0)
F_g6 = np.insert(F_g6, n, 0, axis = 0)
Uh6 = 1/2 * np.transpose(U6) @ F_g6
error_2_6 = abs(uex_E_2 - Uh6)
print("\n")
print("-------------------------- Number of elements: 6 --------------------------")
print("Global stiffness matrix(K_g):")
print(K_g6)
print("\nGlobal load vector(F_g):")
print(F_g6)
print("\nThe finite element approximation (uh) coefficients:")
print(U6)
print("\nEnergy in the finite element solution (1/2*uT*F):")
print(Uh6)
print("\n||eE||^2:")
print(error_2_6)
#---------------------------------------------------- Error ---------------------------------------------------#
error = [error_2_2**0.5, error_2_4**0.5, error_2_6**0.5]
logerror = np.log10(error)
logerror = logerror.ravel()   # Convert to 1D array
dx = [L/2, L/4, L/6]
logdx = -np.log10(dx)
#---------------------------------------------------- Plot ----------------------------------------------------#
plt.figure(1)
plt.plot(x,u_exact,label="Exact Solution",linewidth=1.5)
x2 = np.linspace(0,1,3)
plt.plot(x2,U2,label="Weak Formulation, n = 2",linewidth=0.8)
x4 = np.linspace(0,1,5)
plt.plot(x4,U4,label="Weak Formulation, n = 4",linewidth=0.8)
x6 = np.linspace(0,1,7)
plt.plot(x6,U6,label="Weak Formulation, n = 6",linewidth=0.8)
plt.xlabel("x / L")
plt.ylabel("Finite element approximation (uh)")
plt.title("Comparing Exact Solution with Finite Element Approximation(-u'' = x , u(0) = 0 , u(1) = 0)")
plt.legend()
plt.tight_layout()
plt.show()
#-----------------------------------------------------------
plt.figure(2)
plt.plot(logdx, logerror,linewidth=1.5)
plt.xlabel("-Log h")
plt.ylabel("Log ||error||")
plt.title("Log-log plots of the error as a function of h")
plt.tight_layout()
plt.show()
#-----------------------------------------------------------------------------------------
###-------------- Question 2: -u'' + u = x     |     u(0) = 0 , u(1) = 0 --------------###
#-----------------------------------------------------------------------------------------
#----- Exact Solution -----#
x = np.linspace(0,1,100)
u_exact = x - np.sinh(x)/np.sinh(1)
L = 1
#----- Weak Formulation -----#
# Global K Matrix
def K_global(K_g, n):
    for row in range (1,n):
        if (row == 1):
            K_g[0][0] = Ke[1][1]
        else:
            K_g[row -2][row - 2] += Ke[0][0]
            K_g[row-1][row-1] = Ke[1][1]
            K_g[row -2][row-1] = Ke[0][1]
            K_g[row-1][row-2] = Ke[1][0]
        if(row == n-1):
            K_g[row-1][row-1] += Ke[0][0]
# Global F vector
def F_global(F_g, n):
    for row in range (1,n):
        if (row == 1):
            Xa = 0
            Xb = h
            fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
            F_g[0][0] = fe[1][0]
        else:
            Xa += h
            Xb += h
            fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
            F_g[row-2][0] += fe[0][0]
            F_g[row-1][0] += fe[1][0] 
        if(row == n-1):
            Xa += h
            Xb += h
            fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])            
            F_g[row-1][0] += fe[0][0]
#---------------------------------------------------------------------------------------------------------------------------------------------#
print("--------------------------------------- Question 2: -u'' + u = x     |     u(0) = 0 , u(1) = 0 ---------------------------------------")
#-----------------------------------------------||Uex||^2 -------------------------------------------------#
print("\n")
def d_uex_2(xx):
    return (1 - np.cosh(xx) / np.sinh(1))**2
uex_E_2, _ = spi.quad(d_uex_2, 0, 1)
print("\n||Uex||^2:")
print(uex_E_2)
# ------------------------------------------------ n = 2 ------------------------------------------------ #
n = 2
h = L/n
# Local element stiffness and load vectors
Xa = 0
Xb = h
Ke = np.array([[1/h + h/3, -1/h + h/6], [-1/h + h/6, 1/h + h/3]])
fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
#--------------------------------------
K_g2 = np.zeros((n-1,n-1))
F_g2 = np.zeros((n-1,1))
K_global(K_g2, n)
F_global(F_g2, n)
U2 = np.linalg.inv(K_g2) @ F_g2
U2 = np.insert(U2, 0, 0, axis = 0) # Boundary condition
U2 = np.insert(U2, n, 0, axis = 0)
F_g2 = np.insert(F_g2, 0, 0, axis = 0)
F_g2 = np.insert(F_g2, n, 0, axis = 0)
Uh2 = 1/2 * np.transpose(U2) @ F_g2
error_2_2 = abs(uex_E_2 - Uh2)
print("\n")
print("-------------------------- Number of elements: 2 --------------------------")
print("Global stiffness matrix(K_g):")
print(K_g2)
print("\nGlobal load vector(F_g):")
print(F_g2)
print("\nThe finite element approximation (uh) coefficients:")
print(U2)
print("\nEnergy in the finite element solution (1/2*uT*F):")
print(Uh2)
print("\n||eE||^2:")
print(error_2_2)
# ------------------------------------------------ n = 4 ------------------------------------------------ #
n = 4
h = L/n
# Local element stiffness and load vectors
Xa = 0
Xb = h
Ke = np.array([[1/h + h/3, -1/h + h/6], [-1/h + h/6, 1/h + h/3]])
fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
#-----------------------------------
K_g4 = np.zeros((n-1,n-1))
F_g4 = np.zeros((n-1,1))
K_global(K_g4, n)
F_global(F_g4, n)
U4 = np.linalg.inv(K_g4) @ F_g4
U4 = np.insert(U4, 0, 0, axis = 0) # Boundary condition
U4 = np.insert(U4, n, 0, axis = 0)
F_g4 = np.insert(F_g4, 0, 0, axis = 0)
F_g4 = np.insert(F_g4, n, 0, axis = 0)
Uh4 = 1/2 * np.transpose(U4) @ F_g4
error_2_4 = abs(uex_E_2 - Uh4)
print("\n")
print("-------------------------- Number of elements: 4 --------------------------")
print("Global stiffness matrix(K_g):")
print(K_g4)
print("\nGlobal load vector(F_g):")
print(F_g4)
print("\nThe finite element approximation (uh) coefficients:")
print(U4)
print("\nEnergy in the finite element solution (1/2*uT*F):")
print(Uh4)
print("\n||eE||^2:")
print(error_2_4)
# ------------------------------------------------ n = 6 ------------------------------------------------ #
n = 6
h = L/n
# Local element stiffness and load vectors
Xa = 0
Xb = h
Ke = np.array([[1/h + h/3, -1/h + h/6], [-1/h + h/6, 1/h + h/3]])
fe = h/6 * np.array([[2*Xa + Xb], [Xa + 2*Xb]])
#------------------------------
K_g6 = np.zeros((n-1,n-1))
F_g6 = np.zeros((n-1,1))
K_global(K_g6, n)
F_global(F_g6, n)
U6 = np.linalg.inv(K_g6) @ F_g6
U6 = np.insert(U6, 0, 0, axis = 0) # Boundary condition
U6 = np.insert(U6, n, 0, axis = 0)
F_g6 = np.insert(F_g6, 0, 0, axis = 0)
F_g6 = np.insert(F_g6, n, 0, axis = 0)
Uh6 = 1/2 * np.transpose(U6) @ F_g6
error_2_6 = abs(uex_E_2 - Uh6)
print("\n")
print("-------------------------- Number of elements: 6 --------------------------")
print("Global stiffness matrix(K_g):")
print(K_g6)
print("\nGlobal load vector(F_g):")
print(F_g6)
print("\nThe finite element approximation (uh) coefficients:")
print(U6)
print("\nEnergy in the finite element solution (1/2*uT*F):")
print(Uh6)
print("\n||eE||^2:")
print(error_2_6)
#---------------------------------------------------- Error ---------------------------------------------------#
error = [error_2_2**0.5, error_2_4**0.5, error_2_6**0.5]
logerror = np.log10(error)
# Convert to 1D array
logerror = logerror.ravel()
dx = [L/2, L/4, L/6]
logdx = -np.log10(dx)
#---------------------------------------------------- Plot ----------------------------------------------------#
plt.figure(1)
plt.plot(x,u_exact,label="Exact Solution",linewidth=1.5)
x2 = np.linspace(0,1,3)
plt.plot(x2,U2,label="Weak Formulation, n = 2",linewidth=0.8)
x4 = np.linspace(0,1,5)
plt.plot(x4,U4,label="Weak Formulation, n = 4",linewidth=0.8)
x6 = np.linspace(0,1,7)
plt.plot(x6,U6,label="Weak Formulation, n = 6",linewidth=0.8)
plt.xlabel("x / L")
plt.ylabel("Finite element approximation (uh)")
plt.title("Comparing Exact Solution with Finite Element Approximation(-u'' + u = x , u(0) = 0 , u(1) = 0)")
plt.legend()
plt.tight_layout()
plt.show()
#-----------------------------------------------------------
plt.figure(2)
plt.plot(logdx, logerror,linewidth=1.5)
plt.xlabel("-Log h")
plt.ylabel("Log ||error||")
plt.title("Log-log plots of the error as a function of h")
plt.tight_layout()
plt.show()
