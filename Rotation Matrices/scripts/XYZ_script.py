import numpy as np
import sympy as sp
import math

sp.init_printing(use_unicode=True)

def sym_euler_rot():
    
    alpha , beta , gamma = sp.symbols(u"alpha beta gamma")

    rotX = sp.Matrix([[1,0,0],
                     [0,sp.cos(gamma),-sp.sin(gamma)],
                     [0,sp.sin(gamma),sp.cos(gamma)]])
    
    rotY = sp.Matrix([[sp.cos(beta),0,sp.sin(beta)],
                      [0,1,0],
                      [-sp.sin(beta),0,sp.cos(beta)]])
    
    rotZ = sp.Matrix([[sp.cos(alpha),-sp.sin(alpha),0],
                      [sp.sin(alpha),sp.cos(alpha),0],
                      [0,0,1]])
    
    return rotZ * rotY * rotX

def rot_zyx(m11,m31,m33):
    
    try:
        beta = -math.asin(m31)
        alpha = math.acos(m11/math.cos(beta))
        gamma = math.acos(m33/math.cos(beta))
    except:
        print("CAN'T CONSTRUCT A ROTATION MATRIX FROM THOSE PARAMETERS")
        print("USE VALID M11,M31,M33 VALUES OF A ROTATION MATRIX")
        return
    
    rotZYX = np.array([[math.cos(alpha)*math.cos(beta) , -math.cos(gamma)*math.sin(alpha) + math.cos(alpha)*math.sin(beta)*math.sin(gamma) , math.cos(alpha)*math.cos(gamma)*math.sin(beta) + math.sin(alpha)*math.sin(gamma)],
              [math.cos(beta)*math.sin(alpha) , math.cos(alpha)*math.cos(gamma) + math.sin(alpha)*math.sin(beta)*math.sin(gamma) , math.cos(gamma)*math.sin(alpha)*math.sin(beta) - math.cos(alpha)*math.sin(gamma)],
              [-math.sin(beta) , math.cos(beta)*math.sin(gamma) , math.cos(beta)*math.cos(gamma)]])
    
    return np.array([alpha,beta,gamma]),rotZYX

###############################################
# WRITING
    
try:
    file = open("rijxcolumnes","r")
    text = file.read()
    file.close();
except:
    print("A FILE NAMED 'rijxcolumnes' HAS NOT BEEN FOUND")
    quit()

num_list = text.split()
rij_cols = [float(n[0:-1]) for n in num_list[0:-1]]
rij_cols.append(float(num_list[-1]))

c1 = np.array(rij_cols[0:3])
c2 = np.array(rij_cols[3:6])
c3 = np.array(rij_cols[6:])

R = np.column_stack((c1,c2,c3))

angles , mat = rot_zyx(R[0][0],R[2][0],R[2][2])

angles = angles * 180/math.pi
output = "{} , {} , {}".format(angles[0],angles[1],angles[2])

file = open("fisef.out", "w")
file.write(output)
file.close()
###############################################
