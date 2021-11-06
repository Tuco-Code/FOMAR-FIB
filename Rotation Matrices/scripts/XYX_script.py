import numpy as np
import sympy as sp
import math

sp.init_printing(use_unicode=True)

def alt_sym_euler_rot():
    
    alpha , beta , gamma = sp.symbols(u"alpha beta gamma")

    rotX2 = sp.Matrix([[1,0,0],
                     [0,sp.cos(gamma),-sp.sin(gamma)],
                     [0,sp.sin(gamma),sp.cos(gamma)]])
    
    rotY = sp.Matrix([[sp.cos(beta),0,sp.sin(beta)],
                      [0,1,0],
                      [-sp.sin(beta),0,sp.cos(beta)]])
    
    rotX = sp.Matrix([[1,0,0],
                     [0,sp.cos(alpha),-sp.sin(alpha)],
                     [0,sp.sin(alpha),sp.cos(alpha)]])
    
    return rotX * rotY * rotX2

def rot_xyx2(m11,m31,m33):
    
    try:
        beta = math.acos(m11)
        alpha = math.acos(-m31/math.sin(beta))

        A = -math.sin(alpha)
        B = math.cos(beta) * math.cos(alpha)
        gamma = -math.asin(m33/math.sqrt(A*A + B*B)) - math.atan(B/A)
    except:
        print("CAN'T CONSTRUCT A ROTATION MATRIX FROM THOSE PARAMETERS")
        print("USE VALID M11,M31,M13 VALUES OF A ROTATION MATRIX")
        return
    
    rotXYX = np.array([[math.cos(beta) , math.sin(beta)*math.sin(gamma) , math.cos(gamma)*math.sin(beta) ],
              [math.sin(alpha)*math.sin(beta) , -math.sin(alpha)*math.sin(gamma)*math.cos(beta) + math.cos(alpha)*math.cos(gamma) , -math.sin(alpha)*math.cos(beta)*math.cos(gamma) - math.sin(gamma)*math.cos(alpha) ],
              [-math.sin(beta)*math.cos(alpha) , math.sin(alpha)*math.cos(gamma) + math.sin(gamma)*math.cos(alpha)*math.cos(beta) , -math.sin(alpha)*math.sin(gamma) + math.cos(alpha)*math.cos(beta)*math.cos(gamma) ]])
    
    return np.array([alpha,beta,gamma]),rotXYX

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

angles , mat = rot_xyx2(R[0][0],R[2][0],R[2][2])

angles = angles * 180/math.pi
output = "{} , {} , {}".format(angles[0],angles[1],angles[2])

file = open("fisef.out", "w")
file.write(output)
file.close()
###############################################
