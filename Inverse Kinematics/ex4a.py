#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 20:09:07 2021

@author: ivan
"""

import numpy as np
import os

#CREATION OF THE T,D MATRIX THROUGH DH PARAMETERS
#NEEDED TO BE WORKING ON A FUCTION WITH 4 PARAMETRES
# A, ALPHA, D, THETA
def construction(a,al,d,the):
    alpha = np.deg2rad(al)
    theta = np.deg2rad(the)
    t = np.array([[np.cos(theta),-np.sin(theta),0,a] , [np.cos(alpha)*np.sin(theta),np.cos(alpha)*np.cos(theta),-np.sin(alpha),-d*np.sin(alpha)] , [np.sin(alpha)*np.sin(theta), np.sin(alpha)*np.cos(theta),np.cos(alpha),d*np.cos(alpha)] , [0,0,0,1]])
    dm = np.array([[-np.sin(theta),-np.cos(theta),0,0] , [np.cos(alpha)*np.cos(theta),-np.cos(alpha)*np.sin(theta),0,0] , [np.sin(alpha)*np.cos(theta),-np.sin(alpha)*np.sin(theta),0,0] , [0,0,0,0]])
    return t,dm

def jacobian(t,d):
    
    dx1 = d[0] @ t[1] @ t[2] @ t[3] @ t[4] @ t[5] @ t[6] @ t[7] @ t[8] @ t[9]
    dx2 = t[0] @ d[1] @ t[2] @ t[3] @ t[4] @ t[5] @ t[6] @ t[7] @ t[8] @ t[9]
    dx3 = t[0] @ t[1] @ d[2] @ t[3] @ t[4] @ t[5] @ t[6] @ t[7] @ t[8] @ t[9]
    dx4 = t[0] @ t[1] @ t[2] @ d[3] @ t[4] @ t[5] @ t[6] @ t[7] @ t[8] @ t[9]
    #dx5 = t[0] @ t[1] @ t[2] @ t[3] @ d[4] @ t[5] @ t[6] @ t[7] @ t[8] @ t[9]
    dx6 = t[0] @ t[1] @ t[2] @ t[3] @ t[4] @ d[5] @ t[6] @ t[7] @ t[8] @ t[9]
    dx7 = t[0] @ t[1] @ t[2] @ t[3] @ t[4] @ t[5] @ d[6] @ t[7] @ t[8] @ t[9]
    dx8 = t[0] @ t[1] @ t[2] @ t[3] @ t[4] @ t[5] @ t[6] @ d[7] @ t[8] @ t[9]
    dx9 = t[0] @ t[1] @ t[2] @ t[3] @ t[4] @ t[5] @ t[6] @ t[7] @ d[8] @ t[9]
    
    j = np.array([[dx1[0][3],dx2[0][3],dx3[0][3],dx4[0][3],0,dx6[0][3],dx7[0][3],dx8[0][3],dx9[0][3]],
                  [dx1[1][3],dx2[1][3],dx3[1][3],dx4[1][3],0,dx6[1][3],dx7[1][3],dx8[1][3],dx9[1][3]] ,
                  [dx1[0][0],dx2[0][0],dx3[0][0],dx4[0][0],0,dx6[0][0],dx7[0][0],dx8[0][0],dx9[0][0]] ])
    return j

tabla = open("taula-DH","r")
OG = tabla.read(1000)

tabla.close()

for frame in range(0,90):
    actual_x = 13 + (4-13)*(frame/90)
    next_x = 13 + (4-13)*((frame+1)/90)
    
    dx = next_x - actual_x
    dy = 0
    dnx = 0
    
    tabla = open("taula-DH","r")
    sts = tabla.read(1000)
    sts = sts.split(',')
    st = []
    for i in range(0,len(sts)):
        sts[i] = sts[i].strip()
        st.append(float(sts[i]))
    
    tabla.close()
    
    tm = []
    d = []
    for i in range(0,9):
       tmaux,daux = construction(st[(i*4)+1],st[(i*4)],st[(i*4)+2],st[(i*4)+3])
       tm.append(tmaux)
       d.append(daux)
    
    tm.append(np.array([[1,0,0,1.56] , [0,1,0,0], [0,0,1,0], [0,0,0,1]]))
    jac = jacobian(tm,d)
    pseudo_inverse_jac = np.linalg.pinv(jac)
    
    dtheta = np.rad2deg(pseudo_inverse_jac @ np.array([[dx], [dy], [dnx]]))
    
    for i in range(0,8):
        st[(i*4)+3] = float(st[(i*4)+3] + dtheta[i])
    
    tabla = open("taula-DH","w")
    
    dh_data = []
    for data in range(0,36,4):
        aux = [str(st[data]),str(st[data+1]),str(st[data+2]),str(st[data+3])]
        aux = ",".join(aux)
        dh_data.append(aux)
    
    dh_data = ",\n".join(dh_data)
    
    tabla.write(dh_data)
    
    tabla.close()
    
    if frame < 10:
        command = "povray +O./imagesEX4a/0{}frame jcb.pov".format(frame)
    else:
        command = "povray +O./imagesEX4a/{}frame jcb.pov".format(frame)
    
    os.system(command)
    
    

tabla = open("taula-DH","w")
tabla.write(OG)

tabla.close()