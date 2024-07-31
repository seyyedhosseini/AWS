from scipy.special import exp1
import numpy as np
import math
import pandas as pd

def NumInj_build(nWXmin,nWXmax): 
    NumInj = np.zeros(nWXmax,order='F')
    for i in range(nWXmin, nWXmax+1):
        NumInj[i-1] = i*i
    return NumInj


def Distance_Mat_Builder(simprop, rsvrprop):
    nWXmax = simprop.nWXmax
    nWE = simprop.nWE
    XL = rsvrprop.XL
    YL = rsvrprop.YL
    BXL = rsvrprop.BXL
    BYL = rsvrprop.BYL
    Translate4_1 = [34 / 18, 34 / 18, 2 / 3, 4 / 4, 6 / 5, 4 / 6, 6 / 7, 8 / 8, 6 / 9, 8 / 10, 10 / 11, 12 / 12,
                    14 / 13, 16 / 14, 18 / 15, 20 / 16, 22 / 17, 24 / 18, 26 / 19, 28 / 20]
    Translate8_1 = [1, 1 / 2, 2 / 3, 4 / 4, 2 / 5, 4 / 6, 6 / 7, 4 / 8, 6 / 9, 4 / 10, 6 / 11, 8 / 12, 6 / 13, 8 / 14,
                    10 / 15, 8 / 16, 10 / 17, 12 / 18, 10 / 19, 12 / 20]
    Translate8_2 = [3 / 2, 3 / 2, 34 / 18, 34 / 18, 6 / 5, 8 / 6, 10 / 7, 12 / 8, 10 / 9, 12 / 10, 14 / 11, 16 / 12,
                    14 / 13, 16 / 14, 18 / 15, 20 / 16, 22 / 17, 24 / 18, 26 / 19, 28 / 20]
    Translate16_1 = [1, 1 / 2, 10 / 12, 2 / 16, 2 / 5, 4 / 6, 6 / 7, 4 / 8, 6 / 9, 4 / 10, 6 / 11, 4 / 12, 6 / 13,
                     8 / 14, 6 / 15, 8 / 16, 10 / 17, 8 / 18, 10 / 19, 12 / 20]
    Translate16_2 = [3 / 2, 3 / 2, 5 / 3, 4 / 4, 6 / 5, 8 / 6, 10 / 7, 12 / 8, 10 / 9, 12 / 10, 14 / 11, 16 / 12,
                     18 / 13, 20 / 14, 22 / 15, 24 / 16, 22 / 17, 24 / 18, 26 / 19, 28 / 20]
    Translate16_3 = [3 / 2, 3 / 2, 5 / 3, 4 / 4, 6 / 5, 8 / 6, 10 / 7, 12 / 8, 14 / 9, 12 / 10, 14 / 11, 16 / 12,
                     18 / 13, 20 / 14, 22 / 15, 24 / 16, 22 / 17, 24 / 18, 26 / 19, 28 / 20]
    Translate16_4 = [1 / 2, 1 / 2, 8 / 12, 2 / 16, 2 / 5, 4 / 6, 2 / 7, 4 / 8, 6 / 9, 4 / 10, 6 / 11, 4 / 12, 6 / 13,
                     8 / 14, 6 / 15, 8 / 16, 10 / 17, 8 / 18, 10 / 19, 12 / 20]
    nWYmax = nWXmax
    # nWE: Number of Extraction well
    if nWE == 0:
        nWXmin=1
    elif nWE == 4:
        nWXmin = 2
    elif nWE == 8:
        nWXmin = 3
    else:
        nWXmin = 4
        
    nWmax = nWXmax * nWYmax
    nWmin = nWXmin * nWXmin
    nWTmax = nWmax + nWE

    Rwell = np.zeros(shape=(nWTmax, nWTmax, nWXmax),order='F')
    Xwell = np.zeros(shape=(nWTmax, nWXmax),order='F')
    Ywell = np.zeros(shape=(nWTmax, nWXmax),order='F')
    rL1 = np.zeros(shape=(nWTmax, nWXmax),order='F')
    rT1 = np.zeros(shape=(nWTmax, nWXmax),order='F')
    Qsave2 = np.zeros(shape=(nWTmax, nWXmax),order='F')
    Psave2 = np.zeros(shape=(nWTmax, nWXmax),order='F')

    nWXE = 2
    nWYE = nWXE

    """
    Translate1[0:nWXmax] = 0
    Translate2[0:nWXmax] = 0
    Translate3[0:nWXmax] = 0
    Translate1 = [0] * (nWXmax + 1)
    Translate2 = [0] * (nWXmax + 1)
    Translate3 = [0] * (nWXmax + 1)

    Translate1[0] = 1
    Translate1[1] = 3 / 2

    for nWX in range(2, nWXmax + 1):
        Translate1[nWX] = 2 * (nWX - 2) / nWX

    Translate2[0] = 1 / 2
    Translate2[1] = 1 / 2
    Translate2[2] = 5 / 3
    Translate2[3] = 14 / 8

    for nWX in range(4, nWXmax + 1):
        Translate2[nWX] = 2 * (nWX - 4) / nWX
    """
    for nWX in range(1, nWXmax+1):
        nWY = nWX
        for i in range(1, nWY+1):
            for j in range(1, nWX+1):
                Xwell[(i-1)*nWX+j-1, nWX-1] = XL/nWX/2 + (j-1)*XL/nWX - XL/2 + BXL/2
                Ywell[(i-1)*nWX+j-1, nWX-1] = YL/nWY/2 + (i-1)*YL/nWY + BYL/2 - YL/2
        if nWE == 4:
            for i in range(1, nWYE+1):
                for j in range(1, nWXE+1):
                    Xwell[(i-1) * nWXE + j + nWX * nWY - 1, nWX-1] = Translate4_1[nWX-1] * (XL/nWXE / 2+(j-1)*XL/nWXE - XL/2) + BXL/2
                    Ywell[(i-1) * nWXE + j + nWX * nWY - 1, nWX-1] = Translate4_1[nWX-1] * (YL/nWYE / 2+(i-1)*YL/nWYE - YL/2) + BYL/2
        if nWE == 8:
            for i in range(1, nWYE+1):
                for j in range(1, nWXE+1):
                    Xwell[(i-1) * nWXE + j + nWX * nWY - 1, nWX-1] = Translate8_1[nWX-1] * (XL/nWXE / 2+(j - 1)*XL/nWXE-XL / 2) + BXL / 2
                    Ywell[(i-1) * nWXE + j + nWX * nWY - 1, nWX-1] = Translate8_1[nWX-1] * (YL/nWYE / 2+(i - 1)*YL/nWYE-YL / 2) + BYL / 2
                    Xwell[(i-1) * nWXE + j + nWX * nWY + 4 - 1, nWX-1] = Translate8_2[nWX-1] * (XL/nWXE / 2+(j - 1)*XL/nWXE-XL / 2) + BXL / 2
                    Ywell[(i-1) * nWXE + j + nWX * nWY + 4 - 1, nWX-1] = Translate8_2[nWX-1] * (YL/nWYE / 2+(i - 1)*YL/nWYE-YL / 2) + BYL / 2
        if nWE == 16:
            for i in range(1, nWYE + 1):
                for j in range(1, nWXE + 1):
                    Xwell[(i-1) * nWXE + j + nWX * nWY, nWX-1] = Translate16_1[nWX-1] * (XL / nWXE / 2 + (j - 1) * XL / nWXE - XL / 2) + BXL / 2
                    Ywell[(i-1) * nWXE + j + nWX * nWY, nWX-1] = Translate16_1[nWX-1] * (YL / nWYE / 2 + (i - 1) * YL / nWYE - YL / 2) + BYL / 2
                    Xwell[(i-1) * nWXE + j + nWX * nWY + 4, nWX-1] = Translate16_2[nWX-1] * (XL / nWXE / 2 + (j - 1) * XL / nWXE - XL / 2) + BXL / 2
                    Ywell[(i-1) * nWXE + j + nWX * nWY + 4, nWX-1] = Translate16_2[nWX-1] * (YL / nWYE / 2 + (i - 1) * YL / nWYE - YL / 2) + BYL / 2
            for i in range(1, nWYE + 1):
                for j in range(1, nWXE + 1):
                    Xwell[(i - 1) * nWXE + j + nWX * nWY + 8, nWX-1] = Translate16_4[nWX-1] * (XL / nWXE / 2 + (j - 1) * XL / nWXE - XL / 2) + BXL / 2
                    Ywell[(i - 1) * nWXE + j + nWX * nWY + 8, nWX-1] = Translate16_3[nWX-1] * (YL / nWYE / 2 + (i - 1) * YL / nWYE - YL / 2) + BYL / 2
                    Xwell[(i - 1) * nWXE + j + nWX * nWY + 12, nWX-1] = Translate16_3[nWX-1] * (XL / nWXE / 2 + (j - 1) * XL / nWXE - XL / 2) + BXL / 2
                    Ywell[(i - 1) * nWXE + j + nWX * nWY + 12, nWX-1] = Translate16_4[nWX-1] * (YL / nWYE / 2 + (i - 1) * YL / nWYE - YL / 2) + BYL / 2
    for nWX in range(1, nWXmax + 1):
        for i in range(1, nWTmax + 1):
            X = Xwell[i - 1, nWX - 1]
            Y = Ywell[i - 1, nWX - 1]
            for j in range(1, nWTmax + 1):
                Rwell[i - 1, j - 1, nWX - 1] = math.sqrt((X - Xwell[j - 1, nWX - 1]) ** 2 + (Y - Ywell[j - 1, nWX - 1]) ** 2)

    return Xwell, Ywell, Rwell

def A_ConstRate (nWX, simprop, rsvrprop, NumInj_t, rW, chi_BL, chi_dry, 
                     Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, tD,tDE, rED1, rED2, BCValue,geometry):

    # BCValue = 1: closed boundary condition; BCValue = 2: open boundary condition;
    # BCValue = 3: constant pressure boundary condition (To be implemented)

    # line: 1164-1219
    if geometry == 1:
        pass
    else:
        Xwell, Ywell, Rwell = Distance_Mat_Builder(simprop, rsvrprop)
    
    nWE = simprop.nWE
            
#    NumInj_t = int(NumInj[nWX-1])
    Sa = 0.5 * np.log(chi_dry) + (0.5 / F_Lg) * np.log(chi_BL / chi_dry) \
        -0.5 * (Lamda_g / Lamda_w) * np.log(chi_BL / eta_D3) + 0.5 * 0.577215 * (1 - Lamda_g / Lamda_w)
    
    # Construct A matrix
    A = np.zeros(shape = ((NumInj_t+nWE), (NumInj_t+nWE)),order='F')

    # Injectors to Injectors
    for i in range(NumInj_t):
        for j in range(NumInj_t):
                       
#            Sa = 0.5 * np.log(chi_dry[j]) + (0.5 / F_Lg) * np.log(chi_BL[j] / chi_dry[j]) - 0.5 * (
#                        Lamda_g / Lamda_w) * np.log(chi_BL[j] / eta_D3) + 0.5 * 0.577215 * (1 - Lamda_g / Lamda_w)
            if j > i:
                rD = Rwell[i , j , nWX-1] / rW
                if tD <= rD ** 2 / 4 / np.mean(chi_BL):
                    A[i, j] = 0.5 * Lamda_g / Lamda_w * exp1(
                        rD ** 2 / 4 / tD / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
                elif tD >= rD ** 2 / 4 / np.mean(chi_dry):
                    A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD) - 0.5 * exp1(chi_dry[j]) + 1 / 2 / F_Lg * exp1(
                        chi_dry[j] / F_Lg / eta_D2) - 1 / 2 / F_Lg * exp1(
                        chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                        chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
                elif tD < rD ** 2 / 4 / np.mean(chi_dry) and tD > rD ** 2 / 4 / np.mean(chi_BL):
                    A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD / F_Lg / eta_D2) - 0.5 / F_Lg * exp1(
                        chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                        chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
            elif i == j:
                A[j, j] = 0.5 * (math.log(tD) + 0.80908) + Sa[j] - 0.5 * Lamda_g / Lamda_w * exp1(
                    rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                    -rED2 ** 2 / 4 / tD / eta_D3)
            elif i > j:
                A[i, j] = A[j, i].copy()
    # Injectors to Extractors
    for i in range(NumInj_t):
        for j in range(NumInj_t, NumInj_t + nWE):
            rDE = Rwell[i, j, nWX-1] / rW
            if BCValue == 2:
                A[i, j] = 0.5 * exp1(rDE ** 2 / 4 / tDE)
            elif BCValue == 1:
                A[i, j] = (4 * tDE + rDE ** 2) / 2 / rED1 ** 2 - np.log(rDE / rED1) - 3 / 4

    # Extractors to Injectors
    for i in range(NumInj_t, NumInj_t + nWE):
        for j in range(NumInj_t):
            rDE = Rwell[i, j, nWX-1] / rW
            A[i, j] = 0.5 * Lamda_g / Lamda_w * exp1(rDE ** 2 / 4 / tD / eta_D3) - \
                      0.5 * Lamda_g / Lamda_w * exp1(rED1 ** 2 / 4 / tD / eta_D3) + \
                      2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(-rED2 ** 2 / 4 / tD / eta_D3)
    # Extractors to Extractors
    for i in range(NumInj_t, NumInj_t + nWE):
        for j in range(NumInj_t, NumInj_t + nWE):
            rDE = Rwell[i, j, nWX-1] / rW
            if j > i:
                if BCValue == 2:
                    A[i, j] = 0.5 * exp1(rDE ** 2 / 4 / tDE)
                elif BCValue == 1:
                    A[i, j] = (4 * tDE + rDE ** 2) / 2 / rED1 ** 2 - np.log(rDE / rED1) - 3 / 4
            elif i == j:
                if BCValue == 2:
                    A[i, j] = 0.5 * (np.log(tDE) + 0.80908)
                elif BCValue == 1:
                    A[i, j] = 2 * tDE / rED1 ** 2 + np.log(rED1) - 3 / 4
            elif i > j:
                A[i, j] = A[j, i].copy()
                
    return A


def B_ConstRate(A, nWE, NumInj_t, Qinj, Qext, mug, mub, Thickness, k , kra0, krg0, P0, Psave, it):
    # Construct B matrix

    P = np.zeros(NumInj_t+nWE, order='F')
    for i in range(NumInj_t + nWE):
        BInj = 0
        for j in range(NumInj_t):
            BInj += A[i,j]
        BExt = 0
        for j in range(NumInj_t, NumInj_t + nWE):
            BExt += A[i,j]
            
        P[i] = (Qinj * mug / (2 * math.pi * Thickness * k * krg0) * BInj / 1000000 
                - Qext * mub / (2 * math.pi * Thickness * k * kra0) * BExt / 1000000 
                + P0)  # Krg0 to Krs
            
    P_ave = sum(P[:NumInj_t]) / NumInj_t
    
    err = np.linalg.norm(Psave - P, 2)
    Psave = np.copy(P)
    it += 1   
    return P, P_ave, Psave, err, it


def A_Geometry (nWX, nWE, nWI, rW, chi_BL, chi_dry, 
                     Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, tD,tDE, rED1, rED2, BCValue,Rwell):

    # BCValue = 1: closed boundary condition; BCValue = 2: open boundary condition;
    # BCValue = 3: constant pressure boundary condition (To be implemented)


    Sa = 0.5 * np.log(chi_dry) + (0.5 / F_Lg) * np.log(chi_BL / chi_dry) \
        -0.5 * (Lamda_g / Lamda_w) * np.log(chi_BL / eta_D3) + 0.5 * 0.577215 * (1 - Lamda_g / Lamda_w)
    
    # Construct A matrix
    A = np.zeros(shape = ((nWI+nWE), (nWI+nWE)),order='F')

    # Injectors to Injectors
    for i in range(nWI):
        for j in range(nWI):
                       
#            Sa = 0.5 * np.log(chi_dry[j]) + (0.5 / F_Lg) * np.log(chi_BL[j] / chi_dry[j]) - 0.5 * (
#                        Lamda_g / Lamda_w) * np.log(chi_BL[j] / eta_D3) + 0.5 * 0.577215 * (1 - Lamda_g / Lamda_w)
            if j > i:
                rD = Rwell[i , j] / rW
                if tD <= rD ** 2 / 4 / np.mean(chi_BL):
                    A[i, j] = 0.5 * Lamda_g / Lamda_w * exp1(
                        rD ** 2 / 4 / tD / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
                elif tD >= rD ** 2 / 4 / np.mean(chi_dry):
                    A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD) - 0.5 * exp1(chi_dry[j]) + 1 / 2 / F_Lg * exp1(
                        chi_dry[j] / F_Lg / eta_D2) - 1 / 2 / F_Lg * exp1(
                        chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                        chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
                elif tD < rD ** 2 / 4 / np.mean(chi_dry) and tD > rD ** 2 / 4 / np.mean(chi_BL):
                    A[i, j] = 0.5 * exp1(rD ** 2 / 4 / tD / F_Lg / eta_D2) - 0.5 / F_Lg * exp1(
                        chi_BL[j] / F_Lg / eta_D2) + 1 / 2 * Lamda_g / Lamda_w * exp1(
                        chi_BL[j] / eta_D3) - 0.5 * Lamda_g / Lamda_w * exp1(
                        rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                        -rED2 ** 2 / 4 / tD / eta_D3)
            elif i == j:
                A[j, j] = 0.5 * (math.log(tD) + 0.80908) + Sa[j] - 0.5 * Lamda_g / Lamda_w * exp1(
                    rED1 ** 2 / 4 / tD / eta_D3) + 2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(
                    -rED2 ** 2 / 4 / tD / eta_D3)
            elif i > j:
                A[i, j] = A[j, i].copy()
    # Injectors to Extractors
    for i in range(nWI):
        for j in range(nWI, nWI + nWE):
            rDE = Rwell[i, j] / rW
            if BCValue == 2:
                A[i, j] = 0.5 * exp1(rDE ** 2 / 4 / tDE)
            elif BCValue == 1:
                A[i, j] = (4 * tDE + rDE ** 2) / 2 / rED1 ** 2 - np.log(rDE / rED1) - 3 / 4
    # Extractors to Injectors
    for i in range(nWI, nWI + nWE):
        for j in range(nWI):
            rDE = Rwell[i, j] / rW
            A[i, j] = 0.5 * Lamda_g / Lamda_w * exp1(rDE ** 2 / 4 / tD / eta_D3) - \
                      0.5 * Lamda_g / Lamda_w * exp1(rED1 ** 2 / 4 / tD / eta_D3) + \
                      2 * Lamda_g / Lamda_w * eta_D3 * tD / rED2 ** 2 * np.exp(-rED2 ** 2 / 4 / tD / eta_D3)
    # Extractors to Extractors
    for i in range(nWI, nWI + nWE):
        for j in range(nWI, nWI + nWE):
            rDE = Rwell[i, j] / rW
            if j > i:
                if BCValue == 2:
                    A[i, j] = 0.5 * exp1(rDE ** 2 / 4 / tDE)
                elif BCValue == 1:
                    A[i, j] = (4 * tDE + rDE ** 2) / 2 / rED1 ** 2 - np.log(rDE / rED1) - 3 / 4
            elif i == j:
                if BCValue == 2:
                    A[i, j] = 0.5 * (np.log(tDE) + 0.80908)
                elif BCValue == 1:
                    A[i, j] = 2 * tDE / rED1 ** 2 + np.log(rED1) - 3 / 4
            elif i > j:
                A[i, j] = A[j, i].copy()
                
    return A

def B_Geometry(A, nWE, nWI, Qinj, Qext, mug, mub, Thickness, k , kra0, krg0, P0, Psave, it):
    # Construct B matrix

    P = np.zeros(nWI+nWE, order='F')
    for i in range(nWI + nWE):
        BInj = 0
        for j in range(nWI):
            BInj += A[i,j] *Qinj[j]
        BExt = 0
        for j in range(nWI, nWI + nWE):
            BExt += A[i,j] * Qext[j-nWI]
            
        P[i] = mug / (2 * math.pi * Thickness * k * krg0) * BInj / 1000000 \
            - mub / (2 * math.pi * Thickness * k * kra0) * BExt / 1000000 + P0  # Krg0 to Krs
            
    P_ave = sum(P[:nWI]) / nWI
    
    err = np.linalg.norm(Psave - P, 2)
    Psave = np.copy(P)
    it += 1   
    return P, P_ave, Psave, err, it






