
import numpy as np
import math
import streamlit as st
import sys
# Calculate Reservoir Fluid Properties

class Fluid_Prop:



    def __init__(self, p, temp, salinity, Sar, Sgc, m, n, kra0, krg0):
        self.p = p
        self.temp = temp
        self.salinity = salinity
        self.Sar = Sar
        self.Sgc = Sgc
        self.m = m
        self.n = n
        self.kra0 = kra0
        self.krg0 = krg0

    def Spycher(self):
    # Define constants
        const_r = 8.314472 # J/mole/K
        temp = self.temp
        salinity = self.salinity
        temp_k = self.temp + 273.15 # convert to Kelvin
        p_bar = self.p * 10 # convert to bar
        p = self.p
        m_co2 = 44
        m_h2o = 18.016
        m_nacl = 58.448
        kco2g = 10**(1.189+(1.304/10**2)*temp-(5.446/10**5)*temp**2)
        kco2l = 10**(1.169+(1.368/10**2)*temp-(5.380/10**5)*temp**2)
        kh2o = 10**(-2.209+(3.097/10**2)*temp-(1.098/10**4)*temp**2+(2.048/10**7)*temp**3)
        vco2 = 32.6
        vh2o = 18.1
        aco2 = 7.54*10**7 - 4.13*10**4 * temp_k
        bco2 = 27.80
        bh2o = 18.18
        ah2oco2 = 7.89 * 10 ** 7

        a = 1
        b = -(const_r*temp_k/p)
        c = -(const_r*temp_k*bco2/p) + aco2/(p_bar*np.sqrt(temp_k)) - bco2**2
        d = -aco2*bco2/(p_bar*np.sqrt(temp_k))
        vg = 10
        vg_save = vg
        error = 1

        while error > 0.001:
            f = a*vg**3 + b*vg**2 + c*vg + d
            f_prime = 3*a*vg**2 + 2*b*vg + c
            vg = vg - f/f_prime
            error = abs(vg - vg_save)
            vg_save = vg # cm^3/mole

        self.z = p*vg/(const_r*temp_k)
        self.z2 = (p + 0.001)*vg/(const_r*temp_k) # p in MPa
        try: 
            t1 = math.log(vg/(vg-bco2))
            t2 = bh2o/(vg-bco2)
            t3 = (2*ah2oco2/(83.1447*(temp_k)*np.sqrt(temp_k)*bco2))*math.log((vg+bco2)/vg)
            t4 = (aco2*bh2o/(83.1447*(temp_k)*np.sqrt(temp_k)*bco2*bco2*2))*(math.log(((vg+bco2)/vg)-(bco2/(vg+bco2))))
            t5 = math.log(p*vg/const_r/temp_k)
            phi_water = math.exp(t1+t2-t3+t4-t5)
        except (KeyError, ValueError, ZeroDivisionError) as e:
            st.write('### Maybe selected wrong unit system. Please make sure your inputs are consistent with your Unit System selection:exclamation:') 
            sys.exit()
            
        tt1 = math.log(vg/(vg-bco2))
        tt2 = bco2/(vg-bco2)
        tt3 = (2*aco2/(83.1447*(temp_k)*np.sqrt(temp_k)*bco2))*math.log((vg+bco2)/vg)
        tt4 = (aco2/(83.1447*(temp_k)*np.sqrt(temp_k)*bco2))*(math.log((vg+bco2)/vg)-(bco2/(vg+bco2)))
        tt5 = math.log(p*vg/const_r/temp_k)

        phi_co2 = math.exp(tt1+tt2-tt3+tt4-tt5)
        self.wsb = (salinity*m_nacl)/(1000+salinity*m_nacl)
        self.wwb = 1-self.wsb
        self.wsa = self.wsb

        # By DuanSun
        lamda = -0.411370585 + (6.07632013 / 10 ** 4) * temp_k + (97.5347708 / temp_k) + \
                (-0.0237622469 * p_bar / temp_k) + (0.0170656236 * p_bar / (630 - temp_k)) + \
                (1.41335834 / 10 ** 5) * temp_k * math.log(p_bar)
        xi = (3.36389723 / 10 ** 4) + (-1.98298980 / 10 ** 5) * temp_k + \
                (2.12220830 / 10 ** 3) * p_bar / temp_k + (-5.2487330 / 10 ** 3) * p_bar / (630 - temp_k)
        gam = math.exp(2 * lamda * salinity + xi * salinity * salinity)
        a_a = (kh2o / phi_water / p_bar) * math.exp((p_bar - 1) * vh2o / (83.1447 * temp_k))
        b_b = (phi_co2 * p_bar / 55.508 / gam / kco2g) * math.exp(-(p_bar - 1) * vco2 / (83.1447 * temp_k))

        self.yh2o = (1 - b_b) * 55.508 / ((1 / a_a - b_b) * (2 * salinity + 55.508) + 2 * salinity * b_b)
        self.xco2 = b_b * (1 - self.yh2o)


    def CO2Prop(self):
        # Define CO2 properties
        z = self.z
        z2 = self.z2
        xco2 = self.xco2
        yh2o = self.yh2o
        salinity = self.salinity
        p = self.p
        const_r = 8.314472  # Universal gas constant
        temp_k = self.temp + 273.15
        p_bar = self.p * 10
        Tc = 304.25  # K
        p_bar = 73.9  # bar
        p_mpa = 73.9 / 10  # MPa
        m_co2 = 44
        m_h2o = 18.016
        m_nacl = 58.448
        mco2 = xco2 * (2 * salinity + 1000 / m_h2o) / (1 - xco2)
        self.wca = (mco2 * m_co2) / (1000 + salinity * m_nacl + mco2 * m_co2)
        self.wwa = 1 - self.wsa - self.wca
        self.wwg = (yh2o * m_h2o) / (yh2o * m_h2o + (1 - yh2o) * m_co2)
        self.wcg = 1 - self.wwg

        tr = temp_k / Tc
        vc = const_r * Tc / p_mpa
        pr = p / p_mpa
        pr2 = (p + 0.001) / p_mpa  # 0.001MPa = 0.01 Bar
        rho_co2 = (m_co2) / (vc * tr * z / pr) * 1000  # kg / m ^ 3
        self.rhoc = rho_co2
        rhoc = self.rhoc
        self.rhog = self.rhoc
        muc = (16.485 + (0.0094870 * rhoc) ** 2 - (0.0025939 * rhoc) ** 4 + (0.0019815 * rhoc) ** 6) * 1e-6
        self.mug = muc
        rhoc2 = (m_co2) / (vc * tr * z2 / pr2) * 1000  # kg / m ^ 3
        self.cg = (1 / rhoc) * abs(rhoc - rhoc2) / (0.001 * 1e-5)  # 1 / Pa
        return self.mug, self.cg, self.rhog, self.rhoc

    def BrineProp(self):
        # Define brine properties
        global wca
        pstd = 0.101325
        tstd = 15.56
        temp = self.temp
        salinity = self.salinity
        p = self.p

        m_co2 = 44
        m_h2o = 18.016
        m_nacl = 58.448
        arr = 0.01 * (temp - 150)
        muw = 1 / ((5.38 + 3.80 * arr - 0.26 * arr ** 3) * 1000)
        self.mub = muw * (1 + 0.0816 * salinity - 0.0122 * salinity ** 2 + 0.000128 * salinity ** 3 + 0.000629 * temp * (
                1 - np.exp(-0.7 * salinity)))
        mub = self.mub
        self.mua = mub
        mua = self.mua
        rhow = (1 + (1e-6) * (
                    -80 * temp - 3.3 * temp ** 2 + 0.00175 * temp ** 3 + 489 * p - 2 * temp * p + 0.016 * temp ** 2 * p - (
                1.3e-5) * temp ** 3 * p - 0.333 * p ** 2 - 0.002 * temp * p ** 2)) * 1000
        rhowstd = (1 + (1e-6) * (
                -80 * tstd - 3.3 * tstd ** 2 + 0.00175 * tstd ** 3 + 489 * pstd - 2 * tstd * pstd + 0.016 * tstd ** 2 * pstd - (
            1.3e-5) * tstd ** 3 * pstd - 0.333 * pstd ** 2 - 0.002 * tstd * pstd ** 2)) * 1000

        s = salinity * m_nacl * 0.001
        self.rhob = (rhow * 0.001 + s * (0.668 + 0.44 * s + 1e-6 * (
                300 * p - 2400 * p * s + temp * (80 - 3 * temp - 3300 * s - 13 * p + 47 * p * s)))) * 1000
        rhobstd = (rhowstd * 0.001 + s * (0.668 + 0.44 * s + 1e-6 * (
                300 * pstd - 2400 * pstd * s + tstd * (80 - 3 * tstd - 3300 * s - 13 * pstd + 47 * pstd * s)))) * 1000

        drhowdp = 1e-6 * (489 - 2 * temp + 0.016 * temp ** 2 - 1.3e-5 * temp ** 3 - 0.666 * p - 0.004 * temp * p)
        drhobdp = drhowdp + s * 1e-6 * (300 - 2400 * s + temp * (-13 + 47 * s))
        self.cw = drhobdp / self.rhob / 1000
        vphi = 37.51 + (-9.585e-2) * temp + (8.740e-4) * temp ** 2 + (-5.044e-7) * temp ** 3
        rhocAq = m_co2 / vphi * 1000
        self.rhoa = 1 / ((1 - self.wca) / self.rhob + self.wca / rhocAq)

        return self.cw, self.mua, self.mub, rhobstd, self.rhob



    def Cvalues(self):
     
        wca = self.wca
        wcg = self.wcg
        wwa = self.wwa
        wwb = self.wwb
        wwg = self.wwg
        rhoa = self.rhoa
        rhob = self.rhob
        rhoc = self.rhoc
        rhog = self.rhog
        self.cca = wca * rhoa / rhoc
        self.ccg = wcg * rhog / rhoc
        self.cwa = wwa * rhoa / rhoc
        self.cwg = wwg * rhog / rhoc
        self.cwb = wwb * rhob / rhoc
        
        return self.cca, self.ccg, self.cwa, self.cwg, self.cwb


    def PowerLawFun(self,Sg):

        mua = self.mua
        mug = self.mug

        Sar = self.Sar
        Sgc = self.Sgc
        m = self.m
        n = self.n
        kra0 = self.kra0
        krg0= self.krg0

        kra = kra0 * ((1 - Sg - Sar) / (1 - Sgc - Sar)) ** n
        krg = krg0 * ((Sg - Sgc) / (1 - Sgc - Sar)) ** m

        # No-Hysteresis
        if Sg <= Sgc:
            kra = kra0
            krg = 0

        # **********
        # With-Hysteresis
        # Sgc1=0.02;
        # kra(Sg<=Sgc1)=kra0;
        # krg(Sg<=Sgc1)=0;
        # **********

        if Sg >= 1 - Sar:
            kra = 0
            krg = krg0


 #       try:
 #           if mua * krg ==0:
 #               fg = 0
 #           else:
        fg = 1 / (1 + kra / mua * (mug / krg))
 #       except ZeroDivisionError:
 #           fg = 0  # Set fg to 0 if a ZeroDivisionError occurs


 #       try:
 #           if (Sg - Sgc) * (1 - Sg - Sar) == 0:
 #               dfgdSg = 1e38
 #           else:
 #       dfgdSg = fg * (1 - fg) * (m * (1 - Sg - Sar) + n * (Sg - Sgc)) / (Sg - Sgc) / (1 - Sg - Sar)
        dfgdSg = np.where((Sg - Sgc)*(1 - Sg - Sar) == 0, np.inf, fg * (1 - fg) * (m * (1 - Sg - Sar) + n * (Sg - Sgc)) / (Sg - Sgc) / (1 - Sg - Sar))

        ##! HEADS UP - THIS MAY LEAD TO A PROBLEM OF ERROR ACCUMULATION LATER!

 #       except ZeroDivisionError:
 #           dfgdSg = 1e38  # Set dfgdSg to 0 if a ZeroDivisionError occurs?

        return dfgdSg, fg, kra, krg

    def ShockFun(self,Sg,ca,cg):


        # Objective function to be minimised for locations of zT and zL
        dfgdSg,fg, _, _ = Fluid_Prop.PowerLawFun(self,Sg)
        a = ca * (1 - fg) + cg * fg
        g = ca * (1 - Sg) + cg * Sg

#        try:
#            if (g*dfgdSg == 0) or (a == 0):
#                f = 0
#            else:
#        print(dfgdSg)
#        f = abs(1 / (a / g / dfgdSg) - 1)
        f = np.where((a / g / dfgdSg) == 0, np.inf, abs(1 / (a / g / dfgdSg) - 1))
#        except ZeroDivisionError:
#            f = 0  # Set f to 0 if a ZeroDivisionError occurs?
        return f



    def Gvalues(self,Ss,SgT,SgL):
   
        cca = self.cca
        ccg = self.ccg
        cwa = self.cwa
        cwg = self.cwg
        cwb = self.cwb
        
        Gc1 = 1- Ss
        GcT = cca * (1-SgT) + ccg * SgT
        GcL = cca * (1-SgL) + ccg * SgL
        GwT = cwa * (1-SgT) + cwg * SgT
        GwL = cwa * (1-SgL) + cwg * SgL
        Gw3 = cwb
        return Gc1,GcT,GcL,GwT,GwL,Gw3


    def avalues(self,SgT,SgL):
        global mua,mug

        cca = self.cca
        ccg = self.ccg
        cwa = self.cwa
        cwg = self.cwg
        cwb = self.cwb

        dfgdSgT,fgT,kraT,krgT = self.PowerLawFun(SgT)
        dfgdSgL,fgL,kraL,krgL = self.PowerLawFun(SgL)
        ac1=1
        acT=cca*(1-fgT)+ccg*fgT
        acL=cca*(1-fgL)+ccg*fgL
        awT=cwa*(1-fgT)+cwg*fgT
        awL=cwa*(1-fgL)+cwg*fgL
        aw3=cwb
        return ac1,acT,awT,acL,awL,aw3
    # function F1=F1Fun(alpha,z,zE,zG)
    # F1=1/alpha./zE-3/2+log(zE./z)+(z-zG)./zE;
    # F1(zE>0.5615/alpha)=expint(alpha*z(zE>0.5615/alpha));


    def Fronts(self,k,Porosity,cr,SimTime,rW):
        # Fronts function calculates the location of the front between the gas and water phases
        # and the front between the water and oil phases
        import scipy.optimize as opt
        Sar = self.Sar
        Sgc = self.Sgc

        m = self.m
        n = self.n
        kra0 = self.kra0
        krg0 = self.krg0

        #xco2, yh2o, wsa, wsb, wwb, z, z2 = self.Spycher()
#        mug, cg, rhog, rhoc = self.CO2Prop()
        ## rhoc = 1.03 * rhoc; # to match better
#        cw, mua, mub, rhobstd = self.BrineProp()
        #[cca, ccg, cwa, cwg, cwb] = self.Cvalues(rhoa, rhoc, rhog, rhob, wcg, wca, wwa, wwg, wwb)
        ##options = {'ftol': 1e-8, 'xtol': 1e-6, 'maxfev': 500, 'method': 'interior-point'}
        ##result = minimize(fun, x0, method='SLSQP', options=options)


        # Calculate properties at leading shock / needs a good initial guess
        slope1 = 0
        SgLi = -1
        for sat in np.arange(Sgc + 0.001, 1 - Sar, 0.001):
            dfgdSgL, fgL, _, _ = self.PowerLawFun(sat)
            slope2 = (fgL - 0) / (sat - Sgc)
            if slope2 < slope1:
                SgLi = sat
                break
            slope1 = slope2

        if SgLi < 0:
            # display error message
            st.warning(
                'Error: Convergence did not achieve. '
                'It typically happens when the calculated pressures or rates are too low or too high. Inconsistent relative permeability parameters could be a reason for this error as well.')
            st.stop()
        cca,ccg,_,_,_ = self.Cvalues()
        
        SgL = opt.fmin(self.ShockFun, SgLi, args=(cca, ccg), disp=False)
        SgL = SgL.item()

        dfgdSgL, fgL, _, _ = self.PowerLawFun(SgL)

        # Calculate properties at trailing shock / needs a good initial guess

        slope1 = 0
        SgTi = -1
        for sat in np.arange(1 - Sar - 0.001, SgL, -0.001):
            dfgdSgT, fgT, _, _ = self.PowerLawFun(sat)
            dfgdSgT2, fgT2,_, _ = self.PowerLawFun(1 - Sar)
            slope2 = (fgT - (fgT2 + 0.0008)) / (sat - (1 - Sar))
            if slope2 < slope1:
                SgTi = sat
                break
            slope1 = slope2

        if SgTi < 0:
            # display error message
            st.warning('Error: Convergence did not achieve.'
                             'It typically happens when the calcualted pressures or rates are too low or too high.')
            st.stop()

        cwa = self.cwa
        cwg = self.cwg
        SgT = opt.fmin(self.ShockFun, SgTi, args=(cwa, cwg), disp=False)
        SgT = SgT.item()

        if SgT > 1 - Sar:
            SgT = 1 - Sar - 0.0001
        dfgdSgT, fgT, _, _ = self.PowerLawFun(SgT)

        rhos = 2160   # kg/m3
        Ss = self.wsb * self.rhob * (1-SgT)/rhos
        krs = (1-Ss) **3 / (1 + Porosity/(1+Porosity) * Ss **2) #1-Ss, modify
        Gc1, GcT, GcL, GwT, GwL, Gw3 = self.Gvalues(Ss, SgT, SgL)
        ac1, acT, awT, acL, awL, aw3 = self.avalues(SgT, SgL)
        # Shock locations
        # Bg = 1.842 / rhoc
        # Q = Q/rhoc/86400 # m3/s
        qD2 = ac1 * (GwT - 0) / (awT * (Gc1 - GcT) + acT * (GwT - 0))
        qD3 = qD2 * ((awL * (GcL - 0) + acL * (Gw3 - GwL)) / (aw3 * (GcL - 0) + 0 * (Gw3 - GwL)))
        zT = qD2 * dfgdSgT
        zL = qD2 * dfgdSgL

        # Multi - well injection
        cg = self.cg
        cw = self.cw
        mug = self.mug
        mub = self.mub
        mua = self.mua
        tD = k * krg0 * SimTime / (mug * Porosity * (rW ** 2) * (cg + cr))
        tDE = k * SimTime / (mua * Porosity * (rW ** 2) * (cw + cr))
        Save1 = ((zL - zT) * (0.5 * SgT + 0.5 * SgL) + (1 - Ss) * zT) / zL # considering dryout region


        Krg_ave = krg0 * ((Save1 - Sgc) / (1 - Sgc - Sar)) ** n
        Kra_ave = kra0 * ((1 - Save1 - Sar) / (1 - Sgc - Sar)) ** m


        Lamda_g = krg0 / mug # end point
        Lamda_w = kra0 / mub # end point


        F_Lg = ((Kra_ave / mub) + (Krg_ave / mug)) / Lamda_g
        eta_D2 = (cg + cr) / (cr + cg * Save1 + cw * (1 - Save1))
        eta_D3 = (cg + cr) * Lamda_w / (cw + cr) / Lamda_g
        


        return tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL



"""
#TESTING
Fluid = Fluid_Prop(p=20, temp=65, salinity=2, Sar=0.5, Sgc=0.1,m=3, n=3, kra0=1, krg0=0.3)
Fluid.Spycher()
mug, cg, rhog, rhoc = Fluid.CO2Prop()


cw, mua, mub, rhobstd, rhob = Fluid.BrineProp()
cca,ccg,cwa,cwg,cwb= Fluid.Cvalues()


f=Fluid.ShockFun(Sg=0.3,ca=cwa,cg=cwg)
print(f)

Gc1,GcT,GcL,GwT,GwL,Gw3 = Fluid.Gvalues(Ss=0.5,SgT=0.3,SgL=0.3)
#print(Gc1,GcT,GcL,GwT,GwL,Gw3)


ac1,acT,awT,acL,awL,aw3 = Fluid.avalues(SgT=0.3,SgL=0.3)
#print(ac1,acT,awT,acL,awL,aw3)



tD, tDE,Lamda_g, Lamda_w, F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL= Fluid.Fronts(k=100* 9.869e-16,Porosity=0.2,cr=5e-10,SimTime=20* 365*24*60**2,rW=0.1)
print( F_Lg, eta_D2, eta_D3, dfgdSgT, dfgdSgL, zT, zL)
"""

