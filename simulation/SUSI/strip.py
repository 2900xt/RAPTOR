# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 19:59:31 2019

@author: lauren
"""
import numpy as np
from susi_utils import peat_hydrol_properties, CWTr, wrc

#from susi_utils import gmeanTr, Hadjacent, Amatrix, boundConst, rightSide


class StripHydrology():
    def __init__(self, spara):
        nLyrs = spara['nLyrs']                                                 # number of soil layers
        dz = np.ones(nLyrs)*spara['dzLyr']                                     # thickness of layers, m
        z = np.cumsum(dz)-dz/2.                                                # depth of the layer center point, m 
        
        # Set default peat type if not specified
        if 'peat type' not in spara or len(spara['peat type']) == 0:
            spara['peat type'] = ['A']  # Default to 'A' (all) type
        if 'peat type bottom' not in spara:
            spara['peat type bottom'] = ['A']  # Default to 'A' (all) type
            
        if spara['vonP']:
            lenvp=len(spara['vonP top'])    
            vonP = np.ones(nLyrs)*spara['vonP bottom'] 
            vonP[0:lenvp] = spara['vonP top']                                      # degree of  decomposition, von Post scale
            ptype = spara['peat type bottom']*spara['nLyrs']
            lenpt = len(spara['peat type']); ptype[0:lenpt] = spara['peat type']    
            self.pF, self.Ksat = peat_hydrol_properties(vonP, var='H', ptype=ptype) # peat hydraulic properties after Päivänen 1973    
        else:
            lenbd=len(spara['bd top'])    
            bd = np.ones(nLyrs)*spara['bd bottom'] 
            bd[0:lenbd] = spara['bd top']                                      # degree of  decomposition, von Post scale
            ptype = spara['peat type bottom']*spara['nLyrs']
            lenpt = len(spara['peat type']); ptype[0:lenpt] = spara['peat type']    
            self.pF, self.Ksat = peat_hydrol_properties(bd, var='bd', ptype=ptype) # peat hydraulic properties after Päivänen 1973    
        
        for n in range(nLyrs):
            if z[n] < 0.41: 
                self.Ksat[n]= self.Ksat[n]*spara['anisotropy']
            else:
                self.Ksat[n]= self.Ksat[n]*1.
                
        self.hToSto, self.stoToGwl, self.hToTra, self.C, self.hToRat, self.hToAfp = CWTr(nLyrs, z, dz, self.pF, 
                                                               self.Ksat, direction='negative') # interpolated storage, transmissivity, diff water capacity, and ratio between aifilled porosoty in rooting zone to total airf porosity  functions
        
        """
        import pandas as pd
        out = {}; out['z'] = []; out['sto']=[]
        for kz in np.arange(0.,2.0, 0.05):
            out['z'].append(-kz)
            out['sto'].append(self.hToSto(-kz))            
        dfout = pd.DataFrame(out)
        dfout.to_excel('soil.xlsx')
        import sys; sys.exit()
        """
        self.L= spara['L']                                                     # compartemnt width, m
        self.n= spara['n']                                                     # number of computation nodes
        self.dy = float(self.L/self.n)                                              # node width m
        sl= spara['slope']                                                     # slope %
        lev=1.                                                                 # basic level of soil surface
        self.ele = np.linspace(0,self.L*sl/100., self.n) + lev                      # surface rise in y direction, m
        self.dt= 1                                                             # time step, days
        self.implic = 1.# 0.5                                                  # 0-forward Euler, 1-backward Euler, 0.5-Crank-Nicolson
        self.DrIrr=False
        self.hini = spara['initial h']                                         # h in the compartment
        self.dwt =np.empty(self.n)                                             # depth to water table, negative down m
        print ('Peat strip initialized')

        #print spara['vonP'], self.Ksat
        #import sys; sys.exit()
        
 
    def reset_domain(self):
        self.A=np.zeros((self.n,self.n))                                       # computation matrix 
        self.h = np.ones(self.n)*self.hini                                     # right hand side vector
        self.H=self.ele+self.h                                                 # head with respect to absolute reference level, m
        self.sruno=0.
        self.roff = 0.
        print ('Resetting strip scenario')


    def run_timestep(self,d,h0ts_west, h0ts_east, p, moss):
        """
        IN: 
            d day number
            h0ts boudary (ditch depth, m) in time series
            p rainfall-et m, arrayn n length
            moss as object
        """
        n= self.n        
        self.h[0]= h0ts_west
        self.h[n-1]= h0ts_east                                            # symmetrical boundaries, set water level
        Htmp = self.H.copy(); Htmp1=self.H.copy()
        #S = p/1000.*np.ones(n)                                                # source/sink, in m
        S = p.copy() #*np.ones(n)                                              # source/sink, in m
        airv = self.hToSto(self.ele)-self.hToSto(Htmp-self.ele)                # air volume, in m
        S=np.where(S>airv, airv, S)                        
        exfil = p - S
        surface_runoff = moss.returnflow(exfil)
        self.sruno += surface_runoff
        #self.sruno += np.where(np.ones(len(airv))*(p)/1000. > airv, np.ones(len(airv))*(p)/1000.-airv, 0.0)  #cut the surface water above to runoff
        Tr0 = self.hToTra(self.H-self.ele)                                     # Transmissivity from the previous time step 
        Trminus0, Trplus0 = self.gmeanTr(Tr0)                                  # geometric mean of adjacent node transmissivities
        Hminus, Hplus = self.Hadjacent(self.H)                                 # vector of adjacent node H

        for it in range(100):                                                  # iteration loop for implicit solution
            Tr1 = self.hToTra(Htmp1-self.ele)                                  # transmissivity in new iteration
            Tr0 = np.maximum(self.hToTra(self.H-self.ele),0.0)               
            Tr1 = np.maximum(self.hToTra(Htmp1-self.ele),0.0)  
            CC = self.C(Htmp1-self.ele)                                        # storage coefficient in new iteration
            Trminus1, Trplus1 = self.gmeanTr(Tr1)                              # geometric mean of adjacent node transmissivity                
            alfa = CC*self.dy**2/self.dt
            self.A = self.Amatrix(self.A, n, self.implic, Trminus1, Trplus1, alfa)              # construct tridiaginal matrix
            self.A = self.boundConst(self.A, n)                                # constant head boundaries to A matrix
            hs=self.rightSide(S,self.dt,self.dy,self.implic,alfa, self.H,Trminus0,Hminus,Trplus0,Hplus,\
                self.DrIrr, Htmp1, self.ele, h0ts_west,h0ts_east)                             # right hand side of the equation
            Htmp1 = np.linalg.multi_dot([np.linalg.inv(self.A),hs])            # solve equation                    
            #Htmp1 = np.matmul(np.linalg.inv(self.A),hs)            # solve equation                    

            Htmp1=np.where(Htmp1>self.ele, self.ele,Htmp1)                     # cut the surface water 
            conv = max(np.abs(Htmp1-Htmp))                                     # define convergence
            Htmp=Htmp1.copy()                                                  # new wt to old for new iteration
            if conv < 1.e-7: 
                if d%365==0: print ('  - day #',d, 'iterations', it)                    
                break
        self.H=Htmp1.copy()                     
        self.roff = self.runoff(self.H, Trminus1, Trplus1, self.dt, self.dy, self.L) + np.average(surface_runoff)       
        self.dwt= self.H-self.ele
        self.air_ratio = self.hToRat(self.dwt)
        self.afp = self.hToAfp(self.dwt)
        return self.dwt, self.H, self.roff, self.air_ratio, self.afp

    def Hadjacent(self,H):
        """
        Input:
            H vector, H in each node
        Output:
            Hwest H(i-1), Heast H(i+1)
        """
        n=len(H)
        Hwest = H[0:n-1]; Hwest=np.append(Hwest, 0.0)
        Heast = H[1:]; Heast=np.insert(Heast, 0, 0.0)
        return Hwest, Heast  

    def Amatrix(self, A, n, implic, Trwest, Treast, alfa):
        """
        Construction of tridiagonal matrix
        """     
        i,j = np.indices(A.shape)
        A[i==j]= implic*(Trwest+Treast)+alfa                                   # diagonal element
        A[i==j+1]=-implic*Trwest[:n-1]                                         # West element
        A[i==j-1]=-implic*Treast[1:]                                           # East element     
    
        return A
    
    def boundConst(self, A, n):
        """
        Diriclet (constant head boundary conditions)
        """    
        A[0,0]=1; A[0,1]=0.                                                    # Dirichlet, west boundary
        A[n-1,n-1]=1.; A[n-1, n-2]=0.                                          # Dirichlet, east boundary
        return A
    
    def boundNoFlow(A, n, implic, Trwest, Treast, alfa):
        """
        Diriclet (constant head boundary conditions)
        """    
                                                                                
        A[0,0]= 2.*implic*(Treast[0])+alfa[0]                                  # Diagonal element
        A[0,1]=-2*implic*Trwest[0]                                             # East element     
        A[n-1,n-1]=2.*implic*(Trwest[n-1])+alfa[0]
        A[n-1, n-2]=-2*implic*Treast[n-1]                                      # West element
    
        return A
    
    
    def rightSide(self, S,dt,dy, implic,alfa, H, Trminus0,Hminus,Trplus0,Hplus, DrIrr, Htmp1, ele, h0_west, h0_east):
        hs = S*dt*dy**2 + alfa*H + (1-implic)*(Trminus0*Hminus) - (1-implic)*(Trminus0 + Trplus0)*H  + (1-implic)*(Trplus0*Hplus)
        n=len(Htmp1)
                
        if DrIrr==False:                
            hs[0]=Htmp1[1] if Htmp1[0]>Htmp1[1] else min(ele[0]+h0_west, Htmp1[1])
            hs[n-1]=Htmp1[n-2] if Htmp1[n-1]>Htmp1[n-2] else min(ele[n-1]+h0_east, Htmp1[n-2])    #if wt below canal water level, lower the canal wl to prevent water inflow to compartment
        else:
            hs[0]=ele[0]+h0_west
            hs[n-1]=ele[n-1]+h0_east 
        return hs
    
    
    def gmeanTr(self, Tr):
        """
        Input: 
            Transmissivity vector, tr in node center point
        Output:
            Transmissivity, tr in west surface sqrt(Tr(i-1)*Tr(i)) and east sqrt(Tr(i)*Tr(i+1)) 
        """
        n = len(Tr)
        Trwest = np.sqrt(Tr[:n-1]*Tr[1:]); Trwest=np.append(Trwest, 0.0)
        Treast = np.sqrt(Tr[1:]*Tr[:n-1]); Treast=np.insert(Treast, 0, 0.0)                
        return Trwest, Treast

    def runoff(self, H, Trminus, Trplus, dt, dy,L):
        roffeast = ((H[1]-H[0])/dy*Trminus[0]*dt)/L
        roffwest = ((H[-2]-H[-1])/dy*Trplus[-1]*dt/L)        
        return roffeast+roffwest

def drain_depth_development(length, hdr, hdr20y):
    """
    Computes daily level of drain bottom thru the time of the simulation. Model adjusted from Hannu Hökkä drain model.
    Input:
        - drain depth in the beginning of the simulation (m, negative down)
        - drain depth after 20 yrs (m, negative down)
        - length of simulation in days
    Output:
        - daily drain bottom level (m, negative down) 
    """    
    timeyrs = np.linspace(0,length/365., length)                                     # time vector telling the time elapsed from the simulation beginning time, yrs        
    h0ts = ((-100*hdr20y +100.*hdr) /np.log(20.0)*np.log(timeyrs+1.)-100*hdr)/-100.  # Ditch model Hökkä
    return h0ts
  