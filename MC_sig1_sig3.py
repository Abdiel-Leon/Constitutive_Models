#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 17:04:45 2023
@author: abdiel

### Simple Implementation of MC model for Tx compression test in terms of sig1 and sig3###
### benchmark is correct##
"""
import numpy
import math
import matplotlib.pyplot as plt

class MC_pq:

    def __init__(self, *args):
        # confinement pressure
        self.sig3 = args[0]
        self.sig1 = args[1]
        self.eps1 = 0        
        # Material parameters
        self.phi = args[2]
        self.c = args[3]
        self.E = args[4]        
        self.v = args[5]        
        self.K = self.E/(3*(1-2*self.v))
        self.G = self.E/(2*(1+self.v))        
        self.eta = 6*numpy.sin(self.phi)/(3-numpy.sin(self.phi))


    def Initialize(self):
        tol_F = 1e-6
        self.dF_sig1 = 0
        d_lambda = 0
        d_eps1_p = 0
        d_eps3_p = 0
        return tol_F, d_lambda, d_eps1_p,d_eps3_p

    def Compute_F(self):
        # yield surface
        tol_F,d_lambda,d_eps1_p,d_eps3_p = self.Initialize()
        
        F = -(self.sig1 - self.sig3) + (self.sig1 + self.sig3)*numpy.sin(self.phi) + 2*self.c*numpy.cos(self.phi)


        if F>tol_F:
            # compute plastic multiplier
            print("PLASTICITYYY")
            d_eps1_p = self.Compute_Lambda(F)
        return F, d_eps1_p#, d_eps3_p
    
    def Compute_Lambda(self,F):
        # Plastic multiplier
        self.dF_sig1 = 1-numpy.sin(self.phi)
        self.dF_sig3 = -1-numpy.sin(self.phi)
        
        d_lambda = self.d_eps1
        d_lambda = d_lambda/(self.dF_sig1)
        
        d_eps1_p = d_lambda*self.dF_sig1
       # d_eps3_p = d_lambda*self.dF_sig3
        
        #print((self.sig1 - self.sig3))
        return d_eps1_p#, d_eps3_p 

    def Update_Stress(self,d_eps1_p,d_eps1):
     #   print(d_eps1_p)
        self.dsig1 = self.E*(self.d_eps1 - d_eps1_p)    
        print(self.d_eps1 - d_eps1_p)
        self.sig1 = self.sig1 + self.dsig1
        q = -(self.sig1)#-self.sig3)
        return q
    
   # def Compute_eps_v(self,eps1,eps3,d_eps3_p):
    #     eps3 = eps3 + d_eps3_p
     #    ev = eps1 -self.v*eps1
      #   return ev       
        
    def Compute_test(self):
        n = 2000
        q = numpy.zeros(n)   
        e = numpy.zeros(n)   
        ev = numpy.zeros(n)   

        self.d_eps1 = -0.001/n        
        eps1 =0; #eps3 =0
        for i in range(n):
            #strain increment
            eps1 = eps1 + self.d_eps1
            e[i] = eps1
           # compute yield function
            F, d_eps1_p = self.Compute_F()
            # compute shear stress q
            q[i] = self.Update_Stress(d_eps1_p,self.d_eps1)
            #ev[i] = self.Compute_eps_v(eps1, eps3, d_eps3_p)
        return q,e,ev
    
def main():
 sig3 = -1e3
 sig1_0 = sig3
 phi = 30
 phi = phi*(math.pi/180)
 E_modulus = 1e7
 v_poisson = 0.25
 c = -1e3
 args = (sig3,sig1_0,phi,c,E_modulus,v_poisson)

 Tx_test = MC_pq(*args)
 q,e,ev = Tx_test.Compute_test()
 
# plot results
 plt.ylabel('sig1')
 plt.xlabel('eps1')

 plt.plot(e, q)
 #plt.plot(e, ev)

main()