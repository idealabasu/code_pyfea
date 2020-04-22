# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 15:54:17 2020

@author: danaukes
"""
import numpy
import pyfea.fea as fea

class Material(object):
    def __init__(self,E,nu):
        self.E = E
        self.nu = nu
        self.mu = self.compute_mu(E,nu)
        self.Lambda = self.compute_lambda(E,nu)

    @staticmethod    
    def compute_mu(E,nu):
        mu = E/(2*(1+nu))
        return mu

    @staticmethod    
    def compute_lambda(E,nu):
        if fea.NDIM==3:
            l = E*nu/((1+nu)*(1-2*nu))
        elif fea.NDIM==2:
            l= E*nu/((1+nu)*(1-2*nu))
        return l

    def stiffness_mechanics_3d(self,vertices):
#        E = material.E
#        nu = material.nu
#        mu = material.mu
#        Lambda = material.Lambda
        
        C = numpy.zeros((6,6))
        C[:3,:3] = self.Lambda*numpy.ones((3,3))+2*self.mu*numpy.eye(3)
        C[3:,3:] = self.mu*numpy.eye(3)


        if fea.NDIM==3:
            augmented = fea.augment_3d(vertices)
        elif fea.NDIM==2:
            augmented = fea.augment_2d(vertices)
        PhiGrad = fea.phigrad(augmented)
        
        R = numpy.zeros((6,12))
        R[[0,3,4],0::3] = PhiGrad.T
        R[[3,1,5],1::3] = PhiGrad.T
        R[[4,5,2],2::3] = PhiGrad.T
    
        a = numpy.linalg.det(augmented)
        b = (R.T).dot(C.dot(R))
        stima = a*b/6
        return stima
    
    def stiffness_mechanics_tri_2d(self,vertices):
#        E = material.E
#        nu = material.nu
#        mu = material.mu
#        Lambda = material.Lambda
        
        
        C = numpy.zeros((3,3))
        C[:2,:2] = self.Lambda*numpy.ones((2,2))+2*self.mu*numpy.eye(2)
        C[2,2] = self.mu

        augmented = fea.augment_2d(vertices)
        PhiGrad = fea.phigrad(augmented)
        
        R = numpy.zeros((3,6))
        R[[0,2],0::2] = PhiGrad.T
        R[[2,1],1::2] = PhiGrad.T

        a = numpy.linalg.det(augmented)
        b = (R.T).dot(C.dot(R))
        stima = a*b/2
        return stima
    
    def stiffness_mechanics_quad_2d(self,vertices):
#        E = material.E
#        nu = material.nu
#        mu = material.mu
#        Lambda = material.Lambda
    
        R11 = numpy.array([[2,-2,-1,1],[-2,2,1,-1],[-1,1,2,-2],[1,-1,-2,2]])/6
        R12 = numpy.array([[1,1,-1,-1],[-1,-1,1,1],[-1,-1,1,1],[1,1,-1,-1]])/4
        R22 = numpy.array([[2,1,-1,-2],[1,2,-2,-1],[-1,-2,2,1],[-2,-1,1,2]])/6
        
        F1 = numpy.array([vertices[1,:]-vertices[0,:],vertices[3,:]-vertices[0,:]])
        F = numpy.linalg.inv(F1)
        
        L = numpy.array([self.Lambda + 2*self.mu,self.Lambda,self.mu])
        stima4 = numpy.zeros((8,8))
    
        Lmod = numpy.array([[L[0],0],[0,L[2]]])
        E = (F.T.dot(Lmod)).dot(F)
        stima4[::2,::2] = E[0,0]*R11+E[0,1]*R12 + E[1,0]*(R12.T)+E[1,1]*R22
    
        Lmod = numpy.array([[L[2],0],[0,L[0]]])
        E = (F.T.dot(Lmod)).dot(F)
        stima4[1::2,1::2] = E[0,0]*R11+E[0,1]*R12 + E[1,0]*(R12.T)+E[1,1]*R22
    
        Lmod = numpy.array([[0,L[2]],[L[1],0]])
        E = (F.T.dot(Lmod)).dot(F)
        stima4[1::2,::2] = E[0,0]*R11+E[0,1]*R12 + E[1,0]*(R12.T)+E[1,1]*R22
        stima4[::2,1::2] = stima4[1::2,::2].T
        stima4/=numpy.linalg.det(F)
        return stima4
    
class ThermalMaterial(object):
    def __init__(self,k):
        self.k = k