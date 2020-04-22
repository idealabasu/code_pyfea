# -*- coding: utf-8 -*-
"""
copyright 2016-2017 Dan Aukes
"""


import numpy
numpy.set_printoptions(precision = 3)
import matplotlib.cm as cm
import matplotlib.pyplot as plt
plt.ion()
from mpl_toolkits.mplot3d import Axes3D
import scipy.sparse
import scipy.sparse.linalg
import os
import idealab_tools.data_exchange.dat as dat

NDIM = 3
NBOUNDARYNODES = 3

models = ['mechanics']

def volume_force_empty(x):
    volforce = numpy.zeros((x.shape[0],NDIM))
    return volforce

def surface_force_empty(x,n):
    sforce = numpy.zeros((x.shape[0],NDIM))
    return sforce

def u_d_empty(x):
    mm = x.shape[0]
    M = numpy.zeros((NDIM*mm,NDIM))
    W = numpy.zeros((NDIM*mm,1))
    return W,M

def phigrad(A):
    B = numpy.r_[numpy.zeros((1,NDIM)),numpy.eye(NDIM)]
    PhiGrad = numpy.linalg.solve(A,B)
    return PhiGrad

def stiffness_mechanics_3d(vertices,material):
    E = material.E
    nu = material.nu
    mu = material.mu
    Lambda = material.Lambda
    
    if NDIM==3:
        augmented = augment_3d(vertices)
    elif NDIM==2:
        augmented = augment_2d(vertices)
    PhiGrad = phigrad(augmented)
    
    R = numpy.zeros((6,12))
    R[[0,3,4],0::3] = PhiGrad.T
    R[[3,1,5],1::3] = PhiGrad.T
    R[[4,5,2],2::3] = PhiGrad.T

    C = numpy.zeros((6,6))
    C[:3,:3] = Lambda*numpy.ones((3,3))+2*mu*numpy.eye(3)
    C[3:,3:] = mu*numpy.eye(3)
    a = numpy.linalg.det(augmented)
    b = (R.T).dot(C.dot(R))
    stima = a*b/6
    return stima

def stiffness_mechanics_tri_2d(vertices,material):
    E = material.E
    nu = material.nu
    mu = material.mu
    Lambda = material.Lambda
    
    augmented = augment_2d(vertices)
    PhiGrad = phigrad(augmented)
    
    R = numpy.zeros((3,6))
    R[[0,2],0::2] = PhiGrad.T
    R[[2,1],1::2] = PhiGrad.T
    
    C = numpy.zeros((3,3))
    C[:2,:2] = Lambda*numpy.ones((2,2))+2*mu*numpy.eye(2)
    C[2,2] = mu
    a = numpy.linalg.det(augmented)
    b = (R.T).dot(C.dot(R))
    stima = a*b/2
    return stima

def stiffness_mechanics_quad_2d(vertices,material):
    E = material.E
    nu = material.nu
    mu = material.mu
    Lambda = material.Lambda

    R11 = numpy.array([[2,-2,-1,1],[-2,2,1,-1],[-1,1,2,-2],[1,-1,-2,2]])/6
    R12 = numpy.array([[1,1,-1,-1],[-1,-1,1,1],[-1,-1,1,1],[1,1,-1,-1]])/4
    R22 = numpy.array([[2,1,-1,-2],[1,2,-2,-1],[-1,-2,2,1],[-2,-1,1,2]])/6
    
    F1 = numpy.array([vertices[1,:]-vertices[0,:],vertices[3,:]-vertices[0,:]])
    F = numpy.linalg.inv(F1)
    
    L = numpy.array([Lambda + 2*mu,Lambda,mu])
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

def augment_3d(vertices):
    A = numpy.r_[[[1,1,1,1]],vertices.T]
    return A

def augment_2d(vertices):
    A = numpy.r_[[[1,1,1]],vertices.T]
    return A

def analyze(coordinates,elements):
    '''find elements with negative volume'''
    cp = coordinates[elements].sum(1)/4
    v1 = coordinates[elements][:,1]-coordinates[elements][:,0]
    v2 = coordinates[elements][:,2]-coordinates[elements][:,0]
    n = numpy.cross(v1,v2)
    vcp = cp-coordinates[elements][:,0]
    p = (n*vcp).sum(1)
    a = (p<0).nonzero()[0]
    return a
    
def element_max_stress(coordinates,row,u,material):
    Lambda = material.Lambda
    mu = material.mu
    vertices = coordinates[row,:]
    if NDIM==3:
        augmented = augment_3d(vertices)
    elif NDIM==2:
        augmented = augment_2d(vertices)
    PhiGrad = phigrad(augmented)

    CC = 3*numpy.c_[[1,1,1]].dot(numpy.r_[[row]])
    DD = numpy.c_[[0,1,2]].dot(numpy.c_[1,1,1,1])
    U_Grad = u[CC+DD,0].dot(PhiGrad)
    
    SIGMA = Lambda * U_Grad.trace()*numpy.eye(3)+mu*(U_Grad + U_Grad.T)
    F = ((numpy.linalg.eigvals(SIGMA)**2).sum())**.5
    return F

def max_stress(elements,coordinates,u,material):
    num_elements = elements.shape[0]
    C = numpy.zeros((num_elements,1))

    for jj,row in enumerate(elements):
        F = element_max_stress(coordinates,row,u,material)
        C[jj] = F

    return C

def is_empty(array):
    return len(array)==0

def element_to_nodes(coordinates,elements,val):
    mm = elements.shape[0]
    nn = coordinates.shape[0]
    Area = numpy.zeros((mm,1))
    AreaOmega = numpy.zeros((nn,1))
    AvC = numpy.zeros((nn,1))
    
    for jj,row in enumerate(elements):
        vertices = coordinates[row,:]
        
        Area = numpy.linalg.det(augment_3d(vertices))/6
        AreaOmega[row]+=Area
        AvC[row]+=Area*numpy.c_[[1,1,1,1]]*val[jj]

    AvC = AvC/AreaOmega
    AvC = AvC.squeeze()
    return AvC
        
def show(elements,elements4,tris,coordinates,u,material,factor =100):
    C = max_stress(elements,coordinates,u,material)
    AvC = element_to_nodes(coordinates,elements,C)   
   
    cmap = cm.Spectral

    c_face = (AvC[tris]).sum(1)/3
    c_face -= c_face.min()
    c_face /= c_face.max()
    c_face = 1-c_face
    c_face = cmap(c_face)

    xx = factor*u[0::3,0]+coordinates[:,0]
    yy = factor*u[1::3,0]+coordinates[:,1]
    zz = factor*u[2::3,0]+coordinates[:,2]
    xyz = numpy.c_[xx,yy,zz]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ts = ax.plot_trisurf(xx,yy,zz,triangles = tris)
    ts.set_facecolors(c_face)
    plt.show()
    import idealab_tools.matplotlib_tools
    idealab_tools.matplotlib_tools.equal_axes(ax,xyz)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    return ax

def show23_int(elements,tris,coordinates,u,material,factor =100):    
    C = max_stress(elements,coordinates,u,material)
    AvC = element_to_nodes(coordinates,elements,C)   

    cmap = cm.Spectral

    c_face = (AvC[tris]).sum(1)/3
    c_face -= c_face.min()
    c_face /= c_face.max()
    c_face = [cmap(1-item) for item in c_face]    
    c_face = numpy.array(c_face)

    c_vertex = AvC.copy()
    c_vertex -= c_vertex.min()
    c_vertex /= c_vertex.max()
    c_vertex = 1-c_vertex
    c_vertex = cmap(c_vertex)

    xx = factor*u[0::3,0]+coordinates[:,0]
    yy = factor*u[1::3,0]+coordinates[:,1]
    zz = factor*u[2::3,0]+coordinates[:,2]
    xyz = numpy.c_[xx,yy,zz]
    return xyz,tris,c_face,c_vertex

def show2(elements,tris,coordinates,u,material,factor =100):
    xyz,tris,c_face,c_vertex = show23_int(elements,tris,coordinates,u,material,factor)
    import idealab_tools.plot_tris as pt
    pt.plot_tris(xyz,tris,face_colors = c_face)


def show3(elements,tris,coordinates,u,material,factor =100):
    xyz,tris,c_face,c_vertex = show23_int(elements,tris,coordinates,u,material,factor)
    import idealab_tools.plot_tris as pt
    pt.plot_tris(xyz,tris,verts_colors = c_vertex)

def element_to_global_3d(row):
    I = 3*row[[0,0,0,1,1,1,2,2,2,3,3,3]] + numpy.r_[[0,1,2,0,1,2,0,1,2,0,1,2]]
    return I

def element_to_global_2d(row):
    I = 2*row[[0,0,1,1,2,2]] + numpy.r_[[0,1,0,1,0,1]]
    return I

def quad_element_to_global_2d(row):
    I = 2*row[[0,0,1,1,2,2,3,3]] + numpy.r_[[0,1,0,1,0,1,0,1]]
    return I

def local_to_global_3d(I):
    kk,ll = numpy.c_[numpy.meshgrid(I,I)].transpose(0,2,1).reshape(2,-1)
    return kk,ll


def compute(material,coordinates,elements,elements4,neumann,dirichlet_nodes,f,g,u_d):
    mm = coordinates.shape[0]

    A = numpy.zeros((NDIM*mm,NDIM*mm))
#    A = scipy.sparse.lil_matrix((NDIM*mm,NDIM*mm),dtype = float)
#    A = scipy.sparse.csc_matrix((NDIM*mm,NDIM*mm),dtype = float)
    
    volume_forces = numpy.zeros((NDIM*mm,1))
    
    for jj,row in enumerate(elements):
        vertices = coordinates[row,:]
        if NDIM==3:
            I = element_to_global_3d(row)
            STIMAOUT = stiffness_mechanics_3d(vertices,material)
        elif NDIM==2:
            I = element_to_global_2d(row)
            STIMAOUT = stiffness_mechanics_tri_2d(vertices,material)

        indeces = local_to_global_3d(I)
        A[indeces] += STIMAOUT.flatten()

    #    Volume Forces
        if NDIM==3:
            fs  = f(numpy.r_[[vertices.sum(0)/4]]).T
            AAA = augment_3d(vertices)
            BBB = numpy.r_[fs,fs,fs,fs]/4
            CCC = numpy.linalg.det(AAA)/6*BBB;
        elif NDIM==2:
            fs  = f(numpy.r_[[vertices.sum(0)/3]]).T
            AAA = augment_2d(vertices)
            BBB = numpy.r_[fs,fs,fs]/3
            CCC = numpy.linalg.det(AAA)/2*BBB;
        volume_forces[I] += CCC

    if NDIM==2:
        for jj,row in enumerate(elements4):
            vertices = coordinates[row,:]
            I = quad_element_to_global_2d(row)
        
            STIMAOUT = stiffness_mechanics_quad_2d(vertices,material)
            indeces = local_to_global_3d(I)
            A[indeces] += STIMAOUT.flatten()
    
            fs  = f(numpy.r_[[vertices.sum(0)/4]]).T
            AAA = augment_2d(vertices[:-1,:])
            BBB = numpy.r_[fs,fs,fs,fs]/4
            CCC = numpy.linalg.det(AAA)*BBB;
            volume_forces[I] += CCC
        
    neumann_conditions = numpy.zeros((NDIM*mm,1))
    
    if not is_empty(neumann):
        if NDIM==2:
            n1 = coordinates[neumann[:,1],:]-coordinates[neumann[:,0],:]
            n2 = numpy.array([[0,-1],[1,0]])
            n = n1.dot(n2)
        for ii,row in enumerate(neumann):
            if NDIM==3:
                I = NDIM*row[[0,0,0,1,1,1,2,2,2]] + numpy.r_[[0,1,2,0,1,2,0,1,2]]
                n = numpy.cross( coordinates[row[1],:]-coordinates[row[0],:], coordinates[row[2],:]-coordinates[row[0],:])
                n_norm = numpy.linalg.norm(n)
                gm = g(numpy.r_[[coordinates[row,:].sum(0)/3]],n/n_norm).T
                neumann_conditions[I] += n_norm*(numpy.r_[gm,gm,gm]/6)   
            elif NDIM==2:
                I = NDIM*row[[0,0,1,1]] + numpy.r_[[0,1,0,1]]
                n_norm = numpy.linalg.norm(n[[ii]])
                gm = g(numpy.r_[[coordinates[row,:].sum(0)/2]],n[[ii]]/n_norm).T
                if not not gm.nonzero()[0].tolist():
                    jj=1
                neumann_conditions[I] += n_norm*(numpy.r_[gm,gm]/2)   


    
    W,M = u_d(coordinates[dirichlet_nodes])
    nn = W.shape[0]
    B = numpy.zeros((nn,NDIM*mm))
    
    for kk in range(NDIM):
        for ll in range(NDIM):
            MM = M[ll::NDIM,kk]
            AA = numpy.diag(MM)
            qq = NDIM*dirichlet_nodes+kk
            B[ll::NDIM,qq] = AA
    
    mask2 = abs(B).sum(1)!=0
    pp = mask2.sum()

    top = numpy.c_[A,B[mask2,:].T]
    bottom = numpy.c_[B[mask2,:],numpy.zeros((pp,pp))]
    A2 = numpy.r_[top,bottom]
    b = volume_forces+neumann_conditions
    b2 = numpy.r_[b,W[mask2]]


    A2_sparse = scipy.sparse.csc_matrix(A2)
    
    x = numpy.c_[scipy.sparse.linalg.spsolve(A2_sparse, b2)]
    u = x[:NDIM*mm]
    return x,u

def find_used_elements(*args):
    used_elements = numpy.concatenate([item.flatten() for item in args],0)
    used_elements = numpy.sort(numpy.unique(used_elements))
    return used_elements

def coord_reduce(coordinates,used_coords):
    coordinates2 = coordinates[used_coords]
    mapping = dict([(val,ii) for ii,val in enumerate(used_coords)])
    return coordinates2,mapping

def element_reduce(elements,mapping):
    elements2 = [[mapping[item] for item in row] for row in elements]
    elements3 = numpy.array(elements2,dtype = int)
    return elements3


def plot_tetrahedra(coordinates,elements,jj):
    tris = elements[jj,[[0,1,2],[1,2,3],[2,3,0],[3,0,1]]]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ts = ax.plot_trisurf(coordinates[:,0],coordinates[:,1],coordinates[:,2],triangles = tris)
    ts.set_facecolor((1,0,0,.1))
    ts.set_edgecolor((0,0,0,1))
    plt.show()
    
    import idealab_tools.matplotlib_tools
    idealab_tools.matplotlib_tools.equal_axes(ax,coordinates[tris].reshape((-1,3)))

    for ii,item in enumerate(coordinates[elements[jj]]):
        ax.text3D(*item,s=str(ii))

def fix_tet_order(xyz,quad):        
    a=analyze(xyz,quad)
    quad[a] = quad[a][:,(0,2,1,3)]
    a=analyze(xyz,quad)
    return quad

def remove_zero_volume(xyz,quad):
    T = xyz[quad[:,1:]]-xyz[quad[:,0:1]]
    dt = numpy.array([numpy.linalg.det(item) for item in T])
    quad = quad[dt!=0]
    return quad

def plot_triangles(coordinates,triangles_outer):
    import idealab_tools.matplotlib_tools
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(*coordinates.T,triangles = triangles_outer)
    idealab_tools.matplotlib_tools.equal_axes(ax,coordinates)
    plt.show()
    return ax

def plot_triangles_pyqtgraph(coordinates,triangles_outer):
    import idealab_tools.plot_tris as pt
    pt.plot_tris(coordinates,triangles_outer,verts_colors = (1,0,0,1))

def compute_deformation(coordinates,u,factor = 1):
    uu = coordinates[:,0]+factor*u[0::3,0]
    vv = coordinates[:,1]+factor*u[1::3,0]
    ww = coordinates[:,2]+factor*u[2::3,0]
    return numpy.c_[uu,vv,ww]
    
def plot_nodes(coordinates,triangles,u,ax,factor = 100):
    xyz = compute_deformation(coordinates,u,factor)
    ax.plot3D(xyz[triangles,0],xyz[triangles,1],xyz[triangles,2],'ro')
    
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
        if NDIM==3:
            l = E*nu/((1+nu)*(1-2*nu))
        elif NDIM==2:
            l= E*nu/((1+nu)*(1-2*nu))
        return l
    
class ThermalMaterial(object):
    def __init__(self,k):
        self.k = k