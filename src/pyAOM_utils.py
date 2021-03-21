import os
import numpy as np
import math
import pandas as pd
import sys
import ast
import tarfile
import time
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
from mulliken import *
from anIres import *

class single_molecule:
    def __init__(self,filename,src='file',input_var=None):
        if src=='file':
            self.atoms,self.species,self.x,self.y,self.z=read_xyz(filename)
        elif src=='variable':
            self.species,self.x,self.y,self.z=input_var
            self.atoms=len(self.species)
        self.unique_species=list(set(self.species))
        self.dist_unit='Ang'
        self.AOM_dict={}
        self.orb_compl_dict={}
        self.state={}
    def angtobohr(self):
        AngToBohr=1.8897259886
        if self.dist_unit=='Ang':
            self.x=[i*AngToBohr for i in self.x]
            self.y=[i*AngToBohr for i in self.y]
            self.z=[i*AngToBohr for i in self.z]
            self.dist_unit='Bohr'
    def bohrtoang(self):
        AngToBohr=1.8897259886
        if self.dist_unit=='Bohr':
            self.x=[i/AngToBohr for i in self.x]
            self.y=[i/AngToBohr for i in self.y]
            self.z=[i/AngToBohr for i in self.z]
            self.dist_unit='Ang'
    def calculate_geom_center(self):
        self.X,self.Y,self.Z=[np.array(self.x).mean(),np.array(self.y).mean(),np.array(self.z).mean()]
    def recenter(self,offset=4.0):
        self.calculate_geom_center()
        R=[max(self.x)-min(self.x)+2.0*offset,max(self.y)-min(self.y)+2.0*offset,max(self.z)-min(self.z)+2.0*offset]
        self.x=[i-self.X+R[0]/2 for i in self.x]
        self.y=[i-self.Y+R[1]/2 for i in self.y]
        self.z=[i-self.Z+R[2]/2 for i in self.z]
        self.calculate_geom_center()
        self.supercell=[R[i] for i in range(3)]
    def prep_cp2k_single(self,configuration_dict,label,template,basis,potential):
        #
        self.recenter()
        # read template
        fp=open(template,mode='r')
        qs=fp.readlines()
        fp.close()
        # apply dict values
        for key,value in configuration_dict.items():
            target=f'_{key}_'
            if isinstance(value,str)==True:
                t_value=value
            else:
                t_value=str(value)
            for c,i in enumerate(qs):
                if i.find(target)!=-1:
                    qs[c]=qs[c].replace(target,t_value)
        # write subsys
        subsys=['&SUBSYS\n',' &CELL\n',f'  ABC {self.supercell[0]} {self.supercell[1]} {self.supercell[2]}\n','  ALPHA_BETA_GAMMA 90.0 90.0 90.0\n','  PERIODIC none\n',' &END CELL\n']
        for i in self.unique_species:
            subsys.append(f' &KIND {i}\n')
            subsys.append(f'  BASIS_SET {basis}\n')
            subsys.append(f'  POTENTIAL {potential}\n')
            subsys.append(f' &END KIND\n')
        subsys.append(' &COORD\n')
        if self.dist_unit!='Ang':
            self.bohrtoang()
        for c,i in enumerate(self.species):
            subsys.append(f'  {i} {self.x[c]} {self.y[c]} {self.z[c]}\n')
        subsys.append(' &END COORD\n')
        subsys.append('&END SUBSYS\n')
        target='#subsys\n'
        pos=[c for c,i in enumerate(qs) if i.find(target)!=-1][0]+1
        qs=qs[:pos:]+subsys+qs[pos::]
        # write to file
        if os.path.exists('output')==False:
            os.mkdir('output')
        path=f'output/{label}'
        if os.path.exists(path)==False:
            os.mkdir(path)
        with open(f'{path}/{configuration_dict["PROJECT_NAME"]}.inp',mode='w') as fp:
            for i in qs:
                print(i,end='',file=fp)
        print(f'Generated {path}/{configuration_dict["PROJECT_NAME"]}.inp')
    def get_cp2k_info(self,MO,cp2k_output_file,cp2k_basis_file,basis):
        self.MO=MO
        self.MOcoeffs=get_cp2k_MO(cp2k_output_file,MO)
        self.basis_dict=read_basis(cp2k_basis_file,basis,self.unique_species)
        self.pcoeff,self.palpha,self.pqn,self.bfnPerAtom,self.GTO_depth=read_CP2K_GTOs(self.species,self.basis_dict)
    def initialize_STOs(self,STO_dict):
        self.angtobohr()
        self.STOs,self.STO_id_array,self.STO_type_array,self.STO_mu_array=initialize_STOs(self.species,self.x,self.y,self.z,STO_dict)
        self.Smatrix=calculate_overlap_S_matrix(self.x,self.y,self.z,self.STOs,self.STO_id_array,self.STO_type_array,self.STO_mu_array)
    def resolve_pvecs(self,STO_matrix=None):
        """Calculate normal vectors. If the STO matrix is given as input, calculate pi projection coeffs"""
        # convert ang to bohr
        self.angtobohr()
        # initialize normal vector components list
        self.px,self.py,self.pz=[np.zeros(self.atoms) for i in range(3)]
        aoinum=0
        atomlist=[]
        for j,i in enumerate(self.species):
            if i!='H':
                aoinum+=1
                atomlist.append(j)
        for i in range(len(atomlist)):
            neighbourlist=[[0 for j in range(4)] for i in range(4)]
            j=0
            for k in range(self.atoms):
                if atomlist[i]!=k:
                    r=np.linalg.norm(np.array([self.x[atomlist[i]],self.y[atomlist[i]],self.z[atomlist[i]]])-np.array([self.x[k],self.y[k],self.z[k]]))
                    if r<3.5:
                        neighbourlist[0][j]=k+1
                        neighbourlist[1][j]=self.x[k]
                        neighbourlist[2][j]=self.y[k]
                        neighbourlist[3][j]=self.z[k]
                        j+=1
            if j==2:
                neighbourlist[0][2]=atomlist[i]+1
                neighbourlist[1][2]=self.x[atomlist[i]]
                neighbourlist[2][2]=self.y[atomlist[i]]
                neighbourlist[3][2]=self.z[atomlist[i]]
            if j==1:
                # use as a second column entry the only 1-2 neighbor
                neighbourlist[0][1]=atomlist[i]+1
                neighbourlist[1][1]=self.x[atomlist[i]]
                neighbourlist[2][1]=self.y[atomlist[i]]
                neighbourlist[3][1]=self.z[atomlist[i]]
                for k in range(self.atoms):
                    if neighbourlist[0][0]!=k+1 and neighbourlist[0][1]!=k+1:
                        r=np.linalg.norm(np.array([self.x[neighbourlist[0][0]-1],self.y[neighbourlist[0][0]-1],self.z[neighbourlist[0][0]-1]])-np.array([self.x[k],self.y[k],self.z[k]]))
                        if r<3.5:
                            neighbourlist[0][2]=k+1
                            neighbourlist[1][2]=self.x[k]
                            neighbourlist[2][2]=self.y[k]
                            neighbourlist[3][2]=self.z[k]
                            break
            veca_x,veca_y,veca_z=[neighbourlist[i][0]-neighbourlist[i][2] for i in [1,2,3]]
            vecb_x,vecb_y,vecb_z=[neighbourlist[i][1]-neighbourlist[i][2] for i in [1,2,3]]
            u=np.cross([veca_x,veca_y,veca_z],[vecb_x,vecb_y,vecb_z])
            u/=np.linalg.norm(u)
            self.px[atomlist[i]],self.py[atomlist[i]],self.pz[atomlist[i]]=u
        if STO_matrix is not None:
            return [self.px[j]*i[1]+self.py[j]*i[2]+self.pz[j]*i[3] for j,i in enumerate(STO_matrix)]
    def project(self):
        self.STO_matrix,self.orb_compl,self.V_array=STO_GTO_projection(self.x,self.y,self.z,
                                        self.STO_id_array,self.STO_type_array,self.STO_mu_array,
                                        self.pcoeff,self.palpha,self.pqn,self.bfnPerAtom,self.GTO_depth,self.MOcoeffs,self.Smatrix)
        self.resolve_pvecs()
        self.AOM_pi_coeffs=[self.px[c]*i[1]+self.py[c]*i[2]+self.pz[c]*i[3] for c,i in enumerate(self.STO_matrix)]
        self.AOM_dict[self.MO]=self.AOM_pi_coeffs
        self.orb_compl_dict[self.MO]=self.orb_compl
    def cube_preview(self,STO_dict,filename):
        create_cube_file(self.species,self.x,self.y,self.z,self.STO_matrix,STO_dict,filename)
    def save_AOM(self,label):
        with open(f'output/{label}/AOM_COEFF.dat',mode='w') as fp:
            print(self.AOM_dict,file=fp)
    def save_state(self,label,to_file=True):
        self.state['pvecs']={'px':list(self.px),'py':list(self.py),'pz':list(self.pz)}
        self.state['S_matrix']=[[j for j in i] for i in self.Smatrix]
        self.state['AOM_dict']=self.AOM_dict
        self.state['V_array']=self.V_array
        self.state['compl_dict']=self.orb_compl_dict
        self.state['STO_matrix']=self.STO_matrix
        _,s,_=np.linalg.svd(self.Smatrix)
        self.state['singular_values']=list(s)
        if to_file==True:
            with open(f'output/{label}/state.dat',mode='w') as fp:
                print(self.state,file=fp)

def read_xyz(filename):
    fp=open(filename,mode='r')
    xyz=fp.readlines()
    fp.close()
    x,y,z=[[float(i.split()[j]) for i in xyz[2::]] for j in [1,2,3]]
    species=[i.split()[0] for i in xyz[2::]]
    return len(species),species,x,y,z
def get_cp2k_MO(filename,MO,report=False):
    # open cp2k output file and store MO info: alpha channel
    fp=open(filename,mode='r')
    cp2k_out=fp.readlines()
    fp.close()
    cp2k_out=[i for i in cp2k_out if i!='\n']
    cart_basis_functions=[int(i.strip().split(':')[-1]) for i in cp2k_out if i.find('Cartesian basis functions')!=-1][0]
    orbitals=[int(i.strip().split(':')[-1]) for i in cp2k_out if i.find('Number of molecular orbitals')!=-1]
    header=[j for j,i in enumerate(cp2k_out) if i.find('MO EIGENVALUES')!=-1]
    blocks=math.ceil(orbitals[0]/4)
    alphaMOs=cp2k_out[header[0]+1:header[0]+1+blocks*(cart_basis_functions+3)]
    alphaMOs=[i.strip().split() for i in alphaMOs]
    frames=[[] for i in range(blocks)]
    basis_info_block=[['MO'],['Energy'],['Occupation']]
    for i in range(blocks):
        start=3+i*(cart_basis_functions+3)
        stop=start+cart_basis_functions-1
        for j in range(3):
            frames[i].append(alphaMOs[start-3+j].copy())
        for j in range(start,stop+1):
            c=alphaMOs[j].copy()
            if i==0:
                basis_info_block.append(' '.join(c[0:4]))
            del c[0:4]
            frames[i].append(c)
    [i for i in frames[0]]
    alpha_records=[]
    for i in range(3,len(frames[0])):
        line=[]
        for j in range(blocks):
            line+=frames[j][i]
        alpha_records.append(line)
    alpha_dataframe=pd.DataFrame(alpha_records)
    MOcoeffs=[float(i) for i in alpha_dataframe[MO-1].to_list()]
    if report==False:
        return MOcoeffs
    else:
        return MOcoeffs,alphaMOs
def read_basis(filename,basis,unique_species,debug=0):
    fp=open(filename,mode='r')
    cp2k_basis_sets=fp.readlines()
    fp.close

    basis_dict={}

    for element in unique_species:
        basis_dict[element]={}
        if debug==1:
            print(f'*** Species: {element}\n*** Basis set breakdown:\n')
        counter=0
        for j,i in enumerate(cp2k_basis_sets):
            line=i.strip().split()
            if line[0]==element and line[2]==basis:
                current=j
                nset=int(cp2k_basis_sets[j+1])
                current+=1
                for k in range(nset):
                    # read set header
                    current+=1
                    if debug==1:
                        print(cp2k_basis_sets[current],end='')
                    lmin,lmax,nexp=[int(cp2k_basis_sets[current].split()[i]) for i in [1,2,3]]
                    nshell=[int(cp2k_basis_sets[current].split()[i+4]) for i in range(lmax-lmin+1)]
                    set_matrix=[]
                    for l in range(nexp):
                        current+=1
                        set_matrix.append([float(i) for i in cp2k_basis_sets[current].split()])
                    if debug==1:
                        for l in set_matrix:
                            for m in l:
                                print(f'{m}\t',end='')
                            print()
                        print(f'lmin={lmin}')
                        print(f'lmax={lmax}')
                    c_counter=0
                    for l in range(lmin,lmax+1):
                        if debug==1:
                            print(f'l={l} has nshell={nshell[l-lmin]} contractions')
                        for n in range(nshell[l-lmin]):
                            c_counter+=1
                            for m in range(-l,l+1):
                                counter+=1
                                if debug==1:
                                    print(f'gto_{element}_{counter} = r^{l} * (',end='')
                                coeff=[]
                                alpha=[]
                                for kk in range(nexp):
                                    if abs(set_matrix[kk][c_counter])>1.0e-9:
                                        if debug==1:
                                            print(f'{set_matrix[kk][c_counter]} * exp(-{set_matrix[kk][0]}*r^2) + ',end='')
                                    coeff.append(set_matrix[kk][c_counter])
                                    alpha.append(set_matrix[kk][0])
                                if debug==1:
                                    print(f') * Y({l},{m})')
                                basis_dict[element][counter]={'l':l,'m':m,'coeff':coeff,'alpha':alpha}
        if debug==1:
            print('\n--------------------------------------------------------------------------\n')
        # GTOs: spherical
        basis_dict[element]['GTOs']=counter
        # CGTOs: cartesian - up to l=2 supported right now...
        basis_dict[element]['CGTOs']=basis_dict[element]['GTOs']
        for i in range(counter):
            if basis_dict[element][i+1]['m']==-2:
                basis_dict[element]['CGTOs']+=1
    return basis_dict
def read_CP2K_GTOs(species,basis_dict):
    atoms=len(species)
    GTOs=sum([basis_dict[i]['GTOs'] for i in species])
    CGTOs=sum([basis_dict[i]['CGTOs'] for i in species])
    GTO_depth=[]
    pcoeff=[]
    palpha=[]
    pqn=[]
    bfnPerAtom=[]
    for i in species:
        bfnPerAtom.append(basis_dict[i]['CGTOs'])
        for j in range(basis_dict[i]['GTOs']):
            if basis_dict[i][j+1]['m']==0:
                # s type
                if basis_dict[i][j+1]['l']==0:
                    GTO_depth.append(len(basis_dict[i][j+1]['alpha']))
                    for l in basis_dict[i][j+1]['coeff']:
                        pcoeff.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        palpha.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,0,0])
                # p type
                if basis_dict[i][j+1]['l']==1:
                    for k in range(3):
                        GTO_depth.append(len(basis_dict[i][j+1]['alpha']))
                        for l in basis_dict[i][j+1]['coeff']:
                            pcoeff.append(l)
                        for l in basis_dict[i][j+1]['alpha']:
                            palpha.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([1,0,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,1,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,0,1])
                # d type
                if basis_dict[i][j+1]['l']==2:
                    for k in range(6):
                        GTO_depth.append(len(basis_dict[i][j+1]['alpha']))
                        for l in basis_dict[i][j+1]['coeff']:
                            pcoeff.append(l)
                        for l in basis_dict[i][j+1]['alpha']:
                            palpha.append(l)
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([2,0,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([1,1,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([1,0,1])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,2,0])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,1,1])
                    for l in basis_dict[i][j+1]['alpha']:
                        pqn.append([0,0,2])
    return pcoeff,palpha,pqn,bfnPerAtom,GTO_depth
def initialize_STOs(species,x,y,z,STO_dict,debug=0):
    one_species=[]
    two_species=[]
    three_species=[]
    for key,value in STO_dict.items():
        for i in value:
            if i.find('1')!=-1:
                one_species.append(key)
                break
        for i in value:
            if i.find('2')!=-1:
                two_species.append(key)
                break
        for i in value:
            if i.find('3')!=-1:
                three_species.append(key)
                break
    atoms=len(species)
    STO_type_string={1:'1s',2:'2s',3:'2px',4:'2py',5:'2pz',6:'3s',7:'3px',8:'3py',9:'3pz'}
    STOs=sum([STO_dict[i]['STOs'] for i in species])
    STO_id_array=[]
    STO_type_array=[]
    STO_mu_array=[]
    j=-1
    for i in range(atoms):
        if species[i] in one_species:
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(1)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[1][0:2]])
        if species[i] in two_species:
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(2)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[2][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(3)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[3][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(4)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[4][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(5)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[5][0:2]])
        if species[i] in three_species:
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(6)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[6][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(7)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[7][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(8)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[8][0:2]])
            j+=1
            STO_id_array.append(i+1)
            STO_type_array.append(9)
            STO_mu_array.append(STO_dict[species[i]][STO_type_string[9][0:2]])
    if debug==1:
        print(f'Number of atoms: {atoms}')
        print('[Atomic id], species and coordinates (Ang|Bohr):')
        for i in range(atoms):
            print(f'[{i+1}] {species[i]} {x[i]/1.8897259886:>10.6f} {y[i]/1.8897259886:>10.6f} {z[i]/1.8897259886:>10.6f} | {x[i]:>10.6f} {y[i]:>10.6f} {z[i]:>10.6f}')
        print('--------------------------------------------------------------------------')
        print(f'Number of total STO basis functions: {STOs}\nAtomic decomposition:\n[id]\tSTO\torbital\tmu\ttype')
        for j,i in enumerate(STO_id_array):
            print(f'[{i}]\t{j+1}\t{STO_type_string[STO_type_array[j]]}\t{STO_mu_array[j]}\t{STO_type_array[j]}')
    return STOs,STO_id_array,STO_type_array,STO_mu_array
def calculate_overlap_S_matrix(x,y,z,STOs,STO_id_array,STO_type_array,STO_mu_array):
    Smatrix=np.identity(STOs)

    for i in range(STOs):
        for j in range(STOs):
            if STO_id_array[i]!=STO_id_array[j]:
                Smatrix[i][j]=overlap(x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],STO_mu_array[i],STO_mu_array[j],STO_type_array[i],STO_type_array[j])
    
    
    return Smatrix
def STO_GTO_projection(x,y,z,STO_id_array,STO_type_array,STO_mu_array,pcoeff,palpha,pqn,bfnPerAtom,GTO_depth,MOcoeffs,Smatrix):
    projection_dict={
        1:{
            'oohatac':[0.0154199, 0.203216, 0.401002, 0.313305, 0.144647, 0.049609, 0.0139344, 0.00318156, 0.000622761, 0.0000886225],
            'oohataa':[0.0372356, 0.0848076, 0.191073, 0.458596, 1.18834, 3.34954, 10.514, 37.9276, 156.411, 1188.35],
            'Lx':0,
            'Ly':0,
            'Lz':0,
        },
        2:{
            'oohatac':[0.00578359, 0.225661, 0.527807, 0.265263, 0.0457229, -0.0488578, -0.0223909, -0.00600006, -0.00120626, -0.000163319],
            'oohataa':[0.0214081, 0.0497851, 0.0991734, 0.2010660, 0.263514, 1.3073, 3.98351, 14.1623, 60.6176, 432.099],
            'Lx':0,
            'Ly':0,
            'Lz':0,
        },
        3:{
            'oohatac':[0.162114, 0.490328, 0.362480, 0.118961, 0.0237179, 0.00306629],
            'oohataa':[0.0554971, 0.127302, 0.31, 0.834323, 2.56, 10.3],
            'Lx':1,
            'Ly':0,
            'Lz':0,
        },
        4:{
            'oohatac':[0.162114, 0.490328, 0.362480, 0.118961, 0.0237179, 0.00306629],
            'oohataa':[0.0554971, 0.127302, 0.31, 0.834323, 2.56, 10.3],
            'Lx':0,
            'Ly':1,
            'Lz':0,
        },
        5:{
            'oohatac':[0.162114, 0.490328, 0.362480, 0.118961, 0.0237179, 0.00306629],
            'oohataa':[0.0554971, 0.127302, 0.31, 0.834323, 2.56, 10.3],
            'Lx':0,
            'Ly':0,
            'Lz':1,
        },
        6:{
            'oohatac':[0.00621804, 0.434003, 0.678264, 0.00507343, -0.163049, -0.049145, -0.00574937, -0.000273434],
            'oohataa':[0.014179, 0.0408735, 0.0803231, 0.201715, 0.376225, 0.995934, 3.2689, 15.8],
            'Lx':0,
            'Ly':0,
            'Lz':0,
        },
        7:{
            'oohatac':[0.0196080, 0.324519, 0.521925, 0.207248, 0.00504818, -0.0139732, -0.00394290, -0.000523609],
            'oohataa':[0.0247242, 0.0510310, 0.100108, 0.207164, 0.424917, 1.24038, 3.81650, 15.8360],
            'Lx':1,
            'Ly':0,
            'Lz':0,
        },
        8:{
            'oohatac':[0.0196080, 0.324519, 0.521925, 0.207248, 0.00504818, -0.0139732, -0.00394290, -0.000523609],
            'oohataa':[0.0247242, 0.0510310, 0.100108, 0.207164, 0.424917, 1.24038, 3.81650, 15.8360],
            'Lx':0,
            'Ly':1,
            'Lz':0,
        },
        9:{
            'oohatac':[0.0196080, 0.324519, 0.521925, 0.207248, 0.00504818, -0.0139732, -0.00394290, -0.000523609],
            'oohataa':[0.0247242, 0.0510310, 0.100108, 0.207164, 0.424917, 1.24038, 3.81650, 15.8360],
            'Lx':0,
            'Ly':0,
            'Lz':1,
        },
    }
    V_array=[]
    for isto in range(len(STO_id_array)):
        S=0
        mm=-1
        for ll in range(len(projection_dict[STO_type_array[isto]]['oohatac'])):
            mm+=1
            k=-1
            m=-1
            for i in range(len(x)):
                for j in range(bfnPerAtom[i]):
                    k+=1
                    for l in range(GTO_depth[k]):
                        m+=1
                        integral=anIres(
                                x[STO_id_array[isto]-1],
                                y[STO_id_array[isto]-1],
                                z[STO_id_array[isto]-1],
                                STO_mu_array[isto]**2*projection_dict[STO_type_array[isto]]['oohataa'][mm],
                                projection_dict[STO_type_array[isto]]['Lx'],
                                projection_dict[STO_type_array[isto]]['Ly'],
                                projection_dict[STO_type_array[isto]]['Lz'],
                                x[i],
                                y[i],
                                z[i],
                                palpha[m],
                                pqn[m][0],
                                pqn[m][1],
                                pqn[m][2])
                        S+=projection_dict[STO_type_array[isto]]['oohatac'][mm]*MOcoeffs[k]*pcoeff[m]*integral
        V_array.append(S)
    STOs=len(STO_id_array)
    res=np.linalg.solve(Smatrix,V_array).round(decimals=12)
    orb_compl=sum([res[i]*res[j]*Smatrix[i][j] for i in range(STOs) for j in range(STOs)])
    res/=math.sqrt(abs(orb_compl))
    STO_matrix=[[0,0,0,0] for i in range(len(x))]
    for i in range(STOs):
        if STO_type_array[i]==1:
            STO_matrix[STO_id_array[i]-1][0]=res[i]
        else:
            STO_matrix[STO_id_array[i]-1][(STO_type_array[i]-2)%4]=res[i]
    return STO_matrix,orb_compl,V_array
def create_cube_file(species,x,y,z,STO_matrix,STO_dict,filename='test.cube',cube_grid=0.5,offset=5.0,print_thres=1.0e-20):
    def atomic_contrib(x,y,z,X,Y,Z,species,line,STO_matrix,STO_dict):
        c1s,c2s,c2px,c2py,c2pz,c3s,c3px,c3py,c3pz=[0 for i in range(9)]
        mu1s,mu2s,mu2p,mu3s,mu3p=[0 for i in range(5)]
        if '1s' in STO_dict[species]:
            c1s=STO_matrix[line][0]
            mu1s=STO_dict[species]['1s']
        if '2s' in STO_dict[species]:
            c2s,c2px,c2py,c2pz=[i for i in STO_matrix[line]]
            mu2s=STO_dict[species]['2s']
            mu2p=STO_dict[species]['2p']
        if '3s' in STO_dict[species]:
            c3s,c3px,c3py,c3pz=STO_matrix[line]
            mu3s=STO_dict[species]['3s']
            mu3p=STO_dict[species]['3p']
        R=math.sqrt((x-X)**2+(y-Y)**2+(z-Z)**2)
        res=(0.5641895835477563*c1s*mu1s**1.5)/math.exp(mu1s*R) \
        + (0.32573500793527993*c2s*mu2s**2.5*R)/math.exp(mu2s*R) \
        + (0.11894160774351806*c3s*mu3s**3.5*R*R)/math.exp(mu3s*R) \
        + (0.5641895835477563*mu2p**2.5*(c2px*(x - X) + c2py*(y - Y) + c2pz*(z - Z)))/math.exp(mu2p*R) \
        + (0.20601290774570113*mu3p**3.5*R*(c3px*(x - X) + c3py*(y - Y) + c3pz*(z - Z)))/math.exp(mu3p*R)
        return res
    atoms=len(species)
    xmin,ymin,zmin=[min(i)-offset for i in [x,y,z]]
    xmax,ymax,zmax=[max(i)+offset for i in [x,y,z]]
    resx,resy,resz=[int((xmax-xmin)/cube_grid+1),int((ymax-ymin)/cube_grid+1),int((zmax-zmin)/cube_grid+1)]
    with open(filename,mode='w') as fp:
        print('My CUBE\n',file=fp)
        print(f'{atoms}\t{xmin}\t{ymin}\t{zmin}',file=fp)
        print(f'{resx}\t{cube_grid}\t{0.0}\t{0.0}',file=fp)
        print(f'{resy}\t{0.0}\t{cube_grid}\t{0.0}',file=fp)
        print(f'{resz}\t{0.0}\t{0.0}\t{cube_grid}',file=fp)
        species_dict={'H':1,'C':6,'N':7,'O':8,'F':9,'S':16}
        for j,i in enumerate(species):
            print(f'{species_dict[i]}\t{float(species_dict[i]):>.2}\t{x[j]}\t{y[j]}\t{z[j]}',file=fp)
        for ix in range(resx):
            for iy in range(resy):
                for iz in range(resz):
                    X,Y,Z=[xmin+ix*cube_grid,ymin+iy*cube_grid,zmin+iz*cube_grid]
                    mysum=0
                    for i in range(atoms):
                        mysum+=atomic_contrib(X,Y,Z,x[i],y[i],z[i],species[i],i,STO_matrix,STO_dict)
                    if mysum**2>print_thres:
                        if mysum<0.0:
                            print(f'{mysum**2:>.2} ',end='',file=fp)
                        else:
                            print(f'{-mysum**2:>.2} ',end='',file=fp)
                    else:
                        print('0.0 ',end='',file=fp)
                    if iz % 6 == 5:
                        print(file=fp)
                print(file=fp)

def resolve_STO_matrix(atoms,STOs,STO_id_array,STO_type_array,px,py,pz,AOM_array):
    STO_matrix=[[0,0,0,0] for i in range(atoms)]
    for i in range(STOs):
        if STO_type_array[i]==1:
            STO_matrix[STO_id_array[i]-1][0]=0
        else:
            STO_matrix[STO_id_array[i]-1][0]=0
            STO_matrix[STO_id_array[i]-1][1]=px[STO_id_array[i]-1]*AOM_array[STO_id_array[i]-1]
            STO_matrix[STO_id_array[i]-1][2]=py[STO_id_array[i]-1]*AOM_array[STO_id_array[i]-1]
            STO_matrix[STO_id_array[i]-1][3]=pz[STO_id_array[i]-1]*AOM_array[STO_id_array[i]-1]
    return STO_matrix
def create_cube_file(species,x,y,z,STO_matrix,STO_dict,filename='test.cube',cube_grid=0.5,offset=5.0,print_thres=1.0e-20):
    def atomic_contrib(x,y,z,X,Y,Z,species,line,STO_matrix,STO_dict):
        c1s,c2s,c2px,c2py,c2pz,c3s,c3px,c3py,c3pz=[0 for i in range(9)]
        mu1s,mu2s,mu2p,mu3s,mu3p=[0 for i in range(5)]
        if '1s' in STO_dict[species]:
            c1s=STO_matrix[line][0]
            mu1s=STO_dict[species]['1s']
        if '2s' in STO_dict[species]:
            c2s,c2px,c2py,c2pz=[i for i in STO_matrix[line]]
            mu2s=STO_dict[species]['2s']
            mu2p=STO_dict[species]['2p']
        if '3s' in STO_dict[species]:
            c3s,c3px,c3py,c3pz=STO_matrix[line]
            mu3s=STO_dict[species]['3s']
            mu3p=STO_dict[species]['3p']
        R=math.sqrt((x-X)**2+(y-Y)**2+(z-Z)**2)
        res=(0.5641895835477563*c1s*mu1s**1.5)/math.exp(mu1s*R) \
        + (0.32573500793527993*c2s*mu2s**2.5*R)/math.exp(mu2s*R) \
        + (0.11894160774351806*c3s*mu3s**3.5*R*R)/math.exp(mu3s*R) \
        + (0.5641895835477563*mu2p**2.5*(c2px*(x - X) + c2py*(y - Y) + c2pz*(z - Z)))/math.exp(mu2p*R) \
        + (0.20601290774570113*mu3p**3.5*R*(c3px*(x - X) + c3py*(y - Y) + c3pz*(z - Z)))/math.exp(mu3p*R)
        return res
    atoms=len(species)
    xmin,ymin,zmin=[min(i)-offset for i in [x,y,z]]
    xmax,ymax,zmax=[max(i)+offset for i in [x,y,z]]
    resx,resy,resz=[int((xmax-xmin)/cube_grid+1),int((ymax-ymin)/cube_grid+1),int((zmax-zmin)/cube_grid+1)]
    with open(filename,mode='w') as fp:
        print('My CUBE\n',file=fp)
        print(f'{atoms}\t{xmin}\t{ymin}\t{zmin}',file=fp)
        print(f'{resx}\t{cube_grid}\t{0.0}\t{0.0}',file=fp)
        print(f'{resy}\t{0.0}\t{cube_grid}\t{0.0}',file=fp)
        print(f'{resz}\t{0.0}\t{0.0}\t{cube_grid}',file=fp)
        species_dict={'H':1,'C':6,'N':7,'O':8,'F':9,'S':16}
        for j,i in enumerate(species):
            print(f'{species_dict[i]}\t{float(species_dict[i]):>.2}\t{x[j]}\t{y[j]}\t{z[j]}',file=fp)
        for ix in range(resx):
            for iy in range(resy):
                for iz in range(resz):
                    X,Y,Z=[xmin+ix*cube_grid,ymin+iy*cube_grid,zmin+iz*cube_grid]
                    mysum=0
                    for i in range(atoms):
                        mysum+=atomic_contrib(X,Y,Z,x[i],y[i],z[i],species[i],i,STO_matrix,STO_dict)
                    if mysum**2>print_thres:
                        if mysum<0.0:
                            print(f'{mysum**2:>.2} ',end='',file=fp)
                        else:
                            print(f'{-mysum**2:>.2} ',end='',file=fp)
                    else:
                        print('0.0 ',end='',file=fp)
                    if iz % 6 == 5:
                        print(file=fp)
                print(file=fp)
def AOM_overlap_calculation(istart,istop,jstart,jstop,x,y,z,STO_id_array,STO_type_array,STO_mu_array,STO_matrix):
    S=0
    for i in range(istart,istop):
        locali=(STO_type_array[i]-2)%4
        for j in range(jstart,jstop):
            localj=(STO_type_array[j]-2)%4
            if locali>0 and localj>0:
                if STO_id_array[i]!=STO_id_array[j]:
                    S=S+STO_matrix[STO_id_array[i]-1][locali]*STO_matrix[STO_id_array[j]-1][localj]\
                        *overlap(x[STO_id_array[i]-1],y[STO_id_array[i]-1],z[STO_id_array[i]-1],
                                 x[STO_id_array[j]-1],y[STO_id_array[j]-1],z[STO_id_array[j]-1],
                                 STO_mu_array[i],STO_mu_array[j],
                                 STO_type_array[i],STO_type_array[j])
                else:
                    if STO_type_array[i]==STO_type_array[j]:
                        S=S+STO_matrix[STO_id_array[i]-1][locali]*STO_matrix[STO_id_array[j]-1][localj]
    return S
def Sab(dimer_xyz_file,frag1_AOM_file,frag2_AOM_file,frag1_MO,frag2_MO,AOM_dict):
    if isinstance(frag1_AOM_file,dict)==True:
        frag1_AOM_dict=frag1_AOM_file
    else:
        with open(frag1_AOM_file) as fp:
            data=fp.read()
            frag1_AOM_dict=ast.literal_eval(data)
    if isinstance(frag2_AOM_file,dict)==True:
        frag2_AOM_dict=frag2_AOM_file
    else:
        with open(frag2_AOM_file) as fp:
            data=fp.read()
            frag2_AOM_dict=ast.literal_eval(data)
    frag1atoms=len(frag1_AOM_dict[frag1_MO])
    fp=open(dimer_xyz_file)
    data=fp.readlines()
    fp.close()
    x,y,z=[[float(i.split()[j]) for i in data[2::]] for j in [1,2,3]]
    species=[i.split()[0] for i in data[2::]]
    atoms=len(species)
    frag2atoms=atoms-frag1atoms
    frag1=single_molecule(None,'variable',[species[0:frag1atoms],x[0:frag1atoms],y[0:frag1atoms],z[0:frag1atoms]])
    frag2=single_molecule(None,'variable',[species[frag1atoms::],x[frag1atoms::],y[frag1atoms::],z[frag1atoms::]])
    frag1.initialize_STOs(AOM_dict)
    frag2.initialize_STOs(AOM_dict)
    frag1.resolve_pvecs()
    frag2.resolve_pvecs()
    STO_matrix_f1=resolve_STO_matrix(frag1atoms,frag1.STOs,frag1.STO_id_array,frag1.STO_type_array,frag1.px,frag1.py,frag1.pz,frag1_AOM_dict[frag1_MO])
    STO_matrix_f2=resolve_STO_matrix(frag2atoms,frag2.STOs,frag2.STO_id_array,frag2.STO_type_array,frag2.px,frag2.py,frag2.pz,frag2_AOM_dict[frag2_MO])
    S_f1=AOM_overlap_calculation(0,frag1.STOs,
                    0,frag1.STOs,
                    frag1.x+frag2.x,
                    frag1.y+frag2.y,
                    frag1.z+frag2.z,
                    frag1.STO_id_array+[i+frag1atoms for i in frag2.STO_id_array],
                    frag1.STO_type_array+frag2.STO_type_array,
                    frag1.STO_mu_array+frag2.STO_mu_array,
                    STO_matrix_f1+STO_matrix_f2)   
    for ci,i in enumerate(STO_matrix_f1):
        for j in range(1,3+1):
            STO_matrix_f1[ci][j]/=math.sqrt(abs(S_f1))
#     S_f1_norm=AOM_overlap_calculation(0,frag1.STOs,
#                     0,frag1.STOs,
#                     frag1.x+frag2.x,
#                     frag1.y+frag2.y,
#                     frag1.z+frag2.z,
#                     frag1.STO_id_array+[i+frag1atoms for i in frag2.STO_id_array],
#                     frag1.STO_type_array+frag2.STO_type_array,
#                     frag1.STO_mu_array+frag2.STO_mu_array,
#                     STO_matrix_f1+STO_matrix_f2)      
    S_f2=AOM_overlap_calculation(frag1.STOs,frag1.STOs+frag2.STOs,
                    frag1.STOs,frag1.STOs+frag2.STOs,
                    frag1.x+frag2.x,
                    frag1.y+frag2.y,
                    frag1.z+frag2.z,
                    frag1.STO_id_array+[i+frag1atoms for i in frag2.STO_id_array],
                    frag1.STO_type_array+frag2.STO_type_array,
                    frag1.STO_mu_array+frag2.STO_mu_array,
                    STO_matrix_f1+STO_matrix_f2)   
    for ci,i in enumerate(STO_matrix_f2):
        for j in range(1,3+1):
            STO_matrix_f2[ci][j]/=math.sqrt(abs(S_f2))
#     S_f2_norm=AOM_overlap_calculation(frag1.STOs,frag1.STOs+frag2.STOs,
#                     frag1.STOs,frag1.STOs+frag2.STOs,
#                     frag1.x+frag2.x,
#                     frag1.y+frag2.y,
#                     frag1.z+frag2.z,
#                     frag1.STO_id_array+[i+frag1atoms for i in frag2.STO_id_array],
#                     frag1.STO_type_array+frag2.STO_type_array,
#                     frag1.STO_mu_array+frag2.STO_mu_array,
#                     STO_matrix_f1+STO_matrix_f2)  
    Sab=AOM_overlap_calculation(0,frag1.STOs,
                    frag1.STOs,frag1.STOs+frag2.STOs,
                    frag1.x+frag2.x,
                    frag1.y+frag2.y,
                    frag1.z+frag2.z,
                    frag1.STO_id_array+[i+frag1atoms for i in frag2.STO_id_array],
                    frag1.STO_type_array+frag2.STO_type_array,
                    frag1.STO_mu_array+frag2.STO_mu_array,
                    STO_matrix_f1+STO_matrix_f2)
    return Sab

def projection_reg_test(ref_data_init,STO_proj_dict,rtol=1.0e-3,atol=1.0e-6,target=None):
    if target is None:
        ref_data=ref_data_init
    else:
        ref_data={}
        for key,value in ref_data_init.items():
            if key in target:
                ref_data[key]=value
    test_data={i:{} for i in ref_data.keys()}
    total=len(test_data.keys())
    passed=0
    total_time=0
    for counter,(key,value) in enumerate(ref_data.items()):
        tar = tarfile.open(value['cp2k_output_archive'])
        cp2k_out=value["cp2k_output_archive"].split(".tgz")[0]+'.out'
        tar.extractall()
        tar.close()
        tic=time.perf_counter()
        mymol=single_molecule(value['xyz'])
        mymol.get_cp2k_info(value['MO'],cp2k_out,'../cp2k_files/GTH_BASIS_SETS','DZVP-GTH')
        mymol.initialize_STOs(STO_proj_dict)
        mymol.project()
        mymol.save_state(key,to_file=False)
        toc=time.perf_counter()
        os.system(f'rm {cp2k_out}')
        test_data[key]['test_time']=toc-tic
        total_time+=test_data[key]['test_time']
        test_data[key]['test']=mymol.state
        check=[]
        check.append(np.allclose(test_data[key]['test']['pvecs']['px'],ref_data[key]['reference']['pvecs']['px'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data[key]['test']['pvecs']['py'],ref_data[key]['reference']['pvecs']['py'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data[key]['test']['pvecs']['pz'],ref_data[key]['reference']['pvecs']['pz'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data[key]['test']['S_matrix'],ref_data[key]['reference']['S_matrix'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data[key]['test']['V_array'],ref_data[key]['reference']['V_array'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data[key]['test']['STO_matrix'],ref_data[key]['reference']['STO_matrix'],rtol=rtol, atol=atol))
        check.append(np.allclose(test_data[key]['test']['singular_values'],ref_data[key]['reference']['singular_values'],rtol=rtol, atol=atol))
        check.append(np.allclose(list(test_data[key]['test']['AOM_dict'].values()),list(ref_data[key]['reference']['AOM_dict'].values()),rtol=rtol, atol=atol))
        check.append(np.allclose(list(test_data[key]['test']['compl_dict'].values()),list(ref_data[key]['reference']['compl_dict'].values()),rtol=rtol, atol=atol))
        print(f'[{counter+1}/{total}] ',end='')
        if check==[True for i in check]:
            print(f'PASS\t{key}\t{test_data[key]["test_time"]:.2f}s')
            test_data[key]['test_status']='PASS'
            passed+=1
        else:
            print(f'! FAIL\t{key}\t{test_data[key]["test_time"]:.2f}s')
            test_data[key]['test_status']='FAIL'
    print(f'Total tests: {total}; passed: {passed}/{total}; failed {total-passed}/{total}')
    print(f'Execution time: {total_time:.2f}s')
    return test_data