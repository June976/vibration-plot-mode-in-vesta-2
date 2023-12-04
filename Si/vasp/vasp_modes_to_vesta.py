import sys
import os
import re
from pymatgen.io.vasp.inputs import Poscar
import pandas as pd
import numpy as np


def make_result_file_dir(dir_out_file):
    try:
        os.mkdir(dir_out_file)
    except FileExistsError:
        print('Directory already exist!')


def parse_nat_from_poscar(poscar_filename):
    pri_poscar = Poscar.from_file(poscar_filename)
    pri_struct = pri_poscar.structure
    nat = len(pri_struct.species)
    return nat


def write_vesta_file_from_poscar(cmd_vesta_path, poscar_filename, dir_out_file, vesta_filename='poscar.vesta'):
    pri_poscar = Poscar.from_file(poscar_filename)
    pri_struct = pri_poscar.structure
    cif_filename = '{}/struct.cif'.format(dir_out_file)
    pri_struct.to(fmt='cif', filename=cif_filename)
    os.system('{} -nogui -i {} -o {}/{}'.format(cmd_vesta_path, cif_filename, dir_out_file, vesta_filename))


def parse_Modes(outcar, nat):
    eigvals = [ 0.0 for i in range(nat*3) ]
    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms   = [ 0.0 for i in range(nat*3) ]
    atom_cart = [ 0.0 for i in range(nat*3) ]
    atom_norms = [ 0.0 for i in range(nat*3) ]
    outcar.seek(0) # just in case
    while True:
        line = outcar.readline()
        if not line:
            break
        if "Eigenvectors and eigenvalues of the dynamical matrix" in line:
            outcar.readline() # ----------------------------------------------------
            outcar.readline() # empty line
            for i in range(nat*3):
                outcar.readline() # empty line
                p = re.search(r'^\s*(\d+).+?([\.\d]+) cm-1', outcar.readline())
                eigvals[i] = float(p.group(2))
                outcar.readline() # X         Y         Z           dx          dy          dz
                eigvec = []
                cart_xyz = []
                atom_norm = []
                for j in range(nat):
                    tmp = outcar.readline().split()
                    eigvec.append([ float(tmp[x]) for x in range(3,6) ])
                    cart_xyz.append([ float(tmp[x]) for x in range(3) ])
                    atom_norm.append(sum([ float(tmp[x])**2 for x in range(3,6) ]))       
                eigvecs[i] = eigvec
                atom_cart[i] = cart_xyz
                atom_norms[i] = atom_norm
                norms[i] = ( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )**0.5
        if "Eigenvectors and eigenvalues of the dynamical matrix" in line:
            break
    outcar.close()
    eigvals.reverse()
    eigvecs.reverse()
    norms.reverse()
    atom_cart.reverse()
    atom_norms.reverse()   
    mat_eigvecs = np.array(eigvecs)
    return eigvals, mat_eigvecs, norms, atom_cart, atom_norms


def write_Vesta_Mode(eigvecs, vesta_front, vesta_end, nat, scaling_factor, dir_out_file):
    for index in range(nat*3):
        with open("{}/mode_{}.vesta".format(dir_out_file, index+1), 'w') as modef:
            eigvec = eigvecs[index]
            modef.write(vesta_front)
            sf = scaling_factor
            towrite = "VECTR\n"
            for i in range(1,1+nat):
                towrite += "%4d%9.5f%9.5f%9.5f\n"%(i,eigvec[i-1][0]*sf,eigvec[i-1][1]*sf,eigvec[i-1][2]*sf)
                towrite += "%5d  0   0    0    0\n 0 0 0 0 0\n"%i
            towrite += " 0 0 0 0 0\n"
            towrite += "VECTT\n"
            for i in range(1,1+nat):
                towrite += "%4d%6.3f 255   0   0 1\n"%(i,0.5)
                towrite += " 0 0 0 0 0\n"
            if i==0:
                print(towrite)
            modef.write(towrite)
            modef.write(vesta_end)
    return 0


def open_Vesta_Outcar_file(vesta_filename, dir_out_file, outcar_filename='OUTCAR'):
    try:
        vesta  = open('{}/{}'.format(dir_out_file, vesta_filename),'r')
    except:
        print("Cannot find Vesta file in current directory")
        sys.exit(0)
    try:
        outcar = open(outcar_filename, 'r')
    except:
        print("Cannot find OUTCAR in current directory")
        sys.exit(0)
    return vesta, outcar


def get_Vesta_Front_End(vesta):
    vfile = vesta.read()
    vesta_front = vfile.split("VECTR")[0]
    vesta_end   = vfile.split("VECTT\n 0 0 0 0 0")[1]
    vesta.close()
    return vesta_front, vesta_end


def write_Summary_Info_csv(csv_filename, eigvals, eigvecs, atom_norms, atom_cart, norms, dir_out_file):
    frqs_count = len(eigvals)
    df_phon_vector_info = pd.DataFrame(data=eigvals, index=[i+1 for i in range(len(eigvals))], columns=['Freq(cm-1)'])
    for i in range(int(frqs_count/3)):
        df_phon_vector_info['atom_{}_cart'.format(i)] = [ atom_cart[x][i] for x in range(frqs_count) ]
        df_phon_vector_info['atom_{}_norms'.format(i)] = [ atom_norms[x][i] for x in range(frqs_count) ]
        df_phon_vector_info['atom_{}_vector'.format(i)] = [ x[i, :] for x in eigvecs ]
    df_phon_vector_info['total_norms'] = norms
    df_phon_vector_info.to_csv('{}/{}'.format(dir_out_file, csv_filename))


if __name__ == '__main__':
    ############################################  initial parameters  ###################################################
    dir_out_file = 'out'                                        # The directory for saving result file   
    scaling_factor = 3                                          # The scale factor of the length of vector in VESTA file     
    poscar_filename = './POSCAR'                                # The full path of POSCAR for VASP frequency calculation
    outcar_filename = './OUTCAR'                                # The full path of OUTCAR generated in VASP frequency calculation
    cmd_vesta_path = 'D:\VESTA-win64\VESTA-win64\VESTA.exe'     # The full path of VESTA command
    ############################################  initial parameters  ###################################################

    make_result_file_dir(dir_out_file)
    write_vesta_file_from_poscar(cmd_vesta_path='D:\VESTA-win64\VESTA-win64\VESTA.exe', poscar_filename='POSCAR', dir_out_file=dir_out_file)
    vesta, outcar = open_Vesta_Outcar_file(vesta_filename='poscar.vesta', dir_out_file=dir_out_file)
    vesta_front, vesta_end = get_Vesta_Front_End(vesta)
    nat = parse_nat_from_poscar(poscar_filename='./POSCAR')
    eigvals, eigvecs, norms, atom_cart, atom_norms = parse_Modes(outcar, nat)
    write_Vesta_Mode(eigvecs, vesta_front, vesta_end, nat, scaling_factor, dir_out_file)
    write_Summary_Info_csv('phons_info.csv', eigvals, eigvecs, atom_norms, atom_cart, norms, dir_out_file)