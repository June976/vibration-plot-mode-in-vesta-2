import sys
import os
import pandas as pd
import phonopy
from phonopy.units import THzToCm
from phonopy.interface.calculator import write_crystal_structure
from pymatgen.io.vasp.inputs import Poscar
import numpy as np


def make_result_file_dir(dir_out_file):
    try:
        os.mkdir(dir_out_file)
    except FileExistsError:
        print('Directory already exist!')
    

def load_Phonon_Calculation(supercell_matrix=[1, 1, 1], 
                          unitcell_filename='./CONTCAR',
                          force_sets_filename='./FORCE_SETS'):
    phonon_run = phonopy.load(supercell_matrix=supercell_matrix, 
                              unitcell_filename=unitcell_filename,
                              force_sets_filename=force_sets_filename)
    return phonon_run


def write_Primitive_Poscar(phonon_run, dir_out_file, filename='primitive.vasp'):
    primitive_cell = phonon_run.get_primitive()
    write_crystal_structure(filename='{}/{}'.format(dir_out_file, filename), interface_mode='vasp', cell=primitive_cell)


def get_atom_cartesian_xyz(dir_out_file, primitive_poscar_filename='primitive.vasp'):
    pri_poscar = pri_poscar = Poscar.from_file('{}/{}'.format(dir_out_file, primitive_poscar_filename))
    pri_struct = pri_poscar.structure
    pri_cart_coords = pri_struct.cart_coords
    nat = pri_cart_coords.shape[0]
    atom_xyz = []
    for i in range(3*pri_cart_coords.shape[0]):
        atom_xyz.append(pri_cart_coords)
    return nat, atom_xyz


def write_Vesta_file_from_Poscar(cmd_vesta_path, dir_out_file, primitive_poscar_filename='primitive.vasp', vesta_filename='primitive.vesta'):
    pri_poscar = Poscar.from_file('{}/{}'.format(dir_out_file, primitive_poscar_filename))
    pri_struct = pri_poscar.structure
    cif_filename = '{}/struct.cif'.format(dir_out_file)
    pri_struct.to(fmt='cif', filename=cif_filename)
    os.system('{} -nogui -i {} -o {}'.format(cmd_vesta_path, cif_filename, '{}/{}'.format(dir_out_file, vesta_filename)))


def open_Vesta_file(dir_out_file, vesta_filename='primitive.vesta'):
    try:
        vesta  = open('{}/{}'.format(dir_out_file, vesta_filename),'r')
    except:
        print("Cannot find poscar.vesta in current directory")
        print("Usage:\n\tpython modes_to_vesta.py <vesta-filename.vesta>")
        sys.exit(0)
    return vesta


def get_Vesta_Front_End(vesta):
    vfile = vesta.read()
    vesta_front = vfile.split("VECTR")[0]
    vesta_end   = vfile.split("VECTT\n 0 0 0 0 0")[1]
    vesta.close()
    return vesta_front, vesta_end


def parse_Modes(phonon_run):
    frqs, eigvs = phonon_run.get_frequencies_with_eigenvectors(q=[0.0, 0.0, 0.0])
    eigvals = frqs*THzToCm
    eigvecs = []
    norms = []
    atom_norms = []
    for i in range(eigvs.shape[1]):
        band_i_eigvs = eigvs[:, i]
        one_band_eigvecs = []
        one_band_norms = []
        for j in range(int(eigvs.shape[0]/3)):
            one_band_one_atom_eigvec = np.real(band_i_eigvs[j*3:(j+1)*3])*(-1.0)
            one_band_eigvecs.append(one_band_one_atom_eigvec)
            one_band_norms.append(np.linalg.norm(one_band_one_atom_eigvec))
        eigvecs.append(one_band_eigvecs)
        atom_norms.append(one_band_norms)
        norms.append(sum(one_band_norms))
    mat_eigvecs = np.array(eigvecs)
    return eigvals, mat_eigvecs, norms, atom_norms


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
    supercell_matrix = [4, 4, 1]                                # The supercell matrix used for force sets calculation
    unitcell_filename = './CONTCAR'                             # The unit cell used for force sets calculation
    force_sets_filename = './FORCE_SETS'                        # The force sets file path calculated by Phonopy
    cmd_vesta_path = 'D:\VESTA-win64\VESTA-win64\VESTA.exe'     # The full path of VESTA command
    ############################################  initial parameters  ####################################################

    make_result_file_dir(dir_out_file)
    phonon_run = load_Phonon_Calculation(supercell_matrix, unitcell_filename, force_sets_filename)
    write_Primitive_Poscar(phonon_run, dir_out_file)
    nat, atom_cart = get_atom_cartesian_xyz(dir_out_file)
    write_Vesta_file_from_Poscar(cmd_vesta_path, dir_out_file)
    vesta = open_Vesta_file(dir_out_file)
    vesta_front, vesta_end = get_Vesta_Front_End(vesta)
    eigvals, eigvecs, norms, atom_norms = parse_Modes(phonon_run)
    write_Vesta_Mode(eigvecs, vesta_front, vesta_end, nat, scaling_factor, dir_out_file)
    write_Summary_Info_csv('phons_info.csv', eigvals, eigvecs, atom_norms, atom_cart, norms, dir_out_file)
