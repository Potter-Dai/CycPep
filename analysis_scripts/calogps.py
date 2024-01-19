import sys
sys.path.append('/share/udata4/daibt/scripts')
import GetDistMatrix as gdm
from Pipeline import DPA
import numpy as np
import mdtraj as md


def getFinalDistanceMatrix(pdbid):
    trajdir = f'/share/udata4/daibt/work/cycpep/{pdbid}/700K/RSFF2C'
    t = md.load(f'{trajdir}/analysis/{pdbid}_rsff2c_700K.nc', top=f'{trajdir}/{pdbid}_noionwat.hmr.prmtop')
    N_resi = t.n_residues
    tp = t.topology
    h_C = tp.select('name C and resid 0')
    h_CA = tp.select('name CA and resid 0')
    h_N = tp.select('name N and resid 0')
    t_C = tp.select(f'name C and resid {N_resi-1}')
    t_CA = tp.select(f'name CA and resid {N_resi-1}')
    t_N = tp.select(f'name N and resid {N_resi-1}')
    datas = gdm.get_phipsiomg_matrix(t, h_C, h_CA, h_N, t_C, t_CA, t_N)

    if pdbid == 'FPaFPa' or pdbid == 'GPPGPP' or pdbid == 'VPaVPa':
        distance_matrix = gdm.get_distance_matrix(datas, N_resi)
        distance_2nd_matrix = gdm.get_2nd_distance_matrix_FPaFPa(datas)
        final_distance_matrix = gdm.get_final_distance_matrix_FPaFPa_AFFPP2(distance_matrix, distance_2nd_matrix)
    elif pdbid == 'PGPGPG':
        distance_matrix = gdm.get_distance_matrix(datas, N_resi)
        distance_2nd_matrix = gdm.get_2nd_distance_matrix_PGPGPG(datas)
        distance_3th_matrix = gdm.get_3th_distance_matrix_PGPGPG(datas)
        final_distance_matrix = gdm.get_final_distance_matrix_PGPGPG(distance_matrix, distance_2nd_matrix, distance_3th_matrix)   
    elif pdbid == 'AFFPPAFFPP' or pdbid == 'VPPFFVPPFF':
        distance_matrix = gdm.get_distance_matrix(datas, N_resi)
        distance_2nd_matrix = gdm.get_2nd_distance_matrix_AFFPPAFFPP(datas)
        final_distance_matrix = gdm.get_final_distance_matrix_FPaFPa_AFFPP2(distance_matrix, distance_2nd_matrix) 
    else:
        final_distance_matrix = gdm.get_distance_matrix(datas, N_resi)
    return final_distance_matrix


def caldensity(distance_matrix, Z):
    est = DPA.DensityPeakAdvanced(Z=Z, frac=0.3, metric='precomputed')
    est.fit(distance_matrix)
    density = est.densities_
    err_density = est.err_densities_
    k_hat = est.k_hat_
    distance = est.nn_distances_
    indice = est.nn_indices_
    return density, err_density, k_hat, distance, indice

pdbids = ['APGfP','GPfAP','FPaFPa','GPPGPP','PGPGPG','VPaVPa','PPGPLG','PGLVIY','AGVPVW','WPISFVP', \
          'PLIFSPI','AIPFNSL','IFPYPIP','SFLPVNL','IIILPPTP','IPPFFVIML','PPFFLIILV','ALLLVLVLP', \
          'ILLLVLVLP','PPIFVLPPYI','AFFPPAFFPP','AFFPPFFVPP','VPPFFVPPFF']
for pdbid in pdbids:
    final_distance_matrix = getFinalDistanceMatrix(pdbid)
    for i in [1.5, 2.5, 3.5, 4.5]:
        alldata = caldensity(final_distance_matrix, i)
        density = alldata[0]
        err_density = alldata[1]
        k_hat = alldata[2]
        distance = alldata[3]
        indice = alldata[4]
        np.save(f'./npydata/density.{pdbid}.Z{i}.npy', density)
        np.save(f'./npydata/err_density.{pdbid}.Z{i}.npy', err_density)
        np.save(f'./npydata/k_hat.{pdbid}.Z{i}.npy', k_hat)
        np.save(f'./npydata/distance.{pdbid}.Z{i}.npy', distance)
        np.save(f'./npydata/indice.{pdbid}.Z{i}.npy', indice)
