import sys
sys.path.append('/share/udata4/daibt/scripts')
from GetDistMatrix import getdistmat
from Pipeline import DPA
import numpy as np


def getFinalDistanceMatrix(sequence):
    trajdir = f'/share/udata4/daibt/work/cycpep/{sequence}/700K/RSFF2C'
    gdm = getdistmat(f'{trajdir}/analysis/allreplica_gap100.nc', f'{trajdir}/{sequence}_noionwat.hmr.prmtop')
    data = gdm.get_phipsiomg_matrix()
    final_distance_matrix = gdm.cal_distance_matrix(data, sequence)
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

sequences = ['APGfP','GPfAP','FPaFPa','GPPGPP','PGPGPG','VPaVPa','PPGPLG','PGLVIY','AGVPVW','WPISFVP', \
          'PLIFSPI','AIPFNSL','IFPYPIP','SFLPVNL','IIILPPTP','IPPFFVIML','PPFFLIILV','ALLLVLVLP', \
          'ILLLVLVLP','PPIFVLPPYI','AFFPPAFFPP','AFFPPFFVPP','VPPFFVPPFF']
for sequence in sequences:
    final_distance_matrix = getFinalDistanceMatrix(sequence)
    for i in [1.5, 2.5, 3.5, 4.5]:
        alldata = caldensity(final_distance_matrix, i)
        density = alldata[0]
        err_density = alldata[1]
        k_hat = alldata[2]
        distance = alldata[3]
        indice = alldata[4]
        np.save(f'./npydata/density.{sequence}.Z{i}.npy', density)
        np.save(f'./npydata/err_density.{sequence}.Z{i}.npy', err_density)
        np.save(f'./npydata/k_hat.{sequence}.Z{i}.npy', k_hat)
        np.save(f'./npydata/distance.{sequence}.Z{i}.npy', distance)
        np.save(f'./npydata/indice.{sequence}.Z{i}.npy', indice)
