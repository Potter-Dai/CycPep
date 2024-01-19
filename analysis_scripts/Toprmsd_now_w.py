import numpy as np
from sys import argv
import mdtraj as md
from Pipeline import DPA

pdbids = ['APGfP','GPfAP','FPaFPa','GPPGPP','PGPGPG','VPaVPa','PPGPLG','PGLVIY','AGVPVW','WPISFVP', \
          'PLIFSPI','AIPFNSL','IFPYPIP','SFLPVNL','IIILPPTP','IPPFFVIML','PPFFLIILV','ALLLVLVLP', \
          'ILLLVLVLP','PPIFVLPPYI','AFFPPAFFPP','AFFPPFFVPP','VPPFFVPPFF'] 

npydata_dir = '/share/udata4/daibt/work/cycpep/0.allanalysis/cluster_analysis/npydata/'
Z = int(argv[1])
Top10_RMSD_dict = {}
for i, pdbid in enumerate(pdbids):
    density = np.load(f'./npydata/density.{pdbid}.Z{Z}.npy')
    err_density = np.load(f'./npydata/err_density.{pdbid}.Z{Z}.npy')
    k_hat = np.load(f'./npydata/k_hat.{pdbid}.Z{Z}.npy')
    distance = np.load(f'./npydata/distance.{pdbid}.Z{Z}.npy')
    indice = np.load(f'./npydata/indice.{pdbid}.Z{Z}.npy')
    density = list(density)
    err_density = list(err_density)
    k_hat = list(k_hat)

    dpc = DPA._DensityPeakAdvanced(density, err_density, k_hat, distance, indice, Z) 
    yhalo = dpc[1]
    index = [i for i, x in enumerate(yhalo) if x == -1]  
    density = [i for j, i in enumerate(density) if j not in index]
    yhalo = [i for i in yhalo if i != -1]
    unique1, counts1 = np.unique(yhalo, return_counts=True)

    Top10_cluster_size_nhalo = [counts1[i] for i in np.argsort(counts1)[-10:]]

    dict_nhalo = dict(zip(unique1, counts1))
    Top10_cluster_label_nhalo = []
    for i in range(10):
        for k, v in dict_nhalo.items():
            if v == Top10_cluster_size_nhalo[i]:
                Top10_cluster_label_nhalo.append(k)

    Index = Top10_cluster_label_nhalo.index
    Top10_cluster_label_nhalo = list(set(Top10_cluster_label_nhalo))
    Top10_cluster_label_nhalo.sort(key=Index)

    trajdir = f'/share/udata4/daibt/work/cycpep/{pdbid}/700K/RSFF2C'
    t = md.load(f'{trajdir}/analysis/{pdbid}_rsff2c_700K.nc', top=f'{trajdir}/{pdbid}_noionwat.hmr.prmtop')
    ref_t = md.load(f'/share/udata4/daibt/work/cycpep/nativepdb/{pdbid}_dry.pdb')
    BB = t.topology.select('backbone')
    Centers = dpc[-1]
    Top10_centers = [Centers[i] for i in Top10_cluster_label_nhalo]
    if pdbid not in ['FPaFPa', 'GPPGPP', 'PGPGPG', 'VPaVPa', 'AFFPPAFFPP', 'VPPFFVPPFF']:
        RMSD_list = md.rmsd(t, ref_t, frame=0, atom_indices=BB)
    elif pdbid in ['FPaFPa', 'GPPGPP', 'VPaVPa', 'AFFPPAFFPP', 'VPPFFVPPFF']:
        ref_t2 = md.load(f'/share/udata4/daibt/work/cycpep/nativepdb/new_{pdbid}_dry.pdb')
        RMSD_list1 = md.rmsd(t, ref_t, frame=0, atom_indices=BB)
        RMSD_list2 = md.rmsd(t, ref_t2, frame=0, atom_indices=BB)
        RMSD_list = np.minimum(RMSD_list1, RMSD_list2)
    else:
        ref_t2 = md.load(f'/share/udata4/daibt/work/cycpep/nativepdb/new_{pdbid}_dry.pdb')
        ref_t3 = md.load(f'/share/udata4/daibt/work/cycpep/nativepdb/new2_{pdbid}_dry.pdb')
        RMSD_list1 = md.rmsd(t, ref_t, frame=0, atom_indices=BB)
        RMSD_list2 = md.rmsd(t, ref_t2, frame=0, atom_indices=BB)
        RMSD_list3 = md.rmsd(t, ref_t3, frame=0, atom_indices=BB)
        RMSD_list = np.minimum(RMSD_list1, RMSD_list2, RMSD_list3)
    Top10_RMSD = RMSD_list[Top10_centers]*10
    # print('Top10 RMSD of {} is {}'.format(pdbid, ['{:.3f}'.format(i) for i in Top10_RMSD]))
    Top10_RMSD_dict[pdbid] = Top10_RMSD  
np.save(f'./npytoprmsd/Top10RMSD_noweight_rm-1_Z{Z}.npy', Top10_RMSD_dict)