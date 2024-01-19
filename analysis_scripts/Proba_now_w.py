import numpy as np
from Pipeline import DPA
import matplotlib.pyplot as plt
from matplotlib import cm


def cal_weight(logps, sim_temp, tar_temp):
    logps = logps - np.min(logps)
    dlogps = logps * (sim_temp / tar_temp - 1)
    dps = 1 / np.exp(dlogps)
    return dps


def getProbReweight(clusterlabels,TopClusterlabel,weight):
    Top_index = []
    for i, label in enumerate(clusterlabels):
        if label == TopClusterlabel:
            Top_index.append(i)
        else:
            continue
    Top_reweight = np.empty(len(Top_index))
    for i, index in enumerate(Top_index):
        Top_reweight[i] = weight[index]
    probability_reweight = np.sum(Top_reweight) / np.sum(weight)
    return probability_reweight


def getProbNoReweight(clusterlabels,TopClusterlabel,noweight):
    Top_index = []
    for i, label in enumerate(clusterlabels):
        if label == TopClusterlabel:
            Top_index.append(i)
        else:
            continue
    Top_noweight = np.empty(len(Top_index))
    for i, index in enumerate(Top_index):
        Top_noweight[i] = noweight[index]
    probability_noweight = np.sum(Top_noweight) / np.sum(noweight)
    return probability_noweight


def getTop10Prob(clusterlabels, logps, Top10_label):
    temps = [300, 700]
    logps = - np.array(logps)
    weight = cal_weight(logps, temps[1], temps[0])
    noweight = np.ones(np.shape(logps)[0])

    ###calcluate the sum weight of Top cluster###
    Prob_weight = np.empty(10)
    Prob_noweight = np.empty(10)
    for i in range(10):
        Prob_weight[i] = getProbReweight(clusterlabels,Top10_label[i],weight)
        Prob_noweight[i] = getProbNoReweight(clusterlabels,Top10_label[i],noweight)
    return Prob_weight, Prob_noweight


pdbids = ['APGfP','GPfAP','FPaFPa','GPPGPP','PGPGPG','VPaVPa','PPGPLG','PGLVIY','AGVPVW','WPISFVP', \
          'PLIFSPI','AIPFNSL','IFPYPIP','SFLPVNL','IIILPPTP','IPPFFVIML','PPFFLIILV','ALLLVLVLP', \
          'ILLLVLVLP','PPIFVLPPYI','AFFPPAFFPP','AFFPPFFVPP','VPPFFVPPFF']

Z = int(input('please input the value of Z \n'))
# Z = float(input('please input the value of Z \n'))
fig, axs = plt.subplots(4, 6, figsize=(24, 16), sharex=True, sharey=True)
ax = axs.flatten() 
width = 1/10
x = np.arange(1)
colors = [cm.tab20c(i) for i in np.linspace(0, 0.6, 10)]
Prob_weight_dict = {}
Prob_noweight_dict = {}
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
    unique, counts = np.unique(yhalo, return_counts=True)
    Top10_cluster_size_nhalo = [counts[i] for i in np.argsort(counts)[-10:]]  # 从小到大排序

    dict_nhalo = dict(zip(unique, counts))
    Top10_cluster_label_nhalo = []
    for m in range(10):
        for k, v in dict_nhalo.items():
            if v == Top10_cluster_size_nhalo[m]:
                Top10_cluster_label_nhalo.append(k)
                
    Index = Top10_cluster_label_nhalo.index
    Top10_cluster_label_nhalo = list(set(Top10_cluster_label_nhalo))
    Top10_cluster_label_nhalo.sort(key=Index)

    Prob_weight, Prob_noweight = getTop10Prob(yhalo, density, Top10_cluster_label_nhalo)
    
    # Save the calculated Top probability
    Prob_weight_dict[pdbid] = Prob_weight
    Prob_noweight_dict[pdbid] = Prob_noweight

    for j in range(9, -1, -1):
        ax[i].bar(x - j*width-0.05, Prob_noweight[j], width, color=colors[9-j], ec=colors[9-j], fill=True, lw=1.5, label=f'C{10-j}_Prob now')
    for j in range(10):
        ax[i].bar(x + j*width+0.05, Prob_weight[9-j], width, ec=colors[j], hatch='/', fill=False, lw=1.5, label=f'C{j+1}_Prob w')

    ax[i].set_xticks([])
    ax[i].set_yticks(np.arange(0, 1.1, 0.2))
    ax[i].tick_params(labelsize=18)
    ax[i].set_title(f'{pdbid}', fontsize=18)
np.save(f'./npydata/Prob_weight_rm-1.Z{Z}.npy', Prob_weight_dict)
np.save(f'./npydata/Prob_noweight_rm-1.Z{Z}.npy', Prob_noweight_dict)

fig.text(0.5, 0.92, 'Probability of Top 10 Clusters', fontsize=20, ha='center')
fig.text(0.08, 0.5, 'Probability', va='center', rotation='vertical', fontsize=20)
ax[23].remove()
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, ncol=1, fontsize=14, loc='center right')
plt.savefig(f'Prob_noW_W_top10_halo_rm-1_Z{Z}.jpg', dpi=300)