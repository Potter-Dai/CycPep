import numpy as np
import pytraj as pt


class GetDistMatrix:
    def __init__(self, traj, top):
        self.traj = pt.load(traj, top)
        self.top =  self.traj.topology
        
    def get_phipsiomg_matrix(self):
        data = pt.multidihedral(self.traj, dihedral_types='phi psi omega')  #Datalist of phi, psi, omega
        data = data.to_ndarray()
        data = data.T
        return data
    
    def period_dihe(self, data_i, angle_data, sequence):
        N_resi = self.top.n_residues
        L = len(sequence)
        if L%2 == 0:
            pointer = L//2
            pointer2 = L//3         
            if sequence[:pointer] == sequence[pointer:]:  # If the first half of the sequence is the same as the second half
                new_data_i = np.concatenate((data_i[pointer:], data_i[:pointer]))  # Swap the first and second halves of the data_i
                data_i = np.abs(data_i - angle_data)
                data_i = np.minimum(data_i, 360 - data_i)
                data_i2 = np.abs(new_data_i - angle_data)
                data_i2 = np.minimum(data_i2, 360 - data_i2)
                data_i = np.sqrt(((data_i**2).sum(axis=1))/N_resi)
                data_i2 = np.sqrt(((data_i2**2).sum(axis=1))/N_resi)
                data_i = np.minimum(data_i, data_i2)
            
            elif sequence[:pointer2] == sequence[pointer2:pointer2*2] == sequence[pointer2*2:]: 
                new_data_i = np.concatenate((data_i[pointer2:pointer2*2], data_i[pointer2*2:], data_i[:pointer2]))
                new2_data_i = np.concatenate((data_i[pointer2*2:], data_i[:pointer2], data_i[pointer2:pointer2*2]))
                data_i = np.abs(data_i - angle_data)
                data_i = np.minimum(data_i, 360 - data_i)
                data_i2 = np.abs(new_data_i - angle_data)
                data_i2 = np.minimum(data_i2, 360 - data_i2)
                data_i3 = np.abs(new2_data_i - angle_data)
                data_i3 = np.minimum(data_i3, 360 - data_i3)
                data_i = np.sqrt(((data_i**2).sum(axis=1))/N_resi)
                data_i2 = np.sqrt(((data_i2**2).sum(axis=1))/N_resi)
                data_i3 = np.sqrt(((data_i3**2).sum(axis=1))/N_resi)
                data_i = np.minimum(data_i, data_i2, data_i3)
            
            else:
                data_i = np.abs(data_i - angle_data)
                data_i = np.minimum(data_i, 360 - data_i)
                data_i = np.sqrt(((data_i**2).sum(axis=1))/N_resi)
                data_i = data_i
        else:
            data_i = np.abs(data_i - angle_data)
            data_i = np.minimum(data_i, 360 - data_i)
            data_i = np.sqrt(((data_i**2).sum(axis=1))/N_resi)
            data_i = data_i
        return data_i
    
    def cal_distance_matrix(self, data, sequence):
        distance_matrix = np.zeros((len(data), len(data)))
        for i in range(len(data)):
            data_i = data[i, :]
            distance_matrix[:, i] = self.period_dihe(data_i, data, sequence)
        return distance_matrix  

'''
useage:
get_distance_matrix = GetDistMatrix('traj.nc', 'amberparm.prmtop')
'''