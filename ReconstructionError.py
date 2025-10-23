import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
import os


analysis_type = 'N170'
ori_betas_mat = loadmat(analysis_type + '_data.mat')
learned_betas_mat = loadmat(analysis_type + '_learned_betas.mat')


ori_betas = ori_betas_mat['data']
learned_betas = learned_betas_mat['reconstructed']
ori_betas = np.array(ori_betas)
learned_betas = np.array(learned_betas)


# shape [nBeta, nSubject, nChan, nTime]
ori_betas = np.transpose(ori_betas, (0, 3, 2, 1))
mae = np.abs((ori_betas - learned_betas) ** 2)
chan_time_mae = mae.mean(axis=(0,1))
nBeta, nSubj, nChan, nTime = mae.shape
print("MAE mean:", np.mean(mae))

yticks = np.arange(1, nChan+1)
time = np.arange(nTime)


# heat-map of all beta and subject
plt.figure(figsize=(16,13))
plt.imshow(chan_time_mae, aspect='auto', origin='lower', cmap='Reds')
plt.colorbar(label='MAE')
plt.yticks(np.arange(nChan), yticks)
plt.xlabel('Time')
plt.ylabel('Channel')
plt.title(analysis_type+'Channel * Time MAE')
plt.tight_layout()
plt.savefig(analysis_type+'_Channel_Time MAE.png')
#plt.show()
plt.close()

save_dir1 = os.path.join(analysis_type, 'Each subject MAE')
os.makedirs(save_dir1, exist_ok=True)

# heat-map per subject
for iSubj in range(nSubj):
    chan_time_mae = mae[:, iSubj].mean(axis=0)
    plt.figure(figsize=(16,13))
    plt.imshow(chan_time_mae, aspect='auto', origin='lower', cmap='Reds')
    plt.colorbar(label='MAE')
    plt.yticks(np.arange(nChan), yticks)
    plt.xlabel('Time')
    plt.ylabel('Channel')
    plt.title(f'{analysis_type}_Subject_{iSubj}_Channel * Time MAE')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir1, f'Subject_{iSubj}_Channel_Time_MAE.png'))
    #plt.show()
    plt.close()

save_dir2 = os.path.join(analysis_type, 'Each beta MAE')
os.makedirs(save_dir2, exist_ok=True)
# heat-map per beta
for iBeta in range(nBeta):
    chan_time_mae = mae[iBeta].mean(axis=0)
    plt.figure(figsize=(16,13))
    plt.imshow(chan_time_mae, aspect='auto', origin='lower', cmap='Reds')
    plt.colorbar(label='MAE')
    plt.yticks(np.arange(nChan), yticks)
    plt.xlabel('Time')
    plt.ylabel('Channel')
    plt.title(f'{analysis_type}_Beta_{iBeta}_Channel * Time MAE')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir2, f'Beta_{iBeta}_Channel_Time_MAE.png'))
    #plt.show()
    
print("ori_betas mean:", np.mean(ori_betas))
print("learned_betas mean:", np.mean(learned_betas))
save_dir3 = os.path.join(analysis_type, 'Reconstruction vs original beta time series')
os.makedirs(save_dir3, exist_ok=True)
# reconstruction and original per beta
for iBeta in range(nBeta):
    ori_mean = ori_betas[iBeta].mean(axis=(0, 1))
    recon_mean = learned_betas[iBeta].mean(axis=(0, 1))
    plt.figure(figsize=(16,13))
    plt.plot(time, ori_mean, label='original', color='blue')
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.title(f'{analysis_type} | Beta {iBeta+1}')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir3, f'Beta{iBeta+1}_original_TimeSeries.png'))
    plt.close()
    
    plt.figure(figsize=(16,13))
    plt.plot(time, recon_mean, label='reconstructed', color='red')
    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.title(f'{analysis_type} | Beta {iBeta+1}')
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir3, f'Beta{iBeta+1}_reconstructed_TimeSeries.png'))
    plt.close()
    

for iBeta in range(nBeta):
    save_subdir = os.path.join(save_dir3, f'Beta_{iBeta+1}_per_channel')
    os.makedirs(save_subdir, exist_ok=True)
    for iChan in range(nChan):
        ori_mean = ori_betas[iBeta,:,iChan,:].mean(axis=0)
        recon_mean = learned_betas[iBeta,:,iChan,:].mean(axis=0)
        plt.figure(figsize=(16,13))
        plt.plot(time, ori_mean, label='original', color='blue')
        plt.xlabel('Time')
        plt.ylabel('Amplitude')
        plt.title(f'Beta {iBeta+1} | Channel {iChan+1}')
        plt.tight_layout()
        plt.savefig(os.path.join(save_subdir, f'Beta{iBeta+1}_Channel{iChan+1}_original_TimeSeries.png'))
        plt.close()
        
        plt.figure(figsize=(16,13))
        plt.plot(time, recon_mean, label='reconstructed', color='red')
        plt.xlabel('Time')
        plt.ylabel('Amplitude')
        plt.title(f'Beta {iBeta+1} | Channel {iChan+1}')
        plt.tight_layout()
        plt.savefig(os.path.join(save_subdir, f'Beta{iBeta+1}_Channel{iChan+1}_reconstructed_TimeSeries.png'))
        plt.close()
