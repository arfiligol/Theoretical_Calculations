import matplotlib.pyplot as plt
import numpy as np
from qutip import *

def hamiltonian(Ec, Ej, N, ng):

    m = np.diag(4 * Ec * (np.arange(-N,N+1)-ng)**2) + 0.5 * Ej * (np.diag(-np.ones(2*N), 1) + 
                                                               np.diag(-np.ones(2*N), -1))
    return Qobj(m)

def plot_energies(ng_vec, energies, ymax=(8, 0)):

    fig, ax0 = plt.subplots(nrows =1, figsize=(12,8),dpi =500)

    for n in range(len(energies[0,:])):
        ax0.plot(ng_vec, energies[:,n],lw=10)
    for n in range(len(energies1[0,:])):
        ax0.plot(ng_vec, energies1[:,n],'k--',lw=5)    
    ax0.set_xlim(min(ng_vec),max(ng_vec))
    ax0.set_ylim(-0.3, ymax[0])
    ax0.set_xlabel(r"$N_{g}$",size ='45')
    ax0.set_ylabel(r'$E_{n}$',fontsize=45)
    xscale=np.linspace(min(ng_vec),max(ng_vec),5) 
    for axis in ['top','bottom','left','right']:
        ax0.spines[axis].set_linewidth(7)
    ax0.set_xticks(xscale)
    ax0.set_xticklabels([str(xscale[i]) for i in range(len(xscale))])
    ax0.tick_params('both',labelsize='35', length=20, width=7)
    #ax0.set_title(r'$CPB\ energy\ levels $',size ='45',y=1.05)
    fig.tight_layout()

    
   
N = 5
Ec = 0.7
Ej = 1
ng_vec = np.linspace(-2, 2, 400)

energies = np.array([hamiltonian(Ec, Ej, N, ng).eigenenergies() for ng in ng_vec])
energies1 = np.array([hamiltonian(Ec, 0, N, ng).eigenenergies() for ng in ng_vec])
plot_energies(ng_vec, energies)


#%%
#Transmon
import matplotlib.pyplot as plt
import numpy as np
from qutip import *
def hamiltonian(Ec, Ej, N, ng):

    m = np.diag(4 * Ec * (np.arange(-N,N+1)-ng)**2) + 0.5 * Ej * (np.diag(-np.ones(2*N), 1) + 
                                                               np.diag(-np.ones(2*N), -1))
    return Qobj(m)

def plot_energies(ng_vec, energies, ymax=(10, 0)):

    fig, ax0 = plt.subplots(2,2, figsize=(12,6),dpi =200)
    for n in range(len(energies[0][0,:])):
        ax0[0][0].plot(ng_vec, energies[0][:,n],lw=5)
 
    ax0[0][0].set_xlim(min(ng_vec),max(ng_vec))
    ax0[0][0].set_ylim(-1, 8.8)
    ax0[0][0].set_ylabel(r'$E_{n}$',fontsize=20)
    xscale=np.linspace(min(ng_vec),max(ng_vec),9) 
    ax0[0][0].set_xticks(xscale)
    ax0[0][0].set_xticklabels([str(xscale[i]) for i in range(len(xscale))])
    ax0[0][0].tick_params(labelsize='15')
    ax0[0][0].set_title(r'$E_{J}/E_{C}=1$',size ='20')

    for n in range(len(energies[1][0,:])):
        ax0[0][1].plot(ng_vec, energies[1][:,n],lw=5)
        
    ax0[0][1].set_xlim(min(ng_vec),max(ng_vec))
    ax0[0][1].set_ylim(-1, 5)
    ax0[0][1].set_ylabel(r'$E_{n}$',fontsize=20)
    xscale=np.linspace(min(ng_vec),max(ng_vec),9) 
    ax0[0][1].set_xticks(xscale)
    ax0[0][1].set_xticklabels([str(xscale[i]) for i in range(len(xscale))])
    ax0[0][1].tick_params(labelsize='15')
    ax0[0][1].set_title(r'$E_{J}/E_{C}=5$',size ='20')

    for n in range(len(energies[2][0,:])):
        ax0[1][0].plot(ng_vec, energies[2][:,n],lw=5)
       
    ax0[1][0].set_xlim(min(ng_vec),max(ng_vec))
    ax0[1][0].set_ylim(-1, 3)
    ax0[1][0].set_xlabel(r"$N_{g}$",size ='20')
    ax0[1][0].set_ylabel(r'$E_{n}$',fontsize=20)
    xscale=np.linspace(min(ng_vec),max(ng_vec),9) 
    ax0[1][0].set_xticks(xscale)
    ax0[1][0].set_xticklabels([str(xscale[i]) for i in range(len(xscale))])
    ax0[1][0].tick_params(labelsize='15')
    ax0[1][0].set_title(r'$E_{J}/E_{C}=10$',size ='20')

    for n in range(len(energies[3][0,:])):
        ax0[1][1].plot(ng_vec, energies[3][:,n],lw=5)
      
    ax0[1][1].set_xlim(min(ng_vec),max(ng_vec))
    ax0[1][1].set_ylim(-1, 2)
    ax0[1][1].set_xlabel(r"$N_{g}$",size ='20')
    ax0[1][1].set_ylabel(r'$E_{n}$',fontsize=20)
    xscale=np.linspace(min(ng_vec),max(ng_vec),9) 
    ax0[1][1].set_xticks(xscale)
    ax0[1][1].set_xticklabels([str(xscale[i]) for i in range(len(xscale))])
    ax0[1][1].tick_params(labelsize='15')
    ax0[1][1].set_title(r'$E_{J}/E_{C}=50$',size ='20')
    fig.tight_layout()
   
N=10
ng_vec = np.linspace(-2, 2, 400)        
energies = [np.array([hamiltonian(1, 1, N, ng).eigenenergies() for ng in ng_vec]),
            np.array([hamiltonian(0.1, 0.5, N, ng).eigenenergies() for ng in ng_vec]),
            np.array([hamiltonian(0.1, 1, N, ng).eigenenergies() for ng in ng_vec]),
            np.array([hamiltonian(0.1, 5, N, ng).eigenenergies() for ng in ng_vec])]


plot_energies(ng_vec, energies)

