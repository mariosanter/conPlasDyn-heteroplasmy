from collections import defaultdict

from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import binom,hypergeom

from modules.growthPlot import growthPlot
from  modules.growthPlot import paler

class partitioningModel:
    def __init__(self, n, apar=0):
        """
        Initializes the partitioning model

        Args:                                                               (var in MS)
        n (int, positive): plasmid copy number                                      (d)
        apar (float, [0,1]): probability of active par of plasmids at cell div      (a)
        """

        if n % 2 == 1 : # uneven
            raise Exception("Uneven plasmid copy number")
        self.n = n

        self.apar=apar

        self.P = np.zeros((self.n+1,self.n+1))

        for j in range(n+1):
            dist = self.distribution(j)
            for i, p in dist.items():
                self.P[i,j] += p
        
        if not np.allclose(self.P.sum(axis=0), 1.0): 
            raise Exception("M columns do not sum up to 1")

    def distribution(self, j):
        """
        Computes the distribution of daughter cell types produced at cell division of given type

        Args:
        j (int, {0,...,n}) : cell type of dividing cell, no. ofnovel plasmids

        """

        apar=self.apar
        n=self.n
        # Initializes dict of state distribution (p(x)=0 for all x) 

        jn,jw = j,self.n-j
        Jwpar_dist = defaultdict(float)
        for ki, pi in zip(np.arange(0,jw+1),binom.pmf(np.arange(0,jw+1),jw, apar)):
            if pi>0: Jwpar_dist[ki] = pi
        Jnpar_dist = defaultdict(float)
        for ki, pi in zip(np.arange(0,jn+1),binom.pmf(np.arange(0,jn+1),jn, apar)):
            if pi>0: Jnpar_dist[ki] = pi
        
        dist = defaultdict(float)

        for jwpar,pjwpar in list(Jwpar_dist.items()):
            for jnpar,pjnpar in list(Jnpar_dist.items()):
                    jwfree = 2*jw - jwpar
                    jnfree = 2*jn - jnpar
                    ndraw = n - jwpar - jnpar
                    nsucc = jwfree
                    ntotal = jwfree + jnfree
                    for Iwfree in range(0,ndraw+1):
                        pIwfree = hypergeom.pmf(Iwfree,M=ntotal,n=nsucc,N=ndraw) *\
                                    pjwpar * pjnpar
                        I_w = jwpar + Iwfree
                        if pIwfree>0: 
                            I_m = n-I_w
                            dist[I_m] += pIwfree

        return dist
    
    def growth(self,gen,f):
        '''
        analogue to method dimerizationModel.dimerizationModel.growth
        '''

        x_t0 = np.zeros(self.n+1)
        x_t0[0]=1-f
        x_t0[1]=f

        x_t=x_t0
        x_t_list=[x_t]
        for gen_i in range(1,gen+1):
            x_t = self.P @ x_t
            x_t_list.append(x_t.copy())

        self.x_t_arr = np.array(x_t_list)

    def growth_groupPlot(self,legend=False,plot=None,alpha=1.0,zorder=None,vivid=1,brighter=1,hue=0):
        # Grouping function
        def group_state(i, n):
            if i == 0:
                return "Homoplasmic ancestral"
            elif i == n:
                return "Homoplasmic novel"
            else:
                return "Heteroplasmic"
        groups = {i: group_state(i, self.n) for i in range(0,self.n+1)}
        group_names = ["Homoplasmic ancestral",
                    "Homoplasmic novel",
                    "Heteroplasmic"]

        # Build group trajectories
        gen=len(self.x_t_arr)-1
        group_traj = {g: np.zeros(gen+1) for g in group_names}
        for i in range(self.n+1):
            g = groups[i]
            group_traj[g] += self.x_t_arr[:, i]

        if plot==None:
            self.plot=growthPlot(legend=legend)
        else:
            self.plot=plot

        plot=self.plot

        plot_groups = ["Heteroplasmic"]

        for g in plot_groups:
            plt.plot(range(gen+1), group_traj[g], label=g, 
                     lw=plot.group_lw["Heteroplasmic (monomers)"], 
                     c=paler(plot.group_cl[g],vivid=vivid,brighter=brighter,hue=
                             hue),
                     ls=plot.group_ls["Heteroplasmic (monomers)"],
                     alpha=1,
                     zorder=zorder
                     )
            plt.plot(range(gen+1), group_traj[g], label=g, 
                     lw=plot.group_lw[g], 
                     c=paler(plot.group_cl[g],vivid=vivid,brighter=brighter,hue=
                             hue),
                     ls=plot.group_ls[g],
                     alpha=0.5,
                     zorder=zorder
                     )
        plot.format()

        return self
    
    def conjugation(self):
        # Analogue to modules.dimerizationModel.dimerizationModel.conjugation

        # time series (array) of transconjugant (tc) types (t)
        # rows (len x_t_array): time 
        # columns (3): ancestral,novel 
        self.tc_t_arr = np.zeros((len(self.x_t_arr), 2))

        for t in range(0, len(self.x_t_arr)):
            for i in range(0,self.n+1):

                i_freq = self.x_t_arr[t, i]
                if i_freq ==0: continue

                # ancestral replicons
                self.tc_t_arr[t,0] += i_freq * (self.n - i) / self.n
                #novel replicons
                self.tc_t_arr[t,1] += i_freq * i / self.n

    def conjugation_plot(self):
        # Plots the transconjugant frequencies
        plot=growthPlot("Transconjugant cell-type frequency")
        gen=len(self.tc_t_arr)-1

        plt.plot(range(gen+1) , self.tc_t_arr[:,0],
                 label="Homoplasmic ancestral", 
                 c=plot.group_cl["Homoplasmic ancestral"])
        plt.plot(range(gen+1) , self.tc_t_arr[:,1],
                 label="Homoplasmic novel", 
                 c=plot.group_cl["Homoplasmic novel"])
        plot.format()
        return plot.fig
