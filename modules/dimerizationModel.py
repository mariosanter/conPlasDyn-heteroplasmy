from collections import defaultdict
import itertools
from math import factorial
import math
from scipy.stats import binom
import numpy as np
import matplotlib.pyplot as plt
from modules.growthPlot import growthPlot
from modules.growthPlot import paler 
class dimerizationModel:
    def __init__(self, n, ddim=0, dmon=0):
        """
        Initializes the model of dimerization

        Args:                                                           (var in MS)
        n (int, even, positive): plasmid copy number
        ddim (float, [0,1]): probability of dimerization of monomers    (d)
        dmon (float, [0,1]): probability of monomeriziation of dimers   (d')
        """

        if n % 2 == 1 : # uneven
            raise Exception("Uneven plasmid copy number")
        self.n = n
        
        self.ddim=ddim
        self.dmon=dmon

        self.states=[
            (iwm,inm,iwd,ind,iwnd)
            for iwm,inm,iwd,ind,iwnd in itertools.product(range(n+1), repeat=5)
            if iwm + inm + 2*iwd + 2*ind + 2*iwnd  == n
        ]

        self.state_to_index = {s:i for i,s in enumerate(self.states)}
        self.index_to_state = {i:s for i,s in enumerate(self.states)}

        self.k = len(self.states)
        self.P = np.zeros((self.k,self.k))

        for s in self.states:
            j = self.state_to_index[s]
            dist = self.distribution(*s)  # returns dict of {state: prob}
            for s_next, p in dist.items():
                i = self.state_to_index[s_next]
                self.P[i,j] += p

        if not np.allclose(self.P.sum(axis=0), 1.0): 
            raise Exception("M columns do not sum up to 1")

    def distribution(self, iwm, inm, iwd, ind, iwnd):
        """
        Computes the distribution of daughter cell types produced at cell division of given type

        Args: 
        (iwm, inm, iwd, ind, iwnd) : cell type of dividing parental cell 
        iwm : no of wild-type plasmids in monomeric form 
        inm : no of novel-type plasmids in monomeric form 
        iwd : no of wild-type plasmid replicons in dimeric form 
        ind : no of novel-type plasmid replicons in dimeric form 
        ind : no of plasmid replicons in dimeric form with a heterogenic composition of 
                plasmids
        """


        dmon=self.dmon
        ddim=self.ddim

        n = iwm + inm + 2*iwd + 2*ind + 2*iwnd # plasmid copy number

        # Initializes dict of state distribution (p(x)=0 for all x) 
        all_states_dist = defaultdict(float)
        
        # State distribution after plasmid replication (all plasmid replicons replicated) and monomerization (plasmid dimers monomerize at a probability dmon)
        # Computes the joint probability distribution:
        # xwd, xnd, xwnd: no. of wt,novel,het dimers monomerizing

        for (xwd, xnd, xwnd) in itertools.product(
                                    range(0, 2*iwd + 1),
                                    range(0, 2*ind + 1),
                                    range(0, 2*iwnd + 1)):
            p = binom.pmf(xwd, 2*iwd, dmon) * \
                binom.pmf(xnd, 2*ind, dmon) * \
                binom.pmf(xwnd, 2*iwnd, dmon)
            if p>0:
                all_states_dist[(2*iwm + 2*xwd + xwnd, 
                                2*inm + 2*xnd + xwnd,
                                2*iwd - xwd, 
                                2*ind - xnd, 
                                2*iwnd - xwnd)] += p


        # Plasmid dimerization
        # Plasmid monomers form pairs and dimerize at probabilty ddim
        for (iwm, inm, iwd, ind, iwnd), p in list(all_states_dist.items()):
            all_states_dist.pop((iwm, inm, iwd, ind, iwnd), None)

            for (xw,xn) in itertools.product(
            # xw, xw: number of wild-type and novel-type plasmid pairs 
                                range(max(0,int((iwm-inm)/2)), int(iwm/2)+1),
                                range(max(0,int((inm-iwm)/2)), int(inm/2)+1),
                                            ):
                diff = (iwm - 2*xw) - (inm - 2*xn)
                if diff == 0:
                    # perfectly balanced -> heterodimers
                    xwn = iwm - 2*xw
                elif diff == 1 or diff == -1:
                    # one leftover monomer that cannot dimerize
                    xwn = min(iwm - 2*xw, inm - 2*xn)
                else:
                    continue  # too imbalanced, impossible

                # build pairs: calculate probability distribution of pair combinations
                combinationsSpec = \
                    math.comb(iwm,2*xw) * factorial(2*xw) / (2**xw) / factorial(xw)\
                * math.comb(inm,2*xn) * factorial(2*xn) / (2**xn) / factorial(xn)\
                * factorial(xwn)
                
                Nmon = iwm + inm
                combinationsAll = factorial(Nmon) / 2**(int(Nmon/2)) / factorial(int(Nmon/2))

                weight = combinationsSpec / combinationsAll
                
                for (xwp,xnp,xwnp) in itertools.product( # no of pairs ...
                                range(0,xw+1), # ... forming dimers is bin. distr.
                                range(0,xn+1),
                                range(0,xwn+1)
                                            ):
                    p1 = p * weight * \
                        binom.pmf(xwp, xw, ddim) *\
                        binom.pmf(xnp, xn, ddim) *\
                        binom.pmf(xwnp, xwn, ddim)
                    if p1>0:
                        all_states_dist[(iwm - 2 * xwp - xwnp,
                                        inm - 2 * xnp - xwnp,
                                        iwd + xwp,
                                        ind + xnp,
                                        iwnd + xwnp)] += p1

        # Random assortment: calculate probability distribution (m_j->i)
        for (iwm, inm, iwd, ind, iwnd), p in list(all_states_dist.items()):

            focal_state_dist=defaultdict(float)
            focal_state_dist[(iwm, inm, iwd, ind, iwnd)]=1 

            halfplasmidsdistributed=False
            while True:
                halfplasmidsdistributed=True
                for (jwm, jnm, jwd, jnd, jwnd), p1 in list(focal_state_dist.items()):
                    N = jwm + jnm + 2*jwd + 2*jnd + 2*jwnd
                    if N == n:
                        continue
                    if N-1 == n: # only monomers can be removed
                        R = jwm + jnm
                        if jwm>0: focal_state_dist[(jwm-1, jnm, jwd, jnd, jwnd)] += \
                            p1 * jwm / R
                        if jnm>0: focal_state_dist[(jwm, jnm-1, jwd, jnd, jwnd)] += \
                            p1 * jnm / R
                    else:
                        R = jwm + jnm + jwd + jnd + jwnd
                        if jwm>0: focal_state_dist[(jwm-1, jnm, jwd, jnd, jwnd)] += \
                            p1 * jwm / R
                        if jnm>0: focal_state_dist[(jwm, jnm-1, jwd, jnd, jwnd)] += \
                            p1 * jnm / R
                        if jwd>0: focal_state_dist[(jwm, jnm, jwd-1, jnd, jwnd)] += \
                            p1 * jwd / R
                        if jnd>0: focal_state_dist[(jwm, jnm, jwd, jnd-1, jwnd)] += \
                            p1 * jnd / R
                        if jwnd>0: focal_state_dist[(jwm, jnm, jwd, jnd, jwnd-1)] += \
                            p1 * jwnd / R
                    
                    focal_state_dist.pop((jwm, jnm, jwd, jnd, jwnd))

                    halfplasmidsdistributed=False
                
                if halfplasmidsdistributed:
                    break

            for key, value in focal_state_dist.items():
                all_states_dist[key] += p * value

            all_states_dist[(iwm, inm, iwd, ind, iwnd)]-=p
            if all_states_dist[(iwm, inm, iwd, ind, iwnd)]==0:
                all_states_dist.pop((iwm, inm, iwd, ind, iwnd), None)
            # 

        return all_states_dist
    
    # describe the allele dynamics of populations that grow infinitly
    def growth(self,gen,f):
        '''
        Args:
        Computes the relative frequencies of cell types for bacterial growth in discrete generations (SI, section 'mathematical description'/'Baseline model')
        gen (int): Number of generations (time) to be computed
        f (float, [0,1]): Initial frequency of cells carrying the novel allele in 1 copy
        '''

        x_t0 = np.zeros(self.k)
        x_t0[self.state_to_index[(self.n,0,0,0,0)]] = 1-f
        x_t0[self.state_to_index[(self.n-1,1,0,0,0)]] = f

        x_t=x_t0
        x_t_list=[x_t]
        for gen_i in range(1,gen+1):
            x_t = self.P @ x_t
            x_t_list.append(x_t.copy())

        self.x_t_arr = np.array(x_t_list)   # shape: (gen+1, k)

    def growth_groupPlot(self,legend=False,plot=None,alpha=1.0,zorder=None,vivid=1,brighter=1,hue=0):
        # Grouping function
        def group_state(s, n):
            iwm, inm, iwd, ind, iwnd = s
            if iwm + 2*iwd == n:
                return "Homoplasmic ancestral"
            elif inm + 2*ind == n:
                return "Homoplasmic novel"
            elif iwm > 0 and inm > 0 and iwnd == 0:
                return "Heteroplasmic (monomers)"
            else:
                return "Fused plasmid carrier"
            
        # Assign each state to a group
        groups = {s: group_state(s, self.n) for s in self.states}
        group_names = ["Homoplasmic ancestral",
                    "Homoplasmic novel",
                    "Heteroplasmic (monomers)",
                    "Fused plasmid carrier"]
        if plot==None:
            self.plot=growthPlot(legend=legend)
        else:
            self.plot=plot

        # Build group trajectories
        gen=len(self.x_t_arr)-1
        group_traj = {g: np.zeros(gen+1) for g in group_names}
        for i, s in enumerate(self.states):
            g = groups[s]
            group_traj[g] += self.x_t_arr[:, i]

        # plot cell-type groups 
        plot_groups = ["Heteroplasmic (monomers)", 
                    "Fused plasmid carrier"]
        # plot heteroplasmic
        plt.plot(range(gen+1), 
                group_traj["Heteroplasmic (monomers)"]+\
                    group_traj["Fused plasmid carrier"], 
                    label="Heteroplasmic" if legend else "",
                    lw=3, 
                    color=paler(self.plot.group_cl["Heteroplasmic"],vivid,brighter,hue), 
                    ls=self.plot.group_ls["Heteroplasmic"],
                    alpha=0.5,
                    zorder=zorder)
        for g in plot_groups:
            plt.plot(range(gen+1), group_traj[g], 
                     label=g if legend else "",
                     lw=self.plot.group_lw[g], 
                     c=paler(self.plot.group_cl[g],vivid,brighter,hue),
                     ls=self.plot.group_ls[g],
                     alpha=alpha,
                     zorder=zorder)
        
        # plt.title(f"Dimerization model n={self.n}, ddim(dmon) = {self.ddim} ({self.dmon})")
        self.plot.format()
        
        return self



    def conjugation(self,):
        # Transforms the genotype frequencies time series (self.x_t_arr,from growth simulation) into plasmid frequencies (to describe the transconjugant frequencies that resemble the relative frequency of plasmid types involved in conjugation)

        # time series (array) of transconjugant (tc) types (t)
        # rows (len x_t_array): time 
        # columns (3): wild-type, novel-type (monomeric or dimeric), heterozygous (dim.)
        self.tc_t_arr = np.zeros((len(self.x_t_arr), 3))

        for t in range(0, len(self.x_t_arr)):
            for s in self.states:
                i = self.state_to_index[s]
                
                state_freq=self.x_t_arr[t, i]
                if state_freq ==0: continue
                state_replicon_sum=sum(s)
                
                # Assume that all plasmid replicons have the same chance to be selected 
                # for conjugation
                
                # wild-type replicons
                self.tc_t_arr[t,0] += state_freq*(s[0] + s[2]) / state_replicon_sum 
                # novel-type replicons
                self.tc_t_arr[t,1] += state_freq*(s[1] + s[3]) / state_replicon_sum  
                # heterozygous replicons
                self.tc_t_arr[t,2] += state_freq * s[4] / state_replicon_sum


    def conjugation_plot(self):
        # Plots the transconjugant frequencies
        plot=growthPlot("Transconjugant cell-type frequency")
        gen=len(self.tc_t_arr)-1

        plt.plot(range(gen+1), self.tc_t_arr[:,0], 
                 label="Homoplasmic wild-type", 
                 c=plot.group_cl["Homoplasmic wild-type"],lw=2,)
        plt.plot(range(gen+1), self.tc_t_arr[:,1], 
                 label="Homoplasmic novel-type", 
                 c=plot.group_cl["Homoplasmic novel-type"],lw=2,)
        plt.plot(range(gen+1), self.tc_t_arr[:,2], 
                 label="Heterozygous",c=plot.group_cl["Heterozygous"],lw=2,)

        # plt.title(f"Dimerization model n={self.n}, ddim(dmon) = {self.ddim} ({self.dmon})")
        plot.format()
        return plot.fig
