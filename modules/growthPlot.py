import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter
import colorsys

plt.rcParams["figure.figsize"] = (5, 2.5)  # default width and height for all plots
plt.rcParams["figure.dpi"] = 300           # resolution
plt.rcParams.update({'font.size': 7})

class growthPlot():
    def __init__(self, ylabel="Cell-type Frequency", legend=False):
        self.legend=legend
        self.fig,self.ax=plt.subplots(figsize=(9/2.54,6/2.54))
        self.ylabel=ylabel

        self.group_cl = {"Homoplasmic ancestral": "#2d64db",
                    "Homoplasmic novel": "#e58da8",
                    "Heteroplasmic": '#872d8eff',
                    "Heteroplasmic (monomers)":'#872d8eff',
                    "Fused plasmid carrier": '#872d8eff'}
        self.group_ls = {"Homoplasmic ancestral": '-',
                    "Homoplasmic novel": '-',
                    "Heteroplasmic": '-',
                    "Heteroplasmic (monomers)":'dotted',
                    "Fused plasmid carrier": 'dashed'}
        self.group_lw = {"Homoplasmic ancestral": 3,
                    "Homoplasmic novel": 3,
                    "Heteroplasmic": 3,
                    "Heteroplasmic (monomers)":2,
                    "Fused plasmid carrier": 2}
        
        
    def format(self):
        ax=self.ax
        plt.rcParams["font.family"] = "Helvetica"
        plt.rcParams["font.size"] = 8
        plt.xlabel("Number of transfer")
        plt.ylabel(self.ylabel)
        if self.legend: plt.legend()
        plt.yscale('log')
        
        # # x-axis
        x=plt.xlim()
        step = np.log2(100)
        xs = np.arange(0, max(x)+step, step)
        plt.xticks(xs)  # optional, only if you want tick labels
        xlabels=labels = [f"{i}" for i,val in enumerate(xs)]
        ax.set_xticklabels(xlabels)
        
        # y-axis
        plt.ylim((1e-5,.3))
        # Major ticks at 10^n
        ax.yaxis.set_major_locator(LogLocator(base=10))
        # Minor ticks at 2–9 between powers of 10
        ax.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(1, 10)))

        plt.subplots_adjust(left=0.15,
                            right=0.95,
                            bottom=0.2,
                            top=0.9
                        )
        plt.grid(which='both', linestyle='-', color='gray', alpha=0.2, linewidth=.5, zorder=1)

def paler(hex_color, vivid=1., brighter=1., hue=0 ):

    # Remove '#' if present
    hex_color = hex_color.lstrip('#')
    
    # Convert HEX → RGB (0–255)
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    
    # Normalize RGB to 0–1
    r, g, b = r/255.0, g/255.0, b/255.0
    
    # Convert RGB → HSV
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    
    # Adjust saturation
    s = max(0, min(s * vivid, 1))
    v = max(0, min(v * brighter, 1))
    h = h + hue;
    h-=int(h)
    
    # Convert HSV → RGB
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    
    # Convert back to HEX
    r, g, b = int(r*255), int(g*255), int(b*255)
    
    return f'#{r:02x}{g:02x}{b:02x}'