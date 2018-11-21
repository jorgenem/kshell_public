import numpy as np 
import matplotlib.pyplot as plt 
import sys
sys.path.insert(0, '../bin/')
import shellmodelutilities as smutil
from matplotlib.colors import LogNorm

# === Set matplotlib choices ===
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('text', usetex=True)
rc('errorbar', capsize=1.5) # Set error bar style




# Set bin width and range
bin_width = 0.20
Emax = 30
Ex_low = 10
Ex_high = 15
Nbins = int(np.ceil(Emax/bin_width))
Emax_adjusted = bin_width*Nbins
bins = np.linspace(0,Emax_adjusted,Nbins+1)
bins_middle = (bins[0:-1]+bins[1:])/2


# Set name to be used for saving figures
save_name = "Ne20"


# Define list of calculation input files and corresponding label names
inputfiles = [
  "summary_Ne20_usda.txt"

]
names = [
  r"$\mathrm{Ne20}$",
]


# Set a spin window by defining list of allowed initial [spins, parities]. 
Jpi_list = [
              # All calculated spins, even A:
              [0,+1],[2,+1],[4,+1],[6,+1],[8,+1],[10,+1],[12,+1],[14,+1],[16,+1],[18,+1],[20,+1],
               [0,-1],[2,-1],[4,-1],[6,-1],[8,-1],[10,-1],[12,-1],[14,-1],[16,-1],[18,-1],[20,-1],
              # All calculated spins, odd A:
              # [1,+1],[3,+1],[5,+1],[7,+1],[9,+1],[11,+1],[13,+1],[15,+1],[17,+1],[19,+1],[21,+1],
              #  [1,-1],[3,-1],[5,-1],[7,-1],[9,-1],[11,-1],[13,-1],[15,-1],[17,-1],[19,-1],[21,-1],
              # A single spin:
              # [[2,-1]]
              ]

# Initialize figure objects
f_rho, ax_rho =   plt.subplots(1,1)
f_gsf, ax_gsf =   plt.subplots(1,1)




print("Entering main loop to plot KSHELL calculations")
# Extract statistical quantities from calculation files
for radiation_type in ["M1"]:#,"E1"]: # Loop over M1 and E1 if you have calculated both
  for i in range(len(inputfiles)):
    inputfile = inputfiles[i]
    name = names[i]
    print("Name =", name)

    # Read energy levels
    levels = smutil.read_energy_levels(inputfile)
    Egs = levels[0,0] # Read out the absolute ground state energy, so we can get relative energies later

    # Read transition strengths
    transitions = smutil.read_transition_strengths(inputfile, type=radiation_type)

    rho = smutil.total_level_density(levels, bin_width, Emax)
    gsf = smutil.strength_function_average(levels, transitions, Jpi_list, bin_width, Ex_low, Ex_high, type=radiation_type)


    lw = 1.6 + 0.03*i
    dashes = (2+0.5*i, 0.1+0.1*i)
    linestyle = '-'
    # -- rho:
    if radiation_type == "M1": # Only do this once, not for both radiation types:
      ax_rho.step(bins, np.append(0,rho), where='pre', label=name, linewidth=lw, linestyle=linestyle, dashes=dashes)
    # -- gSF:
    ax_gsf.plot(bins_middle[0:len(gsf)], gsf, label=name+r" ${:s}\,\,$".format(radiation_type), linewidth=lw, linestyle=linestyle, dashes=dashes)






# General plot settings

# -- rho
ax_rho.set_yscale('log')
ax_rho.set_xlabel(r'$E_x \, \mathrm{(MeV)}$')
ax_rho.set_ylabel(r'$\rho \, \mathrm{(MeV^{-1})}$')
# ax_rho.set_xlim([-0.5,10])
# ax_rho.set_ylim([1e-1, 1e5])
ax_rho.legend(loc='best')

# -- gsf
ax_gsf.set_yscale('log')
ax_gsf.set_ylabel(r'$f\, \mathrm{(MeV^{-3})}$')
ax_gsf.set_xlabel(r'$E_\gamma\,\mathrm{(MeV)}$')
ax_gsf.set_ylim([1e-10,5e-8])
ax_gsf.set_xlim([0,Ex_high])
ax_gsf.legend(loc='best')



f_rho.savefig(save_name+"_rho.pdf")
f_gsf.savefig(save_name+"_gsf.pdf")

plt.show()
