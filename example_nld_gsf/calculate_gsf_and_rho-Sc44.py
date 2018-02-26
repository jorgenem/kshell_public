from __future__ import division
import numpy as np 
import matplotlib.pyplot as plt 
import sys
sys.path.insert(0, '/home/jorgenem/gitrepos/kshell/bin/')
import shellmodelutilities as smutil
from matplotlib.colors import LogNorm

# === Set matplotlib choices ===
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rc('errorbar', capsize=1.5) # Set error bar style




# Set bin width and range
bin_width = 0.20
Emax = 30
Ex_low = 0
Ex_high = 30
Nbins = int(np.ceil(Emax/bin_width))
Emax_adjusted = bin_width*Nbins
bins = np.linspace(0,Emax_adjusted,Nbins+1)
bins_middle = (bins[0:-1]+bins[1:])/2


# Set name to be used for saving figures
save_name = "Sc44"


# Define list of calculation input files and corresponding label names
inputfiles = [
  "summary_Sc44_sdpf-mu.txt"

]
names = [
  r"$\mathrm{Sc44}$",
]


# Set a spin window by defining list of allowed initial [spins, parities]. 
Jpi_lists = [
              # All spins:
              [
              [0,+1],[2,+1],[4,+1],[6,+1],[8,+1],[10,+1],[12,+1],[14,+1],[16,+1],[18,+1],[20,+1],[22,+1],[24,+1],[26,+1],[28,+1],
               [0,-1],[2,-1],[4,-1],[6,-1],[8,-1],[10,-1],[12,-1],[14,-1],[16,-1],[18,-1],[20,-1],[22,-1],[24,-1],[26,-1],[28,-1],
              ],
              # A single spin:
              # [[2,+1]]
              ]

# Initialize figure objects
f_rho, ax_rho =   plt.subplots(1,1)
f_gsf, ax_gsf =   plt.subplots(1,1)




print "Entering main loop to plot KSHELL calculations"
# Extract statistical quantities from calculation files
for radiation_type in ["M1","E1"]: # Loop over E1, M1
  for i in range(len(inputfiles)):
    inputfile = inputfiles[i]
    name = names[i]
    Jpi_list = Jpi_lists[0]
    print "Name =", name

    # Read energy levels
    levels = smutil.read_energy_levels(inputfile)
    Egs = levels[0,0] # Read out the absolute ground state energy, so we can get relative energies later

    # Read transition strengths
    transitions = smutil.read_transition_strengths(inputfile, type=radiation_type)

    rho = smutil.total_level_density(levels, bin_width, Emax)
    gsf = smutil.strength_function_average_updated_definition(levels, transitions, Jpi_list, bin_width, Ex_low, Ex_high, type=radiation_type)


    lw = 1.6 + 0.03*i
    dashes = (2+0.5*i, 0.1+0.1*i)
    linestyle = '-'
    # -- rho:
    if radiation_type == "M1": # Only do this once, not for every iteration
      ax_rho.step(bins, np.append(0,rho), where='pre', label=name, linewidth=lw, linestyle=linestyle, dashes=dashes)
    # -- gSF:
    ax_gsf.plot(bins_middle[0:len(gsf)], gsf, label=name+r"${:s}\,\,$".format(radiation_type), linewidth=lw, linestyle=linestyle, dashes=dashes)






# General plot settings

# -- rho
ax_rho.set_yscale('log')
ax_rho.set_xlabel(r'$E_x \, \mathrm{(MeV)}$', fontsize=22)
ax_rho.set_ylabel(r'$\rho \, \mathrm{(MeV^{-1})}$', fontsize=22)
ax_rho.set_xlim([-0.5,10])
ax_rho.set_ylim([1e-1, 1e5])
ax_rho.legend(loc='best', fontsize=17)
ax_rho.tick_params(labelsize=16)
ax_rho.grid(True)
f_rho.set_size_inches(8,6)
f_rho.subplots_adjust(left=0.12, right=0.95, top=0.97, bottom=0.12)
f_rho.savefig(save_name+"_rho.pdf")

# -- gsf
ax_gsf.set_yscale('log')
ax_gsf.set_ylabel(r'$f\, \mathrm{(MeV^{-3})}$', fontsize=22)
ax_gsf.set_xlabel(r'$E_\gamma\,\mathrm{(MeV)}$', fontsize=22)
ax_gsf.set_ylim([1e-13,5e-7])
ax_gsf.set_xlim([0,Ex_high])
ax_gsf.tick_params(labelsize=16)
ax_gsf.legend(loc='best', fontsize=15)
ax_gsf.grid(True)
f_gsf.subplots_adjust(left=0.08, right=0.95, top=0.97, bottom=0.12)
f_gsf.subplots_adjust(wspace=0)
f_gsf.set_size_inches(8,6)
f_gsf.savefig(save_name+"_gsf.pdf")



plt.show()
