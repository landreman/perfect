from perfect_simulation import perfect_simulation,normalized_perfect_simulation #to get info about simu1;3409;0clation
from generate_compatible_profiles import generate_compatible_profiles
from generate_psi_of_r import drdpsiN_of_r
import numpy
import shutil #copy to copy files

def simul_manipulation(simul,group,field,change):
    #increase field in group by change
    simul.inputs.inputfile[group][field]=simul.inputs.inputfile[group][field]+change    
    simul.inputs.changevar(group,field,simul.inputs.inputfile[group][field])

def change_species(simul,Z,mHat,index):
    #increase impurity Z and mass (roughly)
    simul.inputs.charges[index]=Z
    simul.inputs.masses[index]=mHat
    #str_of_array=' '.join(map(str,l))
    
    simul.inputs.changevar("speciesParameters","charges",' '.join(map(str,simul.inputs.charges)))
    simul.inputs.changevar("speciesParameters","masses",' '.join(map(str,simul.inputs.masses)))


#simulation initialization
input_filename="input.namelist"
norm_filename="norms.namelist"
species_filename="species"
output_filename="perfectOutput.h5"

simul=normalized_perfect_simulation(input_filename,norm_filename,species_filename,None)
species=simul.species_list
#indicate which index in the above list that corresponds to main ions, impurity ions, and electrons.
#main ions: will have density given by input. Potential calculated from this.
#impurity: will have density given by its eta and the potential
#electron: will be assumed to not violate perfect orderings and thus use inputs.
mI=0

Tped=0.9
a=0.7 #minor radius in meters
xwidth=0.03 #pedestal width in r
imp_conc=0.01 #factor relating impurity and main ion concentration at pedestal top
dxdpsiN=drdpsiN_of_r(a,simul) #generate dr/dpsi_N. Returns a function
upShift_denom=2


TScale_i=1.0
Tped_i=0.9
dTCoredx_i=-0.1*Tped_i/xwidth
dTpeddx_i=-0.1*Tped_i/xwidth
dTSOLdx_i=-0.1*Tped_i/xwidth

nScale_i=1.0
nped_i=0.4
dnCoredx_i=-0.6
dnpeddx_i=-nped_i/xwidth
#print dnpeddx_d
dnSOLdx_i=-0.6

profile_params={"Tped_"+species[mI]:Tped_i,"dTCoredx_"+species[mI]:dTCoredx_i,"dTpeddx_"+species[mI]:dTpeddx_i,"dTSOLdx_"+species[mI]:dTSOLdx_i,"TScale_"+species[mI]:TScale_i,
                "nped_"+species[mI]:nped_i,"dnCoredx_"+species[mI]:dnCoredx_i,"dnpeddx_"+species[mI]:dnpeddx_i,"dnSOLdx_"+species[mI]:dnSOLdx_i,"nScale_"+species[mI]:nScale_i,
                }
#
Z=2
He_conc=numpy.array(["0.01","1","100"])
c0 = 0.01
cf = 0.01  # CHANGE THIS IN SCAN
profile_params["nScale_"+species[mI]]=(Z*c0+1)/(Z*cf+1)
imp_conc=cf

generate_compatible_profiles(simul,mI=mI,imp_conc=imp_conc,xwidth=xwidth,a=a,dxdpsiN=dxdpsiN,upShift_denom=upShift_denom,sameflux=True,m2tanh=True,mode="periodic",**profile_params)
