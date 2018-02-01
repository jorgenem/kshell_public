from __future__ import division
import numpy as np 
# import matplotlib.pyplot as plt 
import sys

# Scripts to sort levels and transition strengths from KSHELL 
# into Ex-Eg matrix energy bins, and to make level density and 
# gamma strength function.
# Plus other useful functions for shell model stuff.



def div0( a, b ):
  """ division function designed to ignore / 0, i.e. div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
  with np.errstate(divide='ignore', invalid='ignore'):
    c = np.true_divide( a, b )
    c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
  return c


def read_energy_levels(inputfile):
  # Reads levels from a KSHELL summary file, returns Nx3 matrix of 
  # [Ei, 2*Ji, parity], where E is absolute energy of level and parity 1=+, 0=-
  levels = []
  with open(inputfile, 'r') as f:
    lines = f.readlines()
    i_start = -1
    for i in range(len(lines)):
      if len(lines[i].split())<1: continue
      if lines[i].split()[0] == "Energy":
        i_start = i+4
        break
    for i in range(i_start, len(lines)):
      if len(lines[i].split())<1: break
      words = lines[i].split()
      num_parity = 1 if words[2] == "+" else -1
      levels.append([float(words[5]), float(words[1]), num_parity])
    return np.asarray(levels) 


def read_transition_strengths(inputfile, type="M1"):
  transitions = []
  with open(inputfile, 'r') as f:
    lines = f.readlines()
    i_start = -1
    for i in range(len(lines)):
      if len(lines[i].split())<1: continue
      if lines[i].split()[0] == "B({:s})".format(type):
        # print "hello"
        i_start = i+2
        break
    for i in range(i_start, len(lines)):
      # print lines[i]
      if len(lines[i].split())<1: break
      line = lines[i]
      # Returns
      # [2Ji, pi, Ei, 2Jf, pf, Ef, Eg, B(M1,i->f)]  --  beware that the summary file has opposite initial/final convention to this!
      # JEM: Changed 20170428 from [2Ji, Ei, 2Jf, Ef, Eg, B(M1,i->f)] 
      # print line[0:3], line[12:22], line[22:25], line[34:43], line[43:51], line[51:67]
      # print line[25]
      pi_str = line[25:27].strip()
      # print pi_str
      if pi_str == "+":
        pi = +1
      elif pi_str == "-":
        pi = -1
      else:
        raise Exception("From function read_transition_strengths: Could not assign initial parity. Read value: "+pi_str)
      pf_str = line[4:5].strip()
      if pf_str == "+":
        pf = +1
      elif pf_str == "-":
        pf = -1
      else:
        raise Exception("From function read_transition_strengths: Could not assign final parity Read value: "+pf_str)
      transitions.append([float(line[22:25]), pi, float(line[34:43]), float(line[0:3]), pi, float(line[12:22]), float(line[43:51]), float(line[67:83])])
    return np.asarray(transitions)


def read_calc_tbme(inputfile, tbme_template):
  # Reads tbme expectation values from a KSHELL summary file, returns NxM matrix of 
  # the M SPE+TBMEs for each level N
  # The tbme_template file is used to compare the ordering of <k1k2|V|k3k4;J> because
  # KSHELL can give either ordering of the k1k2, k3k4 pairs.
  calc_tbme = []
  with open(inputfile, 'r') as f:
    lines = f.readlines()
    i_start = -1
    for i in range(len(lines)):
      if len(lines[i].split())<1: continue
      if lines[i].split()[0] == "SPE&TBME":
        i_start = i+1
        break
    # Start looping through all eigenstates:
    i = i_start
    while True:
      if len(lines[i]) >= 30 and lines[i][0:30] == "    N          <E>    log-file":
        # Start reading all calc_tbme for current eigenstate
        i += 1
        words = lines[i].split()
        Nlevel = int(words[0])
        Elevel = float(words[1])
        i += 1
        words = lines[i].split()
        Nspe = int(words[-1])
        i += 1
        spe = np.zeros((Nspe,3)) # Rows of [N, x_i, <O_i>]
        for j in range(Nspe):
          words = lines[i].split()
          spe[j,:] = (int(words[0]), float(words[1]), float(words[2]))
          i += 1
        # i should now be line starting with TBME
        if not lines[i][0:5] == "TBMEs": raise Exception("Error from read_calc_tbme(). Error in summary file format?")
        words = lines[i].split()
        Ntbme = int(words[-1])
        i += 1
        tbme = np.zeros((Ntbme, 9))
        for j in range(Ntbme):
          words = lines[i].split()
          k1, k2, k3, k4, jj, iprty, ipn, xi, Oi = int(words[0]), int(words[1]), int(words[2]), int(words[3]), int(words[4]), int(words[5]), int(words[6]), float(words[7]), float(words[8])
          found_current_tbme = False
          for i_template in range(Ntbme):
            l1, l2, l3, l4, jj_t = int(tbme_template[i_template,0]), int(tbme_template[i_template,1]), int(tbme_template[i_template,2]), int(tbme_template[i_template,3]), int(tbme_template[i_template,4])
            if l1==k1 and l2==k2 and l3==k3 and l4==k4 and jj==jj_t:
              tbme[j,:] = (k1, k2, k3, k4, jj, iprty, ipn, xi, Oi)
              found_current_tbme = True
            elif l1==k3 and l2==k4 and l3==k1 and l4==k2 and jj==jj_t:
              tbme[j,:] = (k3, k4, k1, k2, jj, iprty, ipn, xi, Oi)
              found_current_tbme = True
          if not found_current_tbme: raise Exception("Error in read_calc_tbme(): Did not find matching quantum numbers in tbme_template for TBME no "+str(j))
          i += 1
        calc_tbme.append({"E":Elevel, "SPE":spe, "TBME":tbme})
      # end reading of one egeinstate
      i += 1
      if i >= len(lines): break

    return calc_tbme



def total_level_density(levels, bin_width, Ex_max):
  # 20170816: This function returns the total level density as a function of Ex.
  Nbins = int(np.ceil(Ex_max/bin_width)) # Make sure the number of bins cover the whole Ex region.
  bin_array = np.linspace(0,bin_width*Nbins,Nbins+1) # Array of lower bin edge energy values
  Egs = levels[0,0]
  rho_total, tmp = np.histogram(levels[:,0]-Egs, bins=bin_array)
  rho_total = rho_total/bin_width # To get density per MeV
  return rho_total


def strength_function_average_updated_definition(levels, transitions, Jpi_list, bin_width, Ex_min, Ex_max, type="M1"):
  # 20171009: Updated the way we average over Ex, J, pi to only count pixels with non-zero gSF.
  # 20170815: This function returns the strength function the way we now think is the correct way:
  # By taking only the partial level density corresponding to the specific (Ex, J, pi) pixel in the
  # calculation of the strength function, and then averaging over all three variables to produce
  # <f(Eg)>.
  # This code was first developed in the script strength_function_individual_Jpi.py
  Nbins = int(np.ceil(Ex_max/bin_width)) # Make sure the number of bins cover the whole Ex region.
  # print "Ex_max =", Ex_max, "Nbins*bin_width =", Nbins*bin_width
  bin_array = np.linspace(0,bin_width*Nbins,Nbins+1) # Array of lower bin edge energy values
  bin_array_middle = (bin_array[0:-1]+bin_array[1:])/2 # Array of middle bin values
  # Find index of first and last bin (lower bin edge) where we put counts.
  # It's important to not include the other Ex bins in the averaging later, because they contain zeros which will pull the average down.
  i_Exmin = int(np.floor(Ex_min/bin_width)) 
  i_Exmax = int(np.floor(Ex_max/bin_width))  

  prefactor = {"M1":  11.5473e-9, "E1": 1.047e-6}

  Egs = levels[0,0] # Read out the absolute ground state energy, so we can get relative energies later

  # Allocate matrices to store the summed B(M1) values for each pixel, and the number of transitions counted
  B_pixel_sum = np.zeros((Nbins,Nbins,len(Jpi_list)))
  B_pixel_count = np.zeros((Nbins,Nbins,len(Jpi_list)))

  # # Update 20170915: Realized a problem with summing vs averaging, adding this list to fix that
  # # Update 20170920: Realized that the fix was wrong, the original was correct.
  # Ex_already_seen = []
  # for i_Ex in range(Nbins):
  #   Ex_already_seen.append([])

  # Loop over all transitions and put in the correct pixel:
  for i_tr in range(len(transitions[:,0])):
    Ex = transitions[i_tr,2] - Egs
    # Check if transition is below Ex_max, skip if not
    if Ex < Ex_min or Ex >= Ex_max:
      continue

    # Get bin index for Eg and Ex (initial). Indices are defined with respect to the lower bin edge.
    i_Eg = int(np.floor(transitions[i_tr,6]/bin_width))
    i_Ex = int(np.floor(Ex/bin_width))

    # Read initial spin and parity of level:
    Ji = int(transitions[i_tr,0])
    pi = int(transitions[i_tr,1])
    # print pi
    # Get index for current [Ji,pi] combination in Jpi_list:
    # print "Tr =", transitions[i_tr,7]
    try:
      i_Jpi = Jpi_list.index([Ji,pi])
    except: 
      # print "no pass"
      continue
    # print i_Ex, i_Eg, i_Jpi

    # Add B(M1) value and increment count to pixel, respectively
    B_pixel_sum[i_Ex,i_Eg,i_Jpi] += transitions[i_tr,7]

    # 20170920: We thought this was more correct, but now think not.
    # if not Ex in Ex_already_seen[i_Eg]:
    #   B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1
    #   Ex_already_seen[i_Eg].append(Ex)

    B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1 # Original
  # print B_pixel_sum[B_pixel_sum>0].size, B_pixel_sum[B_pixel_sum>0]
  # print B_pixel_count[B_pixel_count>0].size, B_pixel_count[B_pixel_count>0]



  # Allocate (Ex,Jpi) matrix to store level density
  rho_ExJpi = np.zeros((Nbins,len(Jpi_list)))
  # Count number of levels for each (Ex, J, pi) pixel.
  for i_l in range(len(levels[:,0])):
    E, J, pi = levels[i_l]
    # Skip if level is outside range:
    if E-Egs >= Ex_max:
      continue
    i_Ex = int(np.floor((E-Egs)/bin_width))
    try:
      i_Jpi = Jpi_list.index([J,pi])
    except:
      continue
    rho_ExJpi[i_Ex,i_Jpi] += 1

  rho_ExJpi /= bin_width # Normalize to bin width, to get density in MeV^-1
  # print rho_ExJpi


  # Calculate gamma strength functions for each Ex, J, pi individually, using the partial level density for each J, pi.
  gSF = np.zeros((Nbins,Nbins,len(Jpi_list)))
  a = prefactor[type] # mu_N^-2 MeV^-2, conversion constant
  for i_Jpi in range(len(Jpi_list)):
    for i_Ex in range(Nbins):
                     # a *           <B(M1; Eg, Ex, J, pi)>                                       * rho(Ex, J, pi)  
      # if rho_ExJpi[i_Ex, i_Jpi] > 0:
      #   print "Ding rho, ", "Ex=", bin_array[i_Ex], "Jpi=", Jpi_list[i_Jpi], rho_ExJpi[i_Ex, i_Jpi]
      # if B_pixel_sum[i_Ex,:, i_Jpi].any() > 0:
      #   print "Ding B, ", "Ex=", bin_array[i_Ex], "Jpi=", Jpi_list[i_Jpi], rho_ExJpi[i_Ex, i_Jpi], B_pixel_sum[i_Ex,:, i_Jpi].nonzero()
      gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * rho_ExJpi[i_Ex, i_Jpi]
      # print gSF[i_Ex,gSF[i_Ex,:,i_Jpi]>0,i_Jpi]
      # # DEBUG:
      # for i_Eg in range(i_Ex):
      #   if gSF[i_Ex,i_Eg,i_Jpi] > 0:
      #     print "Ding, ", i_Ex, i_Eg, i_Jpi, gSF[i_Ex,:,i_Jpi]
      # TODO if the below smoothing is used: Check the form of the level density. When I tried sigma=3 it looked like the gSF dropped a bit, indicating that the level density was lowered. Maybe smoothe it some other way, such as drawing straight lines?
      # gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * smoothe(rho_ExJpi[:, i_Jpi], sigma=1)[i_Ex] # Smoothing the initial level density by a gaussian convolution

  # Return the average gSF(Eg) over all (Ex,J,pi)

  # return gSF[i_Exmin:i_Exmax+1,:,:].mean(axis=(0,2))
  # Update 20171009: Took proper care to only average over the non-zero f(Eg,Ex,J,pi) pixels:
  gSF_currentExrange = gSF[i_Exmin:i_Exmax+1,:,:]
  gSF_ExJpiavg = div0(gSF_currentExrange.sum(axis=(0,2)), (gSF_currentExrange!=0).sum(axis=(0,2)))
  return gSF_ExJpiavg




# def strength_function_average_updated_definition_debug(levels, transitions, Jpi_list, bin_width, Ex_min, Ex_max, type="M1"):
#   # 20171009: Updated the way we average over Ex, J, pi to only count pixels with non-zero gSF.
#   # 20170815: This function returns the strength function the way we now think is the correct way:
#   # By taking only the partial level density corresponding to the specific (Ex, J, pi) pixel in the
#   # calculation of the strength function, and then averaging over all three variables to produce
#   # <f(Eg)>.
#   # This code was first developed in the script strength_function_individual_Jpi.py
#   Nbins = int(np.ceil(Ex_max/bin_width)) # Make sure the number of bins cover the whole Ex region.
#   # print "Ex_max =", Ex_max, "Nbins*bin_width =", Nbins*bin_width
#   bin_array = np.linspace(0,bin_width*Nbins,Nbins+1) # Array of lower bin edge energy values
#   bin_array_middle = (bin_array[0:-1]+bin_array[1:])/2 # Array of middle bin values
#   # Find index of first and last bin (lower bin edge) where we put counts.
#   # It's important to not include the other Ex bins in the averaging later, because they contain zeros which will pull the average down.
#   i_Exmin = int(np.floor(Ex_min/bin_width)) 
#   i_Exmax = int(np.floor(Ex_max/bin_width))  

#   prefactor = {"M1":  11.5473e-9, "E1": 1.047e-6}

#   Egs = levels[0,0] # Read out the absolute ground state energy, so we can get relative energies later

#   # Allocate matrices to store the summed B(M1) values for each pixel, and the number of transitions counted
#   B_pixel_sum = np.zeros((Nbins,Nbins,len(Jpi_list)))
#   B_pixel_count = np.zeros((Nbins,Nbins,len(Jpi_list)))

#   # # Update 20170915: Realized a problem with summing vs averaging, adding this list to fix that
#   # # Update 20170920: Realized that the fix was wrong, the original was correct.
#   # Ex_already_seen = []
#   # for i_Ex in range(Nbins):
#   #   Ex_already_seen.append([])

#   # Loop over all transitions and put in the correct pixel:
#   for i_tr in range(len(transitions[:,0])):
#     Ex = transitions[i_tr,2] - Egs
#     # Check if transition is below Ex_max, skip if not
#     if Ex < Ex_min or Ex >= Ex_max:
#       continue

#     # Get bin index for Eg and Ex (initial). Indices are defined with respect to the lower bin edge.
#     i_Eg = int(np.floor(transitions[i_tr,6]/bin_width))
#     i_Ex = int(np.floor(Ex/bin_width))

#     # Read initial spin and parity of level:
#     Ji = int(transitions[i_tr,0])
#     pi = int(transitions[i_tr,1])
#     print pi
#     # Get index for current [Ji,pi] combination in Jpi_list:
#     print "Tr =", transitions[i_tr,7]
#     try:
#       i_Jpi = Jpi_list.index([Ji,pi])
#     except: 
#       print "no pass"
#       continue
#     print i_Ex, i_Eg, i_Jpi

#     # Add B(M1) value and increment count to pixel, respectively
#     B_pixel_sum[i_Ex,i_Eg,i_Jpi] += transitions[i_tr,7]

#     # 20170920: We thought this was more correct, but now think not.
#     # if not Ex in Ex_already_seen[i_Eg]:
#     #   B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1
#     #   Ex_already_seen[i_Eg].append(Ex)

#     B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1 # Original
#   print B_pixel_sum[B_pixel_sum>0].size, B_pixel_sum[B_pixel_sum>0]
#   print B_pixel_count[B_pixel_count>0].size, B_pixel_count[B_pixel_count>0]



#   # Allocate (Ex,Jpi) matrix to store level density
#   rho_ExJpi = np.zeros((Nbins,len(Jpi_list)))
#   # Count number of levels for each (Ex, J, pi) pixel.
#   for i_l in range(len(levels[:,0])):
#     E, J, pi = levels[i_l]
#     # Skip if level is outside range:
#     if E-Egs >= Ex_max:
#       continue
#     i_Ex = int(np.floor((E-Egs)/bin_width))
#     try:
#       i_Jpi = Jpi_list.index([J,pi])
#     except:
#       continue
#     rho_ExJpi[i_Ex,i_Jpi] += 1

#   rho_ExJpi /= bin_width # Normalize to bin width, to get density in MeV^-1
#   print rho_ExJpi


#   # Calculate gamma strength functions for each Ex, J, pi individually, using the partial level density for each J, pi.
#   gSF = np.zeros((Nbins,Nbins,len(Jpi_list)))
#   a = prefactor[type] # mu_N^-2 MeV^-2, conversion constant
#   for i_Jpi in range(len(Jpi_list)):
#     for i_Ex in range(Nbins):
#                      # a *           <B(M1; Eg, Ex, J, pi)>                                       * rho(Ex, J, pi)  
#       if rho_ExJpi[i_Ex, i_Jpi] > 0:
#         print "Ding rho, ", "Ex=", bin_array[i_Ex], "Jpi=", Jpi_list[i_Jpi], rho_ExJpi[i_Ex, i_Jpi]
#       if B_pixel_sum[i_Ex,:, i_Jpi].any() > 0:
#         print "Ding B, ", "Ex=", bin_array[i_Ex], "Jpi=", Jpi_list[i_Jpi], rho_ExJpi[i_Ex, i_Jpi], B_pixel_sum[i_Ex,:, i_Jpi].nonzero()
#       gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * rho_ExJpi[i_Ex, i_Jpi]
#       # print gSF[i_Ex,gSF[i_Ex,:,i_Jpi]>0,i_Jpi]
#       # # DEBUG:
#       # for i_Eg in range(i_Ex):
#       #   if gSF[i_Ex,i_Eg,i_Jpi] > 0:
#       #     print "Ding, ", i_Ex, i_Eg, i_Jpi, gSF[i_Ex,:,i_Jpi]
#       # TODO if the below smoothing is used: Check the form of the level density. When I tried sigma=3 it looked like the gSF dropped a bit, indicating that the level density was lowered. Maybe smoothe it some other way, such as drawing straight lines?
#       # gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * smoothe(rho_ExJpi[:, i_Jpi], sigma=1)[i_Ex] # Smoothing the initial level density by a gaussian convolution

#   # Return the average gSF(Eg) over all (Ex,J,pi)

#   # return gSF[i_Exmin:i_Exmax+1,:,:].mean(axis=(0,2))
#   # Update 20171009: Took proper care to only average over the non-zero f(Eg,Ex,J,pi) pixels:
#   gSF_currentExrange = gSF[i_Exmin:i_Exmax+1,:,:]
#   gSF_ExJpiavg = div0(gSF_currentExrange.sum(axis=(0,2)), (gSF_currentExrange!=0).sum(axis=(0,2)))
#   return gSF_ExJpiavg






def strength_function_average_updated_definition_brute_avg(levels, transitions, Jpi_list, bin_width, Ex_min, Ex_max, type="M1"):
  # 20171009: This is a check to see that the averaging we do over Ex, J, pi really only counts
  # the bins with a non-zero gSF.
  # 20170815: This function returns the strength function the way we now think is the correct way:
  # By taking only the partial level density corresponding to the specific (Ex, J, pi) pixel in the
  # calculation of the strength function, and then averaging over all three variables to produce
  # <f(Eg)>.
  # This code was first developed in the script strength_function_individual_Jpi.py
  Nbins = int(np.ceil(Ex_max/bin_width)) # Make sure the number of bins cover the whole Ex region.
  # print "Ex_max =", Ex_max, "Nbins*bin_width =", Nbins*bin_width
  bin_array = np.linspace(0,bin_width*Nbins,Nbins+1) # Array of lower bin edge energy values
  bin_array_middle = (bin_array[0:-1]+bin_array[1:])/2 # Array of middle bin values
  # Find index of first and last bin (lower bin edge) where we put counts.
  # It's important to not include the other Ex bins in the averaging later, because they contain zeros which will pull the average down.
  i_Exmin = int(np.floor(Ex_min/bin_width)) 
  i_Exmax = int(np.floor(Ex_max/bin_width))  

  prefactor = {"M1":  11.5473e-9, "E1": 1.047e-6}

  Egs = levels[0,0] # Read out the absolute ground state energy, so we can get relative energies later

  # Allocate matrices to store the summed B(M1) values for each pixel, and the number of transitions counted
  B_pixel_sum = np.zeros((Nbins,Nbins,len(Jpi_list)))
  B_pixel_count = np.zeros((Nbins,Nbins,len(Jpi_list)))

  # # Update 20170915: Realized a problem with summing vs averaging, adding this list to fix that
  # # Update 20170920: Realized that the fix was wrong, the original was correct.
  # Ex_already_seen = []
  # for i_Ex in range(Nbins):
  #   Ex_already_seen.append([])

  # Loop over all transitions and put in the correct pixel:
  for i_tr in range(len(transitions[:,0])):
    Ex = transitions[i_tr,2] - Egs
    # Check if transition is below Ex_max, skip if not
    if Ex < Ex_min or Ex >= Ex_max:
      continue

    # Get bin index for Eg and Ex (initial). Indices are defined with respect to the lower bin edge.
    i_Eg = int(np.floor(transitions[i_tr,6]/bin_width))
    i_Ex = int(np.floor(Ex/bin_width))

    # Read initial spin and parity of level:
    Ji = int(transitions[i_tr,0])
    pi = int(transitions[i_tr,1])
    # Get index for current [Ji,pi] combination in Jpi_list:
    try:
      i_Jpi = Jpi_list.index([Ji,pi])
    except: 
      continue

    # Add B(M1) value and increment count to pixel, respectively
    B_pixel_sum[i_Ex,i_Eg,i_Jpi] += transitions[i_tr,7]

    # 20170920: We thought this was more correct, but now think not.
    # if not Ex in Ex_already_seen[i_Eg]:
    #   B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1
    #   Ex_already_seen[i_Eg].append(Ex)

    B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1 # Original



  # Allocate (Ex,Jpi) matrix to store level density
  rho_ExJpi = np.zeros((Nbins,len(Jpi_list)))
  # Count number of levels for each (Ex, J, pi) pixel.
  for i_l in range(len(levels[:,0])):
    E, J, pi = levels[i_l]
    # Skip if level is outside range:
    if E-Egs >= Ex_max:
      continue
    i_Ex = int(np.floor((E-Egs)/bin_width))
    try:
      i_Jpi = Jpi_list.index([J,pi])
    except:
      continue
    rho_ExJpi[i_Ex,i_Jpi] += 1

  rho_ExJpi /= bin_width # Normalize to bin width, to get density in MeV^-1


  # Calculate gamma strength functions for each Ex, J, pi individually, using the partial level density for each J, pi.
  gSF = np.zeros((Nbins,Nbins,len(Jpi_list)))
  a = prefactor[type] # mu_N^-2 MeV^-2, conversion constant
  for i_Jpi in range(len(Jpi_list)):
    for i_Ex in range(Nbins):
                     # a *           <B(M1; Eg, Ex, J, pi)>                                       * rho(Ex, J, pi)  
      gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * rho_ExJpi[i_Ex, i_Jpi]
      # TODO if the below smoothing is used: Check the form of the level density. When I tried sigma=3 it looked like the gSF dropped a bit, indicating that the level density was lowered. Maybe smoothe it some other way, such as drawing straight lines?
      # gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * smoothe(rho_ExJpi[:, i_Jpi], sigma=1)[i_Ex] # Smoothing the initial level density by a gaussian convolution

  # Return the average gSF(Eg) over all (Ex,J,pi)
  # Test of brute force way, to check:
  gSF_ExJpiavg = np.zeros(Nbins)
  for i_Eg in range(Nbins):
    sum = 0
    counter = 0
    for i_Ex in range(i_Exmin, i_Exmax+1):
      for i_Jpi in range(len(Jpi_list)):
        if gSF[i_Ex, i_Eg, i_Jpi] > 0:
          sum += gSF[i_Ex, i_Eg, i_Jpi]
          counter += 1
    if counter > 0:
      gSF_ExJpiavg[i_Eg] = sum / counter


  return gSF_ExJpiavg




def strength_function_average_updated_definition_naive_avg(levels, transitions, Jpi_list, bin_width, Ex_min, Ex_max, type="M1"):
  # 20171009: Realized that this was a naive way of averaging over Ex, J and pi.
  # 20170815: This function returns the strength function the way we now think is the correct way:
  # By taking only the partial level density corresponding to the specific (Ex, J, pi) pixel in the
  # calculation of the strength function, and then averaging over all three variables to produce
  # <f(Eg)>.
  # This code was first developed in the script strength_function_individual_Jpi.py
  Nbins = int(np.ceil(Ex_max/bin_width)) # Make sure the number of bins cover the whole Ex region.
  # print "Ex_max =", Ex_max, "Nbins*bin_width =", Nbins*bin_width
  bin_array = np.linspace(0,bin_width*Nbins,Nbins+1) # Array of lower bin edge energy values
  bin_array_middle = (bin_array[0:-1]+bin_array[1:])/2 # Array of middle bin values
  # Find index of first and last bin (lower bin edge) where we put counts.
  # It's important to not include the other Ex bins in the averaging later, because they contain zeros which will pull the average down.
  i_Exmin = int(np.floor(Ex_min/bin_width)) 
  i_Exmax = int(np.floor(Ex_max/bin_width))  

  prefactor = {"M1":  11.5473e-9, "E1": 1.047e-6}

  Egs = levels[0,0] # Read out the absolute ground state energy, so we can get relative energies later

  # Allocate matrices to store the summed B(M1) values for each pixel, and the number of transitions counted
  B_pixel_sum = np.zeros((Nbins,Nbins,len(Jpi_list)))
  B_pixel_count = np.zeros((Nbins,Nbins,len(Jpi_list)))

  # # Update 20170915: Realized a problem with summing vs averaging, adding this list to fix that
  # # Update 20170920: Realized that the fix was wrong, the original was correct.
  # Ex_already_seen = []
  # for i_Ex in range(Nbins):
  #   Ex_already_seen.append([])

  # Loop over all transitions and put in the correct pixel:
  for i_tr in range(len(transitions[:,0])):
    Ex = transitions[i_tr,2] - Egs
    # Check if transition is below Ex_max, skip if not
    if Ex < Ex_min or Ex >= Ex_max:
      continue

    # Get bin index for Eg and Ex (initial). Indices are defined with respect to the lower bin edge.
    i_Eg = int(np.floor(transitions[i_tr,6]/bin_width))
    i_Ex = int(np.floor(Ex/bin_width))

    # Read initial spin and parity of level:
    Ji = int(transitions[i_tr,0])
    pi = int(transitions[i_tr,1])
    # Get index for current [Ji,pi] combination in Jpi_list:
    try:
      i_Jpi = Jpi_list.index([Ji,pi])
    except: 
      continue

    # Add B(M1) value and increment count to pixel, respectively
    B_pixel_sum[i_Ex,i_Eg,i_Jpi] += transitions[i_tr,7]

    # 20170920: We thought this was more correct, but now think not.
    # if not Ex in Ex_already_seen[i_Eg]:
    #   B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1
    #   Ex_already_seen[i_Eg].append(Ex)

    B_pixel_count[i_Ex,i_Eg,i_Jpi] += 1 # Original



  # Allocate (Ex,Jpi) matrix to store level density
  rho_ExJpi = np.zeros((Nbins,len(Jpi_list)))
  # Count number of levels for each (Ex, J, pi) pixel.
  for i_l in range(len(levels[:,0])):
    E, J, pi = levels[i_l]
    # Skip if level is outside range:
    if E-Egs >= Ex_max:
      continue
    i_Ex = int(np.floor((E-Egs)/bin_width))
    try:
      i_Jpi = Jpi_list.index([J,pi])
    except:
      continue
    rho_ExJpi[i_Ex,i_Jpi] += 1

  rho_ExJpi /= bin_width # Normalize to bin width, to get density in MeV^-1


  # Calculate gamma strength functions for each Ex, J, pi individually, using the partial level density for each J, pi.
  gSF = np.zeros((Nbins,Nbins,len(Jpi_list)))
  a = prefactor[type] # mu_N^-2 MeV^-2, conversion constant
  for i_Jpi in range(len(Jpi_list)):
    for i_Ex in range(Nbins):
                     # a *           <B(M1; Eg, Ex, J, pi)>                                       * rho(Ex, J, pi)  
      gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * rho_ExJpi[i_Ex, i_Jpi]
      # TODO if the below smoothing is used: Check the form of the level density. When I tried sigma=3 it looked like the gSF dropped a bit, indicating that the level density was lowered. Maybe smoothe it some other way, such as drawing straight lines?
      # gSF[i_Ex,:,i_Jpi] = a * div0(B_pixel_sum[i_Ex,:,i_Jpi], B_pixel_count[i_Ex,:,i_Jpi]) * smoothe(rho_ExJpi[:, i_Jpi], sigma=1)[i_Ex] # Smoothing the initial level density by a gaussian convolution

  # Return the average gSF(Eg) over all (Ex,J,pi)

  return gSF[i_Exmin:i_Exmax+1,:,:].mean(axis=(0,2))






def smoothe(array, sigma):
  from scipy.stats import norm
  # Smoothe a one-dimensional array by convoluting with a gaussian of width sigma
  N = len(array)
  index_array = np.linspace(0,N-1, N)
  matrix = np.zeros((N,N))
  for i in range(N):
      try:
          matrix[i] = array[i]*norm.pdf(index_array, loc=index_array[i], scale=sigma) 
      except IndexError:
          pass
  # plt.matshow(matrix)
  # plt.show()
  array_smoothed = matrix.sum(axis=0)
  return array_smoothed




def strength_function_average_updated_definition_Jpiaveraging(levels, transitions, Jpi_list, bin_width, Ex_min, Ex_max, type="M1"):
  # 20170815: This function returns the strength function the way we now think is the correct way:
  # By taking only the partial level density corresponding to the specific (Ex, J, pi) pixel in the
  # calculation of the strength function, and then averaging over all three variables to produce
  # <f(Eg)>.
  # 20170821: Another update to this version: To blind ourselves to the same degree as in the 
  # data, we use an average partial level density for each Ex bin to multiply the <B>.
  # This code was first developed in the script strength_function_individual_Jpi.py
  Nbins = int(np.ceil(Ex_max/bin_width)) # Make sure the number of bins cover the whole Ex region.
  # print "Ex_max =", Ex_max, "Nbins*bin_width =", Nbins*bin_width
  bin_array = np.linspace(0,bin_width*Nbins,Nbins+1) # Array of lower bin edge energy values
  bin_array_middle = (bin_array[0:-1]+bin_array[1:])/2 # Array of middle bin values
  # Find index of first and last bin (lower bin edge) where we put counts.
  # It's important to not include the other Ex bins in the averaging later, because they contain zeros which will pull the average down.
  i_Exmin = int(np.floor(Ex_min/bin_width)) 
  i_Exmax = int(np.floor(Ex_max/bin_width))  

  Egs = levels[0,0] # Read out the absolute ground state energy, so we can get relative energies later

  # Allocate matrices to store the summed B(M1) values for each pixel, and the number of transitions counted
  B_pixel_sum = np.zeros((Nbins,Nbins))
  B_pixel_count = np.zeros((Nbins,Nbins))

  # Loop over all transitions and put in the correct pixel:
  for i_tr in range(len(transitions[:,0])):
    Ex = transitions[i_tr,2] - Egs
    # Check if transition is below Ex_max, skip if not
    if Ex < Ex_min or Ex >= Ex_max:
      continue

    # Get bin index for Eg and Ex (initial). Indices are defined with respect to the lower bin edge.
    i_Eg = int(np.floor(transitions[i_tr,6]/bin_width))
    i_Ex = int(np.floor(Ex/bin_width))

    # Read initial spin and parity of level:
    Ji = int(transitions[i_tr,0])
    pi = int(transitions[i_tr,1])
    if not [Ji,pi] in Jpi_list:
      continue

    # Add B(M1) value and increment count to pixel, respectively
    B_pixel_sum[i_Ex,i_Eg] += transitions[i_tr,7]
    B_pixel_count[i_Ex,i_Eg] += 1


  # plt.pcolormesh(bin_array_middle, bin_array_middle, B_pixel_sum[:,:,20],norm=LogNorm())
  # plt.show()
  # sys.exit(0)
  # OK -- the B summing seems to work correctly.



  # Allocate Ex array to store level density
  rho_ExJpiavg = np.zeros(Nbins)
  # Count number of levels for each Ex pixel, but only selected (J,pi) values
  for i_l in range(len(levels[:,0])):
    E, J, pi = levels[i_l]
    # Skip if level is outside range:
    if E-Egs >= Ex_max:
      continue
    if not [J,pi] in Jpi_list:
      continue
    i_Ex = int(np.floor((E-Egs)/bin_width))
    rho_ExJpiavg[i_Ex] += 1

  rho_ExJpiavg /= (bin_width*len(Jpi_list)) # Normalize to bin width and number of (J,pi) combinations, to get density in MeV^-1


  # Calculate gamma strength functions for each Ex individually, using the average partial level density
  gSF = np.zeros((Nbins,Nbins))
  a = 11.5473e-9 # mu_N^-2 MeV^-2, conversion constant
  for i_Ex in range(Nbins):
                   # a *           <B(M1; Eg, Ex, J, pi)>                                       * rho(Ex, J, pi)  
    gSF[i_Ex,:] = a * div0(B_pixel_sum[i_Ex,:], B_pixel_count[i_Ex,:]) * rho_ExJpiavg[i_Ex]

  # Return the average gSF(Eg) over all Ex
  return gSF[i_Exmin:i_Exmax+1,:].mean(axis=0)













def level_density_matrix(inputfile, bin_width=0.2, Emax=12, Ex_low=5, Ex_high=8):
  
  levels = read_energy_levels(inputfile)

  # Set bin width and range
  Nbins = int(Emax/bin_width)
  Emax_adjusted = bin_width*Nbins
  bins_Ex = np.linspace(0,Emax_adjusted,Nbins+1)
  bins_Ex_middle = (bins_Ex[0:-1]+bins_Ex[1:])/2

  bins_J = np.linspace(0, levels[:,1].max()/2, int(levels[:,1].max()/2)+1)
  print bins_J

  Egs = levels[0,0]
  rho_total, tmp = np.histogram(levels[:,0]-Egs, bins=bins_Ex)
  rho_total = rho_total/bin_width # To get density per MeV

  matrix, xedges, yedges = np.histogram2d(levels[:,1]/2, levels[:,0]-Egs, bins=[bins_J,bins_Ex])
  return matrix, xedges, yedges



def level_density_matrix_parity_decomposed(inputfile, bin_width=0.2, Emax=12, Ex_low=5, Ex_high=8):
  
  levels = read_energy_levels(inputfile)

  # Set bin width and range
  Nbins = int(Emax/bin_width)
  Emax_adjusted = bin_width*Nbins
  bins_Ex = np.linspace(0,Emax_adjusted,Nbins+1)
  bins_Ex_middle = (bins_Ex[0:-1]+bins_Ex[1:])/2

  bins_J = np.linspace(0, levels[:,1].max()/2, int(levels[:,1].max()/2)+1)
  print bins_J

  Egs = levels[0,0]
  # rho_total, tmp = np.histogram(levels[:,0]-Egs, bins=bins_Ex)
  # rho_total = rho_total/bin_width # To get density per MeV

  matrix_plus, xedges, yedges = np.histogram2d(levels[levels[:,2]>0,1]/2, levels[levels[:,2]>0,0]-Egs, bins=[bins_J,bins_Ex])
  matrix_minus, xedges, yedges = np.histogram2d(levels[levels[:,2]<0,1]/2, levels[levels[:,2]<0,0]-Egs, bins=[bins_J,bins_Ex])
  return matrix_plus, matrix_minus, xedges, yedges



def read_interaction_file(filename):
  """
  Read interaction file in KSHELL .snt format
  returns SPE, TBME, msdict, core
  """
  SPE = []
  TBME = []
  core = {}
  msdict = {}
  # Model space dict formatted like
  # msdict = {1: {"n": 0, "l": 3, "j": 7, "Tz": -1},
  #           2: {"n": 0, "l": 3, "j": 5, "Tz": -1},
  #           3: {"n": 1, "l": 1, "j": 3, "Tz": -1},
  #           4: {"n": 1, "l": 1, "j": 1, "Tz": -1},
  #           5: {"n": 0, "l": 3, "j": 5, "Tz":  1},
  #           6: {"n": 1, "l": 1, "j": 3, "Tz":  1},
  #           7: {"n": 1, "l": 1, "j": 1, "Tz":  1},
  #           8: {"n": 0, "l": 4, "j": 9, "Tz":  1}}
  infile = open(filename, "r")
  while True:
    line = infile.readline()
    if not line:
      break
    if line[0:13] == "! model space":
      line = infile.readline() # Line containing Np, Nn, Z0, A0-Z0
      words = line.split()
      core = {"A":int(words[2])+int(words[3]), "Z":int(words[2])}
      Np = int(words[0])
      Nn = int(words[1])
      # Fill msdict with proton and then neutron orbitals:
      for i in range(Np+Nn):
        line = infile.readline()
        words = line.split()
        msdict[int(words[0])] = {"n":int(words[1]), "l":int(words[2]), "j":int(words[3]), "Tz":int(words[4])}
    if line[0:23] == "!  i  j     <i|H(1b)|j>":
      line = infile.readline() # Line containing Nspe
      words = line.split()
      Nspe = int(words[0])
      for i in range(Nspe):
        line = infile.readline()
        words = line.split()
        SPE.append(float(words[2]))
      SPE = np.array(SPE)
    if line[0:6] == "! TBME":
      line = infile.readline() # Line containing Ntbme
      words = line.split()
      Ntbme = int(words[0])
      for i in range(Ntbme):
        line = infile.readline()
        words = line.split()
        TBME.append([int(words[0]), int(words[1]), int(words[2]), int(words[3]), int(words[4]), float(words[5])])
      TBME = np.array(TBME) 




  return SPE, TBME, msdict, core



def write_interaction_file(filename, SPEs, TBMEs, model_space, core, comments="", mass_scaling=False, scaling_A0=-1, scaling_p=-0.300000):
  outfile = "! interaction file written by write_interaction_file() from statistical.py\n"
  outfile += "! User-added comments, if any:\n"
  outfile += "! " + comments + "\n"

  # Model space formatted like
  #   #      n    l   2j   2tz
  # [[1,     0,   3,   7,  -1],
  #  [2,     0,   3,   5,  -1],
  #  [3,     1,   1,   3,  -1],
  #  [4,     1,   1,   1,  -1],
  #  [5,     0,   3,   5,   1],
  #  [6,     1,   1,   3,   1],
  #  [7,     1,   1,   1,   1],
  #  [8,     0,   4,   9,   1]]
  outfile += "! model space\n" 
  # Write a line with (num proton orbitals)  (num neutron orbitals)  (Z_core)   (N_core)
  outfile += " {:3d} {:3d}  {:3d} {:3d}\n".format(len(np.where(model_space[:,4]<0)[0]), len(np.where(model_space[:,4]>0)[0]), core["Z"], core["A"]-core["Z"])
  for i in range(len(model_space[:,0])):
    outfile += " {:2d}  {:2d}  {:2d}  {:2d}  {:2d}\n".format(model_space[i,0],model_space[i,1],model_space[i,2],model_space[i,3],model_space[i,4])

  # SPEs as one-dim array
  outfile += "! interaction\n"
  outfile += "!  i  j     <i|H(1b)|j>\n"
  outfile += "  {:d}   {:d}\n".format(len(SPEs), 0) # Unsure what the last zero signifies.
  for i in range(len(SPEs)):
    outfile += "{:3d} {:3d} {:16.8f}\n".format(model_space[i,0], model_space[i,0], SPEs[i])

  outfile += "! TBME \n"
  if mass_scaling:
    if scaling_A0 < 0:
      raise Exception("Must specify scaling_A0")
    outfile += "    {:5d}   1  {:d} {:f} \n".format(len(TBMEs[:,0]),scaling_A0, scaling_p)
  else:
    outfile += "     {:5d}   {:3d}\n".format(len(TBMEs[:,0]), 0) # Unsure about last zero
  for i in range(len(TBMEs[:,0])):
    outfile += "{:3d} {:3d} {:3d} {:3d}  {:3d}  {:16.8f}\n".format(int(TBMEs[i,0]),int(TBMEs[i,1]),int(TBMEs[i,2]),int(TBMEs[i,3]),int(TBMEs[i,4]),TBMEs[i,5])


  with open(filename, 'w') as f:
    f.write(outfile)

  return True


def write_interaction_file_msdict(filename, SPEs, TBMEs, model_space, core, comments="", mass_scaling=False, scaling_A0=-1, scaling_p=-0.300000):
  outfile = "! interaction file written by write_interaction_file() from statistical.py\n"
  outfile += "! User-added comments, if any:\n"
  outfile += "! " + comments + "\n"

  # Model space formatted like
  # ms = {1: {"n": 0, "l": 3, "j": 7, "Tz": -1},
  #      2: {"n": 0, "l": 3, "j": 5, "Tz": -1},
  #      3: {"n": 1, "l": 1, "j": 3, "Tz": -1},
  #      4: {"n": 1, "l": 1, "j": 1, "Tz": -1},
  #      5: {"n": 0, "l": 3, "j": 5, "Tz":  1},
  #      6: {"n": 1, "l": 1, "j": 3, "Tz":  1},
  #      7: {"n": 1, "l": 1, "j": 1, "Tz":  1},
  #      8: {"n": 0, "l": 4, "j": 9, "Tz":  1}}

  # Convert it to the format 
  #   #      n    l   2j   2tz
  # [[1,     0,   3,   7,  -1],
  #  [2,     0,   3,   5,  -1],
  #  [3,     1,   1,   3,  -1],
  #  [4,     1,   1,   1,  -1],
  #  [5,     0,   3,   5,   1],
  #  [6,     1,   1,   3,   1],
  #  [7,     1,   1,   1,   1],
  #  [8,     0,   4,   9,   1]]  
  model_space_list = []
  for i in model_space.keys():
    model_space_list.append([i,model_space[i]["n"],model_space[i]["l"],model_space[i]["j"],model_space[i]["Tz"]])
  model_space = np.array(model_space_list)

  outfile += "! model space\n" 
  # Write a line with (num proton orbitals)  (num neutron orbitals)  (Z_core)   (N_core)
  outfile += " {:3d} {:3d}  {:3d} {:3d}\n".format(len(np.where(model_space[:,4]<0)[0]), len(np.where(model_space[:,4]>0)[0]), core["Z"], core["A"]-core["Z"])
  for i in range(len(model_space[:,0])):
    outfile += " {:2d}  {:2d}  {:2d}  {:2d}  {:2d}\n".format(model_space[i,0],model_space[i,1],model_space[i,2],model_space[i,3],model_space[i,4])

  # SPEs as one-dim array
  outfile += "! interaction\n"
  outfile += "!  i  j     <i|H(1b)|j>\n"
  outfile += "  {:d}   {:d}\n".format(len(SPEs), 0) # Unsure what the last zero signifies.
  for i in range(len(SPEs)):
    outfile += "{:3d} {:3d} {:16.8f}\n".format(model_space[i,0], model_space[i,0], SPEs[i])

  outfile += "! TBME \n"
  if mass_scaling:
    if scaling_A0 < 0:
      raise Exception("Must specify scaling_A0")
    outfile += "    {:5d}   1  {:d} {:f} \n".format(len(TBMEs[:,0]),scaling_A0, scaling_p)
  else:
    outfile += "     {:5d}   {:3d}\n".format(len(TBMEs[:,0]), 0) 
  for i in range(len(TBMEs[:,0])):
    outfile += "{:3d} {:3d} {:3d} {:3d}  {:3d}  {:16.8f}\n".format(int(TBMEs[i,0]),int(TBMEs[i,1]),int(TBMEs[i,2]),int(TBMEs[i,3]),int(TBMEs[i,4]),TBMEs[i,5])


  with open(filename, 'w') as f:
    f.write(outfile)

  return True

def write_interaction_file_msdict(filename, SPEs, TBMEs, model_space, core, comments="", mass_scaling=False, scaling_A0=-1, scaling_p=-0.300000):
  outfile = "! interaction file written by write_interaction_file() from statistical.py\n"
  outfile += "! User-added comments, if any:\n"
  outfile += "! " + comments + "\n"

  # Model space formatted like
  # ms = {1: {"n": 0, "l": 3, "j": 7, "Tz": -1},
  #      2: {"n": 0, "l": 3, "j": 5, "Tz": -1},
  #      3: {"n": 1, "l": 1, "j": 3, "Tz": -1},
  #      4: {"n": 1, "l": 1, "j": 1, "Tz": -1},
  #      5: {"n": 0, "l": 3, "j": 5, "Tz":  1},
  #      6: {"n": 1, "l": 1, "j": 3, "Tz":  1},
  #      7: {"n": 1, "l": 1, "j": 1, "Tz":  1},
  #      8: {"n": 0, "l": 4, "j": 9, "Tz":  1}}

  # Convert it to the format 
  #   #      n    l   2j   2tz
  # [[1,     0,   3,   7,  -1],
  #  [2,     0,   3,   5,  -1],
  #  [3,     1,   1,   3,  -1],
  #  [4,     1,   1,   1,  -1],
  #  [5,     0,   3,   5,   1],
  #  [6,     1,   1,   3,   1],
  #  [7,     1,   1,   1,   1],
  #  [8,     0,   4,   9,   1]]  
  model_space_list = []
  for i in model_space.keys():
    model_space_list.append([i,model_space[i]["n"],model_space[i]["l"],model_space[i]["j"],model_space[i]["Tz"]])
  model_space = np.array(model_space_list)

  outfile += "! model space\n" 
  # Write a line with (num proton orbitals)  (num neutron orbitals)  (Z_core)   (N_core)
  outfile += " {:3d} {:3d}  {:3d} {:3d}\n".format(len(np.where(model_space[:,4]<0)[0]), len(np.where(model_space[:,4]>0)[0]), core["Z"], core["A"]-core["Z"])
  for i in range(len(model_space[:,0])):
    outfile += " {:2d}  {:2d}  {:2d}  {:2d}  {:2d}\n".format(model_space[i,0],model_space[i,1],model_space[i,2],model_space[i,3],model_space[i,4])

  # SPEs as one-dim array
  outfile += "! interaction\n"
  outfile += "!  i  j     <i|H(1b)|j>\n"
  outfile += "  {:d}   {:d}\n".format(len(SPEs), 0) # Unsure what the last zero signifies.
  for i in range(len(SPEs)):
    outfile += "{:3d} {:3d} {:16.8f}\n".format(model_space[i,0], model_space[i,0], SPEs[i])

  outfile += "! TBME \n"
  if mass_scaling:
    if scaling_A0 < 0:
      raise Exception("Must specify scaling_A0")
    outfile += "    {:5d}   1  {:d} {:f} \n".format(len(TBMEs[:,0]),scaling_A0, scaling_p)
  else:
    outfile += "     {:5d}   {:3d}\n".format(len(TBMEs[:,0]), 0) 
  for i in range(len(TBMEs[:,0])):
    outfile += "{:3d} {:3d} {:3d} {:3d}  {:3d}  {:16.8f}\n".format(int(TBMEs[i,0]),int(TBMEs[i,1]),int(TBMEs[i,2]),int(TBMEs[i,3]),int(TBMEs[i,4]),TBMEs[i,5])


  with open(filename, 'w') as f:
    f.write(outfile)

  return True




def spider(inputfiles, names, type="M1", threshold=0.1, spinwindow=[], Eg_low=0, Eg_high=1e9, Ex_low=0, Ex_high=1e9, scale=2):

  Nsp = np.ceil(np.sqrt(len(inputfiles))).astype(int)
  f, ax_list = plt.subplots(Nsp,Nsp,squeeze=False,sharex='col', sharey='row')

  for i in range(len(inputfiles)):
    inputfile = inputfiles[i]
    name = names[i]
    ax = ax_list[i%Nsp][int((i-i%Nsp)/Nsp)]


    levels = read_energy_levels(inputfile)
    Egs = levels[0,0]

    Ex_high = min(Ex_high, levels[:,0].max()-Egs)
    Eg_high = min(Eg_high, levels[:,0].max()-Egs)
  
  
    levels_plot_J = []
    levels_plot_Ex = []
    for iEx in range(len(levels[:,0])):
      # print levels[iEx,:]
      J2 = levels[iEx,1]
      par = levels[iEx,2]
      if len(spinwindow) > 0 and not ([J2,par] in spinwindow or [J2-2,par] in spinwindow or [J2+2,par] in spinwindow):
        continue
      Ex = levels[iEx,0]-Egs
      if Ex < Ex_low or Ex > Ex_high:
        continue
  
      levels_plot_J.append(J2/2)
      levels_plot_Ex.append(Ex)
  
    ax.plot(levels_plot_J, levels_plot_Ex, 'o', color='grey', linewidth=0.5)
    ax.set_xlim([levels[:,1].min()/2-1,levels[:,1].max()/2+1])
    ax.set_title(name+r'$\,E_\gamma\in[{:.1f},{:.1f}]$'.format(Eg_low,Eg_high))
    ax.set_ylabel(r'$E_x\,\mathrm{[MeV]}$')
    ax.set_xlabel(r'$J\,\mathrm{[\hbar]}$')
  
  
  
    transitions = read_transition_strengths(inputfile, type=type)
    for iEx in range(len(transitions[:,0])):
      J2i = int(transitions[iEx,0])
      pari = int(transitions[iEx,1])
      if len(spinwindow) > 0 and not [J2i,pari] in spinwindow:
        continue
      B = transitions[iEx,7]
      Eg = transitions[iEx,6]
      if B < threshold or Eg<Eg_low or Eg>Eg_high:
        continue
      Ei = transitions[iEx,2]-Egs
      if Ei < Ex_low or Ei > Ex_high:
        continue
      J2f = int(transitions[iEx,3])
      parf = int(transitions[iEx,4])
      Ef = transitions[iEx,5]-Egs
      ax.plot([J2i/2,J2f/2],[Ei,Ef], color='teal', linewidth=(scale*B))
  
  return f, ax_list












