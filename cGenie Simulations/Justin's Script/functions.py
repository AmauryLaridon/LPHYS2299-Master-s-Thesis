###############################################################################
# THIS CODE CONTAINS ALL THE FUNCTIONS USED TO TRANSFORM THE DATA             #
###############################################################################

import numpy as np
import numpy.ma as ma
import scipy.integrate as integrate
import open_files as fi
from math import factorial

#-----------------------------------------------------------------------------#
# GLOBAL FUNCTIONS                                                            #
#-----------------------------------------------------------------------------#

# Return the min and max from a vector (any size)
def Minmax(x): 
    mini = np.nanmin(x)
    maxi = np.nanmax(x)
    return mini,maxi

# Return the annual mean from monthly data 
def Year_average(x):
    number_month = len(x)
    number_year = max(1, int(number_month/12))
    year_av = ma.copy(x[0:number_year,...])
    for i in range(number_year):
        frac_x = x[i*12:i*12+12,...]
        year_av[i] = ma.mean(frac_x,axis=0)
    return year_av

# Return a vector with the max of every time-slice of a 3D vector
def Year_max_3(x, x1_min, x1_max, x2_min, x2_max): 
    year_max = []
    for i in range(len(x)):
        year_max.append(np.nanmax(x[i,x1_min:x1_max+1,x2_min:x2_max+1]))
    return year_max

# Return a vector with the mean of every time-slice of a 3D vector
def Year_mean_3(x, x1_min, x1_max, x2_min, x2_max):
    year_mean = np.zeros(len(x))
    for i in range(len(x)):
        year_mean[i] = np.nanmean(x[i,x1_min:x1_max+1,x2_min:x2_max+1])
    return year_mean

# Return a vector with the max of every time-slice of a 4D vector
def Year_max_4(x, x1_min, x1_max, x2_min, x2_max, x3_min, x3_max):
    year_max = []
    for i in range(len(x)):
        year_max.append(np.nanmax(x[i,x1_min:x1_max+1,x2_min:x2_max+1,
                                    x3_min:x3_max+1]))
    return year_max

# Return a vector with the mean of every time-slice of a 4D vector
def Year_mean_4(x, x1_min, x1_max, x2_min, x2_max, x3_min, x3_max):
    year_mean = []
    for i in range(len(x)):
        year_mean.append(np.nanmean(x[i,x1_min:x1_max+1,x2_min:x2_max+1,
                                      x3_min:x3_max+1]))
    return year_mean

# Transform a salt forcing (PSU/kg/yr) into a freshwater forcing (in Sv)
def Salt_to_Sv_forcing(x):
    PSU = 34.9 # mean global salinity of the ocean
    base_flux = -PSU*1000 # the flux corresponding to 1 m³/yr of freshwater
    forcing_Sv = x/(base_flux*1e6*31557600) # to have forcing in Sv (10⁶m³/s)
    return forcing_Sv

# The invert function of Salt_to_Sv_forcing
def Sv_to_Salt_forcing(x):
    PSU = 34.9 # mean global salinity of the ocean
    base_flux = -PSU*1000 # the flux corresponding to 1 m³/yr of freshwater
    forcing_Salt = x*(base_flux*1e6*31557600) # to have forcing in Sv (10⁶m³/s)
    return forcing_Salt

print(Salt_to_Sv_forcing(-0.68319e+17))

# Function for smoothing other functions
def Savitzky_golay(y, window_size, order, deriv=0, rate=1):
     try:
         window_size = np.abs(np.int(window_size))
         order = np.abs(np.int(order))
     except ValueError:
         raise ValueError("window_size and order have to be of type int")
     if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
     if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
     order_range = range(order+1)
     half_window = (window_size -1) // 2
     # precompute coefficients
     b = np.mat([[k**i for i in order_range] for k in range(-half_window, 
                                                            half_window+1)])
     m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
     # pad the signal at the extremes with
     # values taken from the signal itself
     firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
     lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
     y = np.concatenate((firstvals, y, lastvals))
     return np.convolve( m[::-1], y, mode='valid')

#-----------------------------------------------------------------------------#
# PRESENT RELATIVE FUNCTIONS                                                  #
#-----------------------------------------------------------------------------#

# 2-dimensional mask for the Atlantic ocean
def Atlantic_2D_mask():
    atlantic_mask = np.ones((36,36))
    for i in range(9):
        for j in range(10):
            if ma.is_masked(fi.grid_mask[2+i,19+j]) == False:
                atlantic_mask[2+i,19+j] = 0
    for i in range(25):
        for j in range(11):
            if ma.is_masked(fi.grid_mask[11+i,16+j]) == False:
                atlantic_mask[11+i,16+j] = 0
            if j == 0 and i < 12:
                atlantic_mask[11+i,16+j] = 1
            if j == 1 and i < 11:
                atlantic_mask[11+i,16+j] = 1
            if j == 2 and i < 2:
                atlantic_mask[11+i,16+j] = 1
    atlantic_mask[28,25] = atlantic_mask[28,26] = atlantic_mask[29,26] = 1
    atlantic_mask[34,27] = atlantic_mask[35,27] = atlantic_mask[35,28] = 0
    atlantic_mask[35,29] = 0
    return ma.masked_array(atlantic_mask,atlantic_mask)+1

# 2-dimensional mask for the Atlantic ocean above a give latitude number
def Atlantic_2D_mask_sup_lat(lat_number):
    atlantic_mask_2D_lat = np.ones((36,36))
    atlantic_mask_2D = Atlantic_2D_mask()
    for i in range(36):
        for j in range(36):
            if ma.is_masked(atlantic_mask_2D[i,j]) == False and lat_number<=i<34:
                atlantic_mask_2D_lat[i,j] = 0
    return ma.masked_array(atlantic_mask_2D_lat,atlantic_mask_2D_lat)+1

aaa_2 = Atlantic_2D_mask_sup_lat(27)

# 3-dimensional mask for the Atlantic ocean
def Atlantic_3D_mask():
    atlantic_mask_3d = np.ones((16,36,36))
    atlantic_mask_2d = Atlantic_2D_mask()
    for z in range(16):
        for i in range(36):
            for j in range(36):        
                if (ma.is_masked(fi.grid_mask_3d[z,i,j]) 
                    == False) and atlantic_mask_2d[i,j] == 1:
                    atlantic_mask_3d[z,i,j] = 0
    return ma.masked_array(atlantic_mask_3d,atlantic_mask_3d)+1

# Return the longitude mean of 4 dimension vector in the Atlantic ocean
def Atlantic_long_mean(x,t):
    atlantic_mean = np.zeros((16,36))
    atlantic_mask_3d = Atlantic_3D_mask()
    masked_x = x[t,:,:,:]*atlantic_mask_3d[:,:,:]
    for i in range(np.size(atlantic_mean[0,:])):
        for z in range(np.size(atlantic_mean[:,0])):
            atlantic_mean[z,i] = ma.mean(masked_x[z,i,:])
    return atlantic_mean

# Return the "relative importance" of ocean layers (number of cells + height)
def Zt_vert_coef():
    zt_vert_coef = np.zeros(16)
    zt_height_coef = np.zeros(16)
    zt_nombre = np.zeros(16)
    for i in range(16):
        zt_height_coef[i] = (fi.zt_edges[i+1]-fi.zt_edges[i])/(fi.zt_edges[1]-
                                                               fi.zt_edges[0])
        zt_nombre[i] = ma.sum(fi.grid_mask_3d[i,:,:])
        zt_vert_coef[i] = zt_height_coef[i]*zt_nombre[i]
    return zt_vert_coef

# Return the "relative importance" of all the Atlantic ocean layers
def Atlantic_zt_vert_coef():
    at_zt_vert_coef = np.zeros(16)
    zt_height_coef = np.zeros(16)
    zt_nombre = np.zeros(16)
    atlantic_3D_mask = Atlantic_3D_mask()
    for i in range(16):
        zt_height_coef[i] = (fi.zt_edges[i+1]-fi.zt_edges[i])/(fi.zt_edges[1]-
                                                               fi.zt_edges[0])
        zt_nombre[i] = ma.sum(atlantic_3D_mask[i,:,:])
        at_zt_vert_coef[i] = zt_height_coef[i]*zt_nombre[i]
    return at_zt_vert_coef

#
# GEOSTROPHIC STREAMFUNCTION (Levang and Schmitt)
#

def Atlantic_edges(lat_number):
    atlantic_3D_mask = Atlantic_3D_mask()[:,lat_number,:]
    atlantic_edges = np.ones((16,36))
    for i in range(16):
        for j in range(35):
            if (ma.is_masked(atlantic_3D_mask[i,j-1]) == True and 
                ma.is_masked(atlantic_3D_mask[i,j]) == False and
               ma.is_masked(atlantic_3D_mask[i,j+1]) == False):
                atlantic_edges[i,j] = 0
            if (ma.is_masked(atlantic_3D_mask[i,j-1]) == False and 
                ma.is_masked(atlantic_3D_mask[i,j]) == False and
               ma.is_masked(atlantic_3D_mask[i,j+1]) == True):
                atlantic_edges[i,j] = 0
    return ma.masked_array(atlantic_edges,atlantic_edges)+1 

def Useful_quantities(lat_number):
    atlantic_3D_mask = Atlantic_3D_mask() 
    H = np.max((fi.grid_topo*Atlantic_2D_mask())[lat_number]) # max depth 
    H_number = 0 # max depth indice 
    for i in range(17):
        if int(fi.zt_moc[i]+0.5) == int(H+0.5):
            H_number = i
            break
    number_box = np.zeros(17) # number of box per depth layer
    for i in range(H_number):
        number_box[i] = np.sum(atlantic_3D_mask[i,lat_number])
    return H, H_number, number_box

def Geo_stream_all_layers(lat_number, rho_field, atlantic_edges, H, H_number, 
                          number_box):
    # Some needed quantities
    lat_deg = fi.lat[lat_number] # Latitude in degree
    g = 9.81 # acceleration constant
    rho_0 = 1027 # reference density in kg/m³
    earth_radius = 6.371e6
    # length of a grid cell (from east to west) : depends on the latitude !
    length = 2*np.pi*earth_radius*np.cos(lat_deg*np.pi/180)/36 
    f_coriolis = 2*7.2921e-5*np.sin(lat_deg*np.pi/180) # Coriolis coef
    alpha = g/(rho_0*f_coriolis)
    # To account for under water mountain
    montagne = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                0,0,0,0,[2,2],[2,2],0,[3,2],0,0,0,0,0,
                [[2,5],[2,5],[2,5]],[[4,3],[4,3],[4,3],[2,3],[2,3]]]
    
    def L():
        L = int(length)*np.ones(17) 
        for i in range(17):
            L[i] = L[i]*number_box[i]
        return np.flip(L)
    
    def Int_bound(z):
        int_bound = -fi.zt_edges
        for i in range(17):
            if -z <= -fi.zt_edges[i]:
                int_bound[i] = -z
            if int_bound[i] < -H:
                int_bound[i] = -H
        return np.flip(int_bound)
    
    def East_west_rho_anom():
        density_edge = atlantic_edges*rho_field[:,lat_number,:]
        difference = np.zeros(16)
        montain_it = 0
        for i in range(16):
            if ma.is_masked(ma.max(density_edge[i])-ma.min(density_edge[i])):
                pass
            else:
                indice = ma.nonzero(density_edge[i])[0]
                if len(indice) == 2:
                    difference[i] = (density_edge[i, indice[0]]-
                                     density_edge[i, indice[1]])
                else: # The mountain problem
                    for k in range(int(len(indice)/2)):
                        if lat_number < 34:
                            difference[i] += (number_box[i]*
                (density_edge[i, indice[2*k]]- density_edge[i, indice[2*k+1]])/
                montagne[lat_number][k])
                        else:
                            difference[i] += (number_box[i]*
                (density_edge[i, indice[2*k]]- density_edge[i, indice[2*k+1]])/
                montagne[lat_number][montain_it][k])
                    montain_it += 1
        return np.flip(difference)

    rho_anomaly = East_west_rho_anom()

    def Zonally_average_velo(z):
        integral = 0
        bound = Int_bound(z)
        L_z = L()
        for i in range(16):
            if L_z[i+1] == 0:
                integral = integral
            else:
                integral += integrate.quad(lambda z: 1/L_z[i+1]*rho_anomaly[i],
                                           bound[i], bound[i+1])[0]
        return alpha*integral
    
    zav_vector = np.zeros(16) 
    for i in range(16):
        zav_vector[i] = Zonally_average_velo(fi.zt_moc[i])
    zav_vector = np.flip(zav_vector)

    def Int_zav(z):
        integral = 0
        bound = Int_bound(0) 
        for i in range(16):
            integral += integrate.quad(lambda z: zav_vector[i], 
                                       bound[i], bound[i+1])[0]
        return integral
    
    int_zav = np.zeros(16)
    for i in range(16):
        int_zav[i] = Int_zav(fi.zt_moc[i])
    zav_min_int_zav = zav_vector - 1/H*int_zav
    
    def Streamfunction():
        integral = np.zeros(17)
        bound = Int_bound(0)
        for i in range(16):
            integral[i+1] = integral[i] + (integrate.quad(lambda z: 
                            zav_min_int_zav[i], bound[i], bound[i+1])[0])
        for i in range(17):
            integral[i] = integral[i]*-length*number_box[-1-i]/1000000
        return np.flip(integral)
    
    return Streamfunction(), np.flip(zav_vector), np.flip(rho_anomaly)

# The Total streamfunction with cGENIE velocity field
def Total_streamfunction(v_field, number_box, dimension=3):
    earth_radius = 6.371e6
    length = 2*np.pi*earth_radius*np.cos(fi.lat[31]*np.pi/180)/36
    total_stream = np.zeros(17)
    if dimension == 3:
        atlantic_velo = (Atlantic_3D_mask()*v_field)[:,31,:]
        for i in range(16):
            if np.ma.is_masked(np.sum(atlantic_velo[-1-i,:])):
                total_stream[i] = np.nan
                a = i
            else:
                total_stream[i+1] = total_stream[i] + np.sum(atlantic_velo[-1-i
                                 ,:])*length*-(fi.zt_moc[-1-i]-fi.zt_moc[-2-i])
    else:
        atlantic_velo = v_field
        a = 0
        for i in range(16):
            if atlantic_velo[-1-i]==0:
                total_stream[i] = np.nan
                a = i
            else:
                total_stream[i+1] = total_stream[i] + atlantic_velo[-1-i
                   ]*number_box[-1-i]*length*-(fi.zt_moc[-1-i]-fi.zt_moc[-2-i])
    total_stream[-1] = total_stream[a+1] = np.nan
    return np.flip(total_stream/1000000)

# The UNESCO formula used to compute the density in cGENIE
def GOLDSTEIN_dens_eq(T, S):
    density = 1000.0 + (0.7968*S - 0.0559*T - 0.0063*T**2 + 3.7315e-05*T**3)
    return density

# Function that computes the max of a given streamfunction
def Max_streamfunction(lat_number, opsia_field, opsia_field_ini, rho_field,
                       temp_field, temp_field_ini, sal_field, sal_field_ini):
    
    H, H_number, number_box = Useful_quantities(lat_number)

    atlantic_edges = Atlantic_edges(lat_number)
    max_opsia = np.zeros(len(rho_field)+1)
    true_max = np.zeros(len(rho_field)+1)
    for i in range(len(rho_field)):
        max_opsia[i+1]=np.nanmax(Geo_stream_all_layers(lat_number,rho_field[i],
                                 atlantic_edges, H, H_number, number_box)[0])
        true_max[i+1] = ((np.nanmax(opsia_field[i,5:H_number-1,lat_number]*np.sin(np.pi/
                        180*fi.lat_moc[lat_number])) + np.nanmax(opsia_field[i,
                        5:H_number-1,lat_number+1]*np.sin(np.pi/180*fi.lat_moc[lat_number
                        +1])))/(2*np.sin(np.pi/180*fi.lat[lat_number])))

    dens_only_T = []
    dens_only_S = []
    for i in range(len(rho_field)):
        dens_only_T.append(GOLDSTEIN_dens_eq(temp_field[i],
                                                 sal_field_ini))
        dens_only_S.append(GOLDSTEIN_dens_eq(temp_field_ini,
                                                 sal_field[i]))

    max_opsia_T_only = np.zeros(len(rho_field)+1)
    max_opsia_S_only = np.zeros(len(rho_field)+1)
    for i in range(len(rho_field)):
        max_opsia_T_only[i+1] = np.nanmax(Geo_stream_all_layers(lat_number, 
                   dens_only_T[i], atlantic_edges, H, H_number, number_box)[0][5:H_number-1])
        max_opsia_S_only[i+1] = np.nanmax(Geo_stream_all_layers(lat_number, 
                   dens_only_S[i], atlantic_edges, H, H_number, number_box)[0][5:H_number-1])

    true_max[0] = (
    (np.nanmax(opsia_field_ini[:,lat_number]*np.sin(np.pi/180*fi.lat_moc[lat_number
    ])) + np.nanmax(opsia_field_ini[:,lat_number+1]*np.sin(np.pi/180*fi.lat_moc[
    lat_number+1])))/(2*np.sin(np.pi/180*fi.lat[lat_number])))
    
    rho_field_ini = GOLDSTEIN_dens_eq(temp_field_ini, sal_field_ini)
    
    max_opsia_T_only[0] = max_opsia_S_only[0] = max_opsia[0] = np.nanmax(
    Geo_stream_all_layers(lat_number,rho_field_ini,atlantic_edges, H, H_number, 
                          number_box)[0])
    return max_opsia, max_opsia_T_only, max_opsia_S_only, true_max

def Max_streamfunction_lat(lat_number, opsia_field, opsia_field_ini, rho_field,
                       temp_field, temp_field_ini, sal_field, sal_field_ini):
    
    H = np.zeros(len(lat_number))
    H_number = np.zeros(len(lat_number))
    number_box = np.zeros((len(lat_number),17))
    for i in range(len(rho_field)):
        if i == 0:
            H[i], H_number[i], number_box[i,:] = Useful_quantities(int(lat_number[i]))
        elif i == i-1:
            H[i], H_number[i], number_box[i,:] = H[i-1], H_number[i-1], number_box[i-1,:]
        else:
            H[i], H_number[i], number_box[i,:] = Useful_quantities(int(lat_number[i]))
        
    max_opsia = np.zeros(len(rho_field)+1)
    true_max = np.zeros(len(rho_field)+1)
    for i in range(len(rho_field)):
        max_opsia[i+1]= np.nanmax(Geo_stream_all_layers(int(lat_number[i]),rho_field[i],
                                 Atlantic_edges(int(lat_number[i])), H[i], H_number[i], number_box[i])[0])
        true_max[i+1] = ((np.nanmax(opsia_field[i,:,int(lat_number[i])][5:int(H_number[i]-1)]*np.sin(np.pi/
                        180*fi.lat_moc[int(lat_number[i])])) + np.nanmax(opsia_field[i,
                        :,int(lat_number[i])+1][5:int(H_number[i]-1)]*np.sin(np.pi/180*fi.lat_moc[int(lat_number[i])
                        +1])))/(2*np.sin(np.pi/180*fi.lat[int(lat_number[i])])))

    dens_only_T = []
    dens_only_S = []
    for i in range(len(rho_field)):
        dens_only_T.append(GOLDSTEIN_dens_eq(temp_field[i],
                                                 sal_field_ini))
        dens_only_S.append(GOLDSTEIN_dens_eq(temp_field_ini,
                                                 sal_field[i]))

    max_opsia_T_only = np.zeros(len(rho_field)+1)
    max_opsia_S_only = np.zeros(len(rho_field)+1)
    for i in range(len(rho_field)):
        max_opsia_T_only[i+1] = np.nanmax(Geo_stream_all_layers(int(lat_number[i]), 
                   dens_only_T[i], Atlantic_edges(int(lat_number[i])), H[i], H_number[i], number_box[i])[0][5:int(H_number[i]-1)])
        max_opsia_S_only[i+1] = np.nanmax(Geo_stream_all_layers(int(lat_number[i]), 
                   dens_only_S[i], Atlantic_edges(int(lat_number[i])), H[i], H_number[i], number_box[i])[0][5:int(H_number[i]-1)])

    true_max[0] = (
    (np.nanmax(opsia_field_ini[:,int(lat_number[0])]*np.sin(np.pi/180*fi.lat_moc[int(lat_number[0])
    ])) + np.nanmax(opsia_field_ini[:,int(lat_number[0])+1]*np.sin(np.pi/180*fi.lat_moc[
    int(lat_number[0])+1])))/(2*np.sin(np.pi/180*fi.lat[int(lat_number[0])])))
        
    rho_field_ini = GOLDSTEIN_dens_eq(temp_field_ini, sal_field_ini)
    
    max_opsia_T_only[0] = max_opsia_S_only[0] = max_opsia[0] = np.nanmax(
        Geo_stream_all_layers(int(lat_number[0]),rho_field_ini,Atlantic_edges(int(lat_number[0])), 
                              H[0], H_number[0], number_box[0])[0])
        
    return max_opsia, max_opsia_T_only, max_opsia_S_only, true_max

# Function for computing the northward advection of salt vs the atmospheric haline contribution
def Advective_vs_atmospheric(lat_number, time_number, v_field, sal_field, precip_field, 
                             evap_field):
    earth_radius = 6.371e6
    lat_deg = fi.lat[lat_number] # Latitude in degree
    
    # length of a grid cell (from east to west) : depends on the latitude !
    length = 2*np.pi*earth_radius*np.cos(lat_deg*np.pi/180)/36
    salt_export = 0 # The amount of salt export (positive northward) in kg/s
    salt_atmos = 0 # The amount of salt from atmosphere (positive 'in') in kg/s
    
    """ We suppose that the output of cGENIE is already in m/s
    # Transorm the P and E field in m/s
    evap_field, precip_field = evap_field*1e-3, precip_field*1e-3
    """
    H, H_number, number_box = Useful_quantities(lat_number)
    
    def L():
        L = int(length)*np.ones(17) 
        for i in range(17):
            L[i] = L[i]*number_box[i]
        return L
    
    L = L()
    
    # The salinity is in PSU = 1g/kg = 1kg/1000kg = 1kg/m³ (freshwater)
    # The computation of the northward transport of salt (in kg/s)
    v_times_salinity = Atlantic_3D_mask()[:,lat_number,:]*((fi.ocn_v_RCP85_100yr_from_2000*
                       fi.ocn_sal_RCP85_100yr_from_2000)[time_number,:,lat_number,:])
    for i in range(16):
        salt_export += integrate.quad(lambda z: L[i]*np.mean(v_times_salinity[i,:]), 
                                   fi.zt_edges[i], fi.zt_edges[i+1])[0]
        
    # The computation of the atmospheric freshwater forcing (salt flux in kg/s)
    atlantic_2D_mask_sup_lat = Atlantic_2D_mask_sup_lat(lat_number)
    P_minus_E = (precip_field-evap_field)[time_number]*atlantic_2D_mask_sup_lat
    
    salt_atmos = -0.035*fi.grid_area[35,35]*np.sum(P_minus_E)
    
    return salt_export, salt_atmos  # Je triche avec le 1e6

"""
Advective_vs_atmospheric(27, 0, fi.ocn_v_RCP85_1000yr_from_2000, 
           fi.ocn_sal_RCP85_1000yr_from_2000, fi.precip_RCP85_1000yr_from_2000, 
           fi.evap_RCP85_1000yr_from_2000)

for i in range(20):
    print(Advective_vs_atmospheric(27, i, fi.ocn_v_RCP85_1000yr_from_2000, 
               fi.ocn_sal_RCP85_1000yr_from_2000, fi.precip_RCP85_1000yr_from_2000, 
               fi.evap_RCP85_1000yr_from_2000))

aaa = (fi.precip_RCP85_1000yr_from_2000-fi.evap_RCP85_1000yr_from_2000)[0]*Atlantic_2D_mask_sup_lat(27)
aaa_2 = (fi.precip_RCP85_1000yr_from_2000-fi.evap_RCP85_1000yr_from_2000)[0]*1e-3*Atlantic_2D_mask_sup_lat(27)
"""
#-----------------------------------------------------------------------------#
# DEVONIAN RELATIVE FUNCTIONS                                                 #
#-----------------------------------------------------------------------------#