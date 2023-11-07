################################################################################
# OPEN FILES CODE                                                              #
################################################################################

from netCDF4 import Dataset
import time
start_time = time.time()

#------------------------------------------------------------------------------#
# Paths                                                                        #
#------------------------------------------------------------------------------#

cofast = "/cofast/jgerard/cgenie_output/"
paleofast = "/cofast/jgerard/cgenie.muffin/genie-paleo/"
biogem_2d = "/biogem/fields_biogem_2d.nc"
biogem_3d = "/biogem/fields_biogem_3d.nc"

Vera_time = ["370", "383", "393", "408", "420"]

Scot_time = ["360","365","370","375","380","385", "390","395",
             "400","405","410","415","420"]

simulations = []

for t in Scot_time:
    simulations.append("Scot_main_config/JG.Scot"+t+"M.main_config")
    simulations.append("Scot_main_config/JG.Scot"+t+"M.main_config.sol_cst")
simulations.append("Scot_main_config/JG.Scot390M.main_config_20000")
simulations.append("Scot_main_config/JG.Scot390M.main_config.highsol")
simulations.append("Scot_main_config/JG.Scot390M.main_config.low_sol")

for t in Vera_time:
    simulations.append("Vera_main_config/JG.Vera"+t+"M.main_config")
simulations.append("Vera_main_config/JG.Vera393_2M.main_config")

def Open_files(location, exp_name, dim_bio, var):
    return Dataset(location+exp_name+dim_bio).variables[var][:]

def Open_biogem_series(location, exp_name, var):
    return open(location+exp_name+"/biogem/biogem_series_"+var+".res", "r")

def Open_paleo(exp_name, var):
    return open(paleofast+exp_name+var, "r")

def Print_info(x, style):
    size = len(x)
    if size%2 == 0:
        dash = int((80-size)/2)*style
        print(dash+x+dash)
    else:
        dash = int((80-size-1)/2)*style
        dash_2 = int((80-size+1)/2)*style
        print(dash+x+dash_2)
    return None

#------------------------------------------------------------------------------#
# 3D variables                                                                 #
#------------------------------------------------------------------------------#

#
# Temperature (°C)
#

Print_info("Start loading temperature", "-")
ocn_temp = {}
for simu in simulations:
    ocn_temp[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_3d,
                                     'ocn_temp')
Print_info("End loading temperature", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Oxygen (mol/kg)
#

Print_info("Start loading dissolved oxygen", "-")
ocn_O2 = {}
for simu in simulations:
    ocn_O2[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_3d,
                                   'ocn_O2')
Print_info("End loading dissolved oxygen", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Dissolved inorganic carbon 13 (/1000)
#

Print_info("Start loading dissolved inorganic carbon 13", "-")
ocn_DIC_13C = {}
for simu in simulations:
    ocn_DIC_13C[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_3d,
                                        'ocn_DIC_13C')
Print_info("End loading dissolved inorganic carbon 13", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Ventilation age (yr)
#

Print_info("Start loading ventilation age", "-")
ocn_vent = {}
for simu in simulations:
    ocn_vent[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_3d,
                                      'misc_col_Dage')
Print_info("End loading ventilation age", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#------------------------------------------------------------------------------#
# 2D variables                                                                 #
#------------------------------------------------------------------------------#

#
# Benthic temperature (°C)
#

Print_info("Start loading benthic temperature", "-")
ocn_ben_temp = {}
for simu in simulations:
    ocn_ben_temp[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                         'ocn_ben_temp')
Print_info("End loading benthic temperature", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Benthic dissolved oxygen (mol/kg)
#

Print_info("Start loading benthic dissolved oxygen", "-")
ocn_ben_O2 = {}
for simu in simulations:
    ocn_ben_O2[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                       'ocn_ben_O2')
Print_info("End loading benthic dissolved oxygen", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Benthic dissolved inorganic carbon 13 (/1000)
#

Print_info("Start loading benthic dissolved inorganic carbon 13", "-")
ocn_ben_DIC_13C = {}
for simu in simulations:
    ocn_ben_DIC_13C[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                            'ocn_ben_DIC_13C')
Print_info("End loading benthic dissolved inorganic carbon 13", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Solar forcing (W/m²)
#

Print_info("Start loading solar forcing", "-")
phys_solfor = {}
for simu in simulations:
    phys_solfor[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                        'phys_solfor')
Print_info("End loading solar forcing", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Streamfunction (Sv)
#

Print_info("Start loading overturning streamfunction", "-")
phys_opsi = {}
for simu in simulations:
    phys_opsi[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                      'phys_opsi')
Print_info("End loading overturning streamfunction", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")


Print_info("Start loading barotropic streamfunction", "-")
phys_psi = {}
for simu in simulations:
    phys_psi[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                     'phys_psi')
Print_info("End loading barotropic streamfunction", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Convective cost (column integrated adjustments per year)
#

Print_info("Start loading convective cost", "-")
phys_cost = {}
for simu in simulations:
    phys_cost[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                      'phys_cost')
Print_info("End loading convective cost", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Zonally averaged wind-stress (N/m^2)
#

Print_info("Start loading zonal wind-stress", "-")
phys_tau_u = {}
for simu in simulations:
    phys_tau_u[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                       'phys_tau_u')
Print_info("End loading zonal wind-stress", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#
# Meridional averaged wind-stress (N/m^2) ==> ALWAYS 0 IN cGENIE !!!
#

Print_info("Start loading meridional wind-stress", "-")
phys_tau_v = {}
for simu in simulations:
    phys_tau_v[simu[20:]] = Open_files(cofast,simu+".SPIN",biogem_2d,
                                       'phys_tau_v')
Print_info("End loading meridional wind-stress", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#------------------------------------------------------------------------------#
# 1D variables                                                                 #
#------------------------------------------------------------------------------#

Print_info("Start loading 1D variables", "-")
var = simulations[0]+".SPIN"
lat = Open_files(cofast,var,biogem_3d,"lat")
lat_edges = Open_files(cofast,var,biogem_3d,"lat_edges")
lat_moc = Open_files(cofast,var,biogem_2d,"lat_moc")
lat_psi = Open_files(cofast,var,biogem_2d,"lat_psi")
lon = Open_files(cofast,var,biogem_3d,"lon")
lon_edges = Open_files(cofast,var,biogem_3d,"lon_edges")
lon_psi = Open_files(cofast,var,biogem_2d,"lon_psi")
zt = Open_files(cofast,var,biogem_3d,"zt")
zt_edges = Open_files(cofast,var,biogem_3d,"zt_edges")
zt_moc = Open_files(cofast,var,biogem_2d,"zt_moc")
time_year = Open_files(cofast,var,biogem_3d,"time")
var = "Scot_main_config/JG.Scot390M.main_config_20000.SPIN"
time_year_2 = Open_files(cofast,var,biogem_3d,"time")
albedo = {}
for simu in simulations:
    if simu[-11:] != "main_config":
        pass 
    elif simu[20:28] == "Vera393_":
        f = Open_paleo("Vera3932/","Vera3932"+".albd.dat")
        albedo[simu[20:28]]= f.read().splitlines()
        f.close()
    else:
        f = Open_paleo(simu[20:28]+"/",simu[20:28]+".albd.dat")
        albedo[simu[20:28]]= f.read().splitlines()
        f.close()
Print_info("End loading 1D variables", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#------------------------------------------------------------------------------#
# Mean values (cGENIE output)                                                  #
#------------------------------------------------------------------------------#

Print_info("Start loading mean values data", "-")
mean_temp_TS = {}
for simu in simulations:
    f = Open_biogem_series(cofast,simu+".SPIN","ocn_temp")
    mean_temp_TS[simu[20:]] = f.readlines()
    f.close()
    
mean_sal_TS = {}
for simu in simulations:
    f = Open_biogem_series(cofast,simu+".SPIN","ocn_sal")
    mean_sal_TS[simu[20:]] = f.readlines()
    f.close()   

mean_O2_TS = {}
for simu in simulations:
    f = Open_biogem_series(cofast,simu+".SPIN","ocn_O2")
    mean_O2_TS[simu[20:]] = f.readlines()
    f.close()   
    
mean_vent_TS = {}
for simu in simulations:
    f = Open_biogem_series(cofast,simu+".SPIN","misc_col_age")
    mean_vent_TS[simu[20:]] = f.readlines()
    f.close()   

mean_DIC_13C_TS = {}
for simu in simulations:
    f = Open_biogem_series(cofast,simu+".SPIN","ocn_DIC_13C")
    mean_DIC_13C_TS[simu[20:]] = f.readlines()
    f.close()         
Print_info("End loading mean values data", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#------------------------------------------------------------------------------#
# Maks                                                                         #
#------------------------------------------------------------------------------#

Print_info("Start loading masks", "-")
grid_mask = {}
grid_mask_3d = {}
grid_topo = {}
grid_area = {}
for simu in simulations:
    if simu[-11:] != "main_config":
        pass
    else:
        grid_mask[simu[20:28]]=Open_files(cofast,simu+".SPIN",biogem_3d,
                                          'grid_mask')
        grid_mask_3d[simu[20:28]]=Open_files(cofast,simu+".SPIN",biogem_3d,
                                             'grid_mask_3d')
        grid_topo[simu[20:28]]=Open_files(cofast,simu+".SPIN",biogem_3d,
                                          'grid_topo')
        grid_area[simu[20:28]]=Open_files(cofast,simu+".SPIN",biogem_3d,
                                          'grid_area')
Print_info("End loading masks", "-")
Print_info("Run time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")

#------------------------------------------------------------------------------#
# All is fine                                                                  #
#------------------------------------------------------------------------------#

print("")
Print_info("", "#")
Print_info(" ALL FILES LOAD CORRECTLY ", "#")
Print_info("", "#")
Print_info("Total time : "+"{:.2f}".format(time.time() - start_time)+"s", " ")
