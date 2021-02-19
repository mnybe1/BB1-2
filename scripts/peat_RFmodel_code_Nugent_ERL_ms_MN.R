# R script to calculate Radiative forcing from continuous, variable
# emissions of GHGs.  
# Follows Frolking and Roulet 2007 Global Change Biology, and Frolking et al. 2006 J. Geophysical Research.
# uses parameterization from Dommain et al.  2018 Global Change Biology, based on Joos et al. Atmos. Chem. Phys. 13, 2793-2825, but with finite lifetime slow pool

model_version = "v15"

#wd = getwd()

wd = "~/Google Drive/Micromet Lab/People/2019-Marion Nyberg/BB1-2"

input_dir = paste(wd,"/peat_RF_model_input_files/",sep="")
output_dir = paste(wd,"/peat_RF_model_output_files/",sep="")



# LOAD DATA FROM CSV FILE

interpolate_input = 0  # 0 if input is annual data (no interpolation), 1 otherwise

LU_name = "restored"
site_name = "BB1"

#infile_name = paste(input_dir,site_name,"_",LU_name,"_3",".csv", sep="")

infile_name = paste(input_dir,site_name,"_",LU_name,".csv", sep="")

#fluxdata_in <- read.table(infile_name, skip = 1, sep = ",", header = TRUE)

fluxdata_in <- read.table(infile_name, sep = ",", header = TRUE)

# flux input values are kg CO2/m2/y and kg CH4/m2/y and kg N2O/m2/y (all N2O fluxes are zero for now)

## VARIABLE AND PARAMETER ARRAYS for Atmospheric CO2 perturbation pools
tau_co2 = array(0,dim=c(5))
frac_co2 = array(0,dim=c(5))   
atmos_co2 = array(0,dim=c(5))

## PARAMETER VALUES  

# radiative efficiencies
rad_eff_ch4 = 0.363e-3     # W/m2/ppb
rad_eff_n2o = 3.00e-3      # W/m2/ppb
rad_eff_co2 = 0.0137e-3    # W/m2/ppb
# from IPCC AR4 7.4.1.1 Biogeochemistry and Budgets of Methane: 2.78 Tg(CH4) per ppb
#   so (2.78e12 / 16) moles CH4 /ppb
rad_eff_ch4 = rad_eff_ch4 / ((2.78e12 / 16) * 16 * 0.001)   # W/m2/ppb to W/m2/kg
#     divided by moles/ppb as  moles CH4/ppb * g CH4/mole * 0.001 kg / g
#rad_eff_n2o = rad_eff_n2o / ((2.78e12 / 16) * 44 * 0.001)   # W/m2/ppb to W/m2/kg
rad_eff_co2 = rad_eff_co2 / ((2.78e12 / 16) * 44 * 0.001)   # W/m2/ppb to W/m2/kg

tau_ch4 = 12.4    # atmos. pool lifetime (years)  from Myrhe et al. 2013 (IPCC AR5 WG1)
tau_n2o = 121.    # atmos. pool lifetime (years)  from Myrhe et al. 2013 (IPCC AR5 WG1) 

## METHANE INDIRECT EFFECT MULTIPLIER 
# factor for CH4 indirect effect on warming 
#     (±0.3; IPCC AR5 Ch. 8, Myhre et al. 2013, see §8.SM.11.3.2 in supplement)
indirect_eff_ch4 = 1.65   
indirect_eff_ch4_high = 1.95
indirect_eff_ch4_low  = 1.35

# two possible paramterizations for radiative forcing model
rad_force_model = 1 # = 1 or 2 
# if 1: model params from Joos et al. 2013; Atmos. Chem. Phys. 13, 2793-2825, 
#       mean model fit, but with finite lifetime slow pool
# if 2: model params from  # Dommain et al GBC 2018;  
#       this may be better for >1000yr simulations; it includes CaCO3 compensation in ocean

if (rad_force_model == 1) {
  rad_force_model_name = "Jea13"    
  tau_co2 = c(4.304, 36.54, 394.4, 2e5, -9999.99) # atmos. pool lifetimes [years] 
  #                               2e5 modified from Jea13 value of 'infinity'
  frac_co2 = c(0.2763, 0.2824, 0.2240, 0.2173, 0.0000) # partitioning of CO2 into atmospheric lifetime pools, pool #5 not used
  
} else if (rad_force_model == 2) {
  rad_force_model_name = "Dea18"     # Dommain et al GBC 2018 
  tau_co2 = c(4.304, 36.54, 394.4, 7000, 2e5) # atmos. pool lifetimes [years]
  #                               7000y pool added by Joos in pers. comm; 2e5 modified from Jea13 value of 'infinity'
  frac_co2 = c(0.2763, 0.2824, 0.2240, 0.1473, 0.0700) # partitioning of CO2 into atmospheric lifetime pools
  
}

## READ IN OR GENERATE ANNUAL FLUXES (POSITIVE TO ATMOSPHERE)

fluxdata_in_matrix = as.matrix(fluxdata_in)

flux_data_years = fluxdata_in_matrix[,1]  
flux_co2 = fluxdata_in_matrix[,3]       # kg CO2/m2/y
flux_ch4 = fluxdata_in_matrix[,4]      # kg CH4/m2/y
#flux_n2o = fluxdata_in_matrix[,4]      # kg N2O/m2/y

#first_year = min(flux_data_years) # CE
#last_year = max(flux_data_years) # CE

first_year = as.numeric(min(flux_data_years)) # CE
last_year = as.numeric(max(flux_data_years)) # CE

num_years = last_year - first_year + 1

years = array(0,dim=c(num_years,1))
years = first_year:last_year

# INTERPOLATE TO ANNUAL FLUXES If NECESSARY 
if (interpolate_input == 1) {
  flux_co2_ann = array(0,dim=c(num_years,1))
  flux_ch4_ann = array(0,dim=c(num_years,1))
  #flux_n2o_ann = array(0,dim=c(num_years,1))
  
  # LINEAR INTERPOLATION BETWEEN VALUES
  flux_co2_ann = approx(flux_data_years, flux_co2, n=num_years)$y 
  flux_ch4_ann = approx(flux_data_years, flux_ch4, n=num_years)$y
  #flux_n2o_ann = approx(flux_data_years, flux_n2o, n=num_years)$y
  
} else {
  flux_co2_ann = as.numeric(flux_co2)
  flux_ch4_ann = as.numeric(flux_ch4)
  #flux_n2o_ann = flux_n2o
  
}

# convert flux per m2 to per ha
flux_co2_ann = flux_co2_ann * 1e4  
flux_ch4_ann = flux_ch4_ann * 1e4  
#flux_n2o_ann = flux_n2o_ann * 1e4

# initialize variable arrays
atmos_co2 = array(0,dim=c(num_years,5))   
atmos_co2_tot = array(0,dim=c(num_years,1))  
atmos_ch4 = array(0,dim=c(num_years,1)) 

#atmos_n2o = array(0,dim=c(num_years,1)) 

rad_force_co2 = array(0,dim=c(num_years,1))    
rad_force_ch4 = array(0,dim=c(num_years,1))    
#rad_force_n2o = array(0,dim=c(num_years,1))  
rad_force_tot = array(0,dim=c(num_years,1)) 
time = array(0,dim=c(num_years,1)) 

## CALCULATE ATMOSPHERIC POOLS

# year 1: initialize perturbation pools with first year fluxes
for (j in 1:5)  {
  atmos_co2[1,j] = flux_co2_ann[1] * frac_co2[j]   
}
atmos_ch4[1] = flux_ch4_ann[1]  
atmos_ch4 = flux_ch4_ann  

#atmos_n2o[1] = flux_n2o_ann[1]   

# loop through all remaining years, removing loss fraction then adding input flux - NOT SURE ABOUT THIS BECAUSE ONLY HAVE 1 YEAR
for (iyear in  2:num_years) {
  
  atmos_ch4[iyear] = atmos_ch4[iyear-1] * (1. - 1/tau_ch4) + flux_ch4_ann[iyear]   
  #atmos_n2o[iyear] = atmos_n2o[iyear-1] * (1. - 1/tau_n2o) + flux_n2o_ann[iyear]   
  for (j in 1:5) {
    atmos_co2[iyear,j] = atmos_co2[iyear-1,j] * (1. - 1/tau_co2[j]) + flux_co2_ann[iyear] * frac_co2[j]  
  }
}

# aggregate atmosphere CO2 perturbation pools into a single total pool
atmos_co2_tot = atmos_co2[,1] + atmos_co2[,2] + atmos_co2[,3] + atmos_co2[,4] + atmos_co2[,5]

## CALCULATE RADIATIVE FORCING (directly proportional to magnitude of atmospheric perturbation pool)
# include methane indirect effect (mean, high, low)
rad_force_co2 = atmos_co2_tot * rad_eff_co2   
rad_force_ch4 = atmos_ch4 * rad_eff_ch4 * indirect_eff_ch4   
rad_force_ch4_high = atmos_ch4 * rad_eff_ch4 * indirect_eff_ch4_high
rad_force_ch4_low = atmos_ch4 * rad_eff_ch4 * indirect_eff_ch4_low
#rad_force_n2o = atmos_n2o * rad_eff_n2o   
rad_force_tot = rad_force_co2 + rad_force_ch4   
rad_force_tot_high = rad_force_co2 + rad_force_ch4_high   
rad_force_tot_low = rad_force_co2 + rad_force_ch4_low  


## WRITE ANNUAL RESULTS TO FILE, INCLUDING METHANE HIGH AND LOW FROM INDIRECT EFFECT RANGE (1.65±0.3)
#temp.out.array = array(0,c(num_years,20))
temp.out.array = array(0,c(num_years,17))


#colnames(temp.out.array) <- c("simulation_year","co2_flx", "ch4_flx", "n2o_flx", "co2_atm", 
 #                             "ch4_atm", "n2o_atm", "co2_atm_1", "co2_atm_2", "co2_atm_3", "co2_atm_4", 
  #                            "co2_atm_5", "rf_co2", "rf_ch4_high", "rf_ch4", "rf_ch4_low", "rf_n2o", 
   #                           "rf_tot_high","rf_tot","rf_tot_low")

colnames(temp.out.array) <- c("simulation_year","co2_flx", "ch4_flx", "co2_atm", 
                              "ch4_atm", "co2_atm_1", "co2_atm_2", "co2_atm_3", "co2_atm_4", 
                              "co2_atm_5", "rf_co2", "rf_ch4_high", "rf_ch4", "rf_ch4_low", 
                              "rf_tot_high","rf_tot","rf_tot_low")


temp.out.array[,1]= years
temp.out.array[,2]= flux_co2_ann
temp.out.array[,3]= flux_ch4_ann
#temp.out.array[,4]= flux_n2o_ann
temp.out.array[,4]= atmos_co2_tot
temp.out.array[,5]= atmos_ch4
#temp.out.array[,7]= atmos_n2o
temp.out.array[,6]= atmos_co2[,1]
temp.out.array[,7]= atmos_co2[,2]
temp.out.array[,8]= atmos_co2[,3]
temp.out.array[,9]= atmos_co2[,4]
temp.out.array[,10]= atmos_co2[,5]
temp.out.array[,11]= rad_force_co2
temp.out.array[,12]= rad_force_ch4_high
temp.out.array[,13]= rad_force_ch4
temp.out.array[,14]= rad_force_ch4_low
#temp.out.array[,15]= rad_force_n2o
temp.out.array[,15]= rad_force_tot_high
temp.out.array[,16]= rad_force_tot
temp.out.array[,17]= rad_force_tot_low

date = Sys.time()
date_str.temp = as.character(date)
date_str = chartr(":",".",date_str.temp)

outfile = paste(output_dir,site_name,"_",model_version, '_output_',rad_force_model_name,'_',date_str,'.csv',sep='')   
write.table(temp.out.array, file=outfile, col.names = TRUE, row.names = FALSE, sep=", ")





