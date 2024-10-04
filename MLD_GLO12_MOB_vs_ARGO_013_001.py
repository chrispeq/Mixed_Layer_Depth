## conda activate MLD_CMT
# export COPERNICUSMARINE_DISABLE_SSL_CONTEXT=True

## this assumes that observations from 
# ftp://my.cmems-du.eu/Core/INSITU_GLO_TS_REP_OBSERVATIONS_013_001_b/CORIOLIS-GLOBAL-EasyCORA-OBS have been downloaded 
# only Argo profiles PR_PF files are used here 
 
import copernicusmarine 
import matplotlib.pyplot as plt
import matplotlib.dates as md
import datetime as dt
import numpy as np
import gsw
from netCDF4 import Dataset, MFDataset, stringtochar
import netCDF4 as nc
import MLD as mld
from holteandtalley import HolteAndTalley
import pickle
plt.ion()



data_request = {
   ...:    "dataset_id_bathy" : "cmems_mod_glo_phy_anfc_0.083deg_static",
   ...:    "longitude" : [-180, 179.91668701171875],
   ...:    "latitude" : [-80, 90],
   ...:    "time" : ["2023-01-01", "2023-01-31"],
   ...:    "variables" : ["deptho", "deptho_lev","mask"]
   ...: }

bathy= copernicusmarine.open_dataset(
    dataset_id = data_request["dataset_id_bathy"],
    minimum_longitude = data_request["longitude"][0],
    maximum_longitude = data_request["longitude"][1],
    minimum_latitude = data_request["latitude"][0],
    maximum_latitude = data_request["latitude"][1],
    minimum_depth=0.49402499198913574,
    maximum_depth=0.49402499198913574,
    variables = data_request["variables"]
)

bathylat=bathy.variables['latitude'][:]
bathylon=bathy.variables['longitude'][:]
bathybat=bathy.variables['deptho'][:]


#smallest model domain
glo_request = {
    "dataset_id_glo" : "cmems_mod_glo_phy_anfc_0.083deg_P1D-m",
    "longitude" : [-180, 179.91668701171875],
    "latitude" : [-80, 90],
    "time" : ["2024-03-12T00:00:00", "2024-03-13T00:00:00"],
    "variables" : ["mlotst"]
}

mob_request = {
    "dataset_id_mob" : "dataset-armor-3d-nrt-weekly",
    "time" : ["2024-03-12T00:00:00", "2024-03-13T00:00:00"],
    "variables" : ["mlotst", "so", "to"]
}



# set parameter for de BoyerMontegut MLD calculation and pressure reference (to match what is in model)
presref=10.
denscrit=0.01
tempcrit=0.2

C4obs_dmld=[]
C4obs_tmld=[]
C4obs_smld=[]
C4GLO12mld=[]
C4MOBmld=[] 
C4time=[]
C4lat=[]
C4lon=[]
GLO12distance=[]
MOBdistance=[]
glo12res=0.08333588
mobres=0.25  
path='/data/users/cpequign/data/CORA/global_013_001/'  # path to observations daily files 

sizedeg=5    # bin size for averaging and produce time series
lonbin=np.linspace(-180,180,int(360/sizedeg)+1)
latbin=np.linspace(-80,90,int(170/sizedeg)+1) 

start_date = dt.datetime(2022, 9, 1)
end_date = dt.datetime(2022, 9, 30)
nbdays=(end_date-start_date).days
tt = np.arange(start_date, end_date, dt.timedelta(days=1)).astype(dt.datetime)

min_lon = glo_request["longitude"][0]
max_lon = glo_request["longitude"][1]
min_lat = glo_request["latitude"][0]
max_lat = glo_request["latitude"][1]

lonbin=np.linspace(-180,180,int(360/sizedeg)+1)
latbin=np.linspace(-80,90,int(170/sizedeg)+1) 

# to keep matching MLD variables
C4obs_dmld=[]
C4obs_tmld=[]
C4obs_smld=[]
C4obs_dmld2=[]
C4obs_tmld2=[]
C4GLO12mld=[]
C4MOBmld=[] 
C4time=[]
C4lat=[]
C4lon=[]
GLO12distance=[]
MOBdistance=[]

# to keep binned values
C4obs_dmld2_binned=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_tmld2_binned=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_dmld_binned=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_tmld_binned=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_smld_binned=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_dmld2_std=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_tmld2_std=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_dmld_std=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_tmld_std=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4obs_smld_std=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4GLO12mld_binned=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4MOBmld_binned=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4GLO12mld_std=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
C4MOBmld_std=np.nan*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))
count=0*np.ones((int((max_lon-min_lon)/sizedeg),int((max_lat-min_lat)/sizedeg),nbdays))   # count number of obs in bin

id_date=f"{start_date.year}-{start_date.month:02d}-{start_date.day:02d}"
# Download loop 
loop_date=start_date
daycount=0

while loop_date < end_date:   # loop for each day
    # observation file
    print(loop_date)
    out_name = f"{loop_date.year}/ECO_DMQCGL01_{loop_date.year}{loop_date.month:02d}{loop_date.day:02d}_PR_PF.nc"
    profile = Dataset(path+out_name,'r')
    temp=profile.variables['TEMP'][:].squeeze()
    sal=profile.variables['PSAL'][:].squeeze()
    pres=profile.variables['PRES'][:].squeeze()
    dens=profile.variables['DENS'][:].squeeze() 
    loop_datep1 = loop_date+dt.timedelta(days=1)   ### remember MOB is weekly  ########
    id_date=f"{loop_date.year}-{loop_date.month:02d}-{loop_date.day:02d}"
    id_datep1=f"{loop_datep1.year}-{loop_datep1.month:02d}-{loop_datep1.day:02d}"
    cm_date= f"{loop_date.year}-{loop_date.month:02d}-{loop_date.day:02d}T00:00:00"    # start_datetime="2024-03-13T00:00:00",
    print(cm_date)
   # model file

    ####### Read product to be assessed
    GLO12= copernicusmarine.open_dataset(
        dataset_id = glo_request["dataset_id_glo"],
        minimum_longitude = glo_request["longitude"][0],
        maximum_longitude = glo_request["longitude"][1],
        minimum_latitude = glo_request["latitude"][0],
        maximum_latitude = glo_request["latitude"][1],
        start_datetime= cm_date,
        end_datetime= cm_date,
        minimum_depth=0.49402499198913574,
        maximum_depth=0.49402499198913574,
        variables = glo_request["variables"]
    )
    MOB= copernicusmarine.open_dataset(
        dataset_id = mob_request["dataset_id_mob"],
        minimum_longitude = min_lon,
        maximum_longitude = max_lon,
        minimum_latitude = min_lat,
        maximum_latitude = max_lat,
        start_datetime= cm_date,
        end_datetime= cm_date,
        minimum_depth=0.49402499198913574,
        maximum_depth=0.49402499198913574,
        variables = mob_request["variables"]
    )


    fit_dmld=np.nan**np.ones(np.shape(profile.variables['LATITUDE'][:]))
    fit_dmld2=np.nan**np.ones(np.shape(profile.variables['LATITUDE'][:]))
    fit_tmld=np.nan**np.ones(np.shape(profile.variables['LATITUDE'][:]))
    fit_tmld2=np.nan**np.ones(np.shape(profile.variables['LATITUDE'][:]))
    fit_smld=np.nan**np.ones(np.shape(profile.variables['LATITUDE'][:]))  
    GLO12_mld_match=np.nan*np.ones(np.shape(profile.variables['LATITUDE'][:]))
    MOB_mld_match=np.nan*np.ones(np.shape(profile.variables['LATITUDE'][:])) 
    distance=np.nan*np.ones(np.shape(profile.variables['LATITUDE'][:])) 
    distance2=np.nan*np.ones(np.shape(profile.variables['LATITUDE'][:])) 
 
    problemprofile=np.array([])
    nbproblem=0

    for ii,ji in np.ndenumerate(profile.variables['LONGITUDE'][:]):   # loop over individual argo profile
        ## matching obs to model with closest point
        latmindiff = abs(GLO12.variables['latitude'][:]-(profile.variables['LATITUDE'][ii]))
        lonmindiff = abs(GLO12.variables['longitude'][:]-(profile.variables['LONGITUDE'][ii]))
        ilatmin = np.argmin(np.array(latmindiff[:]))
        ilonmin = np.argmin(np.array(lonmindiff[:]))
        ilat = GLO12.variables['latitude'][ilatmin].to_numpy()   # float() or .item()
        ilon = GLO12.variables['longitude'][ilonmin].to_numpy()
        distance[ii]=gsw.distance([profile.variables['LONGITUDE'][ii], ilon], [profile.variables['LATITUDE'][ii], ilat], p=0)/1000
        GLO12_mld_match[ii]=GLO12.mlotst[0,ilatmin,ilonmin]

        latmindiff2 = abs(MOB.variables['latitude'][:]-(profile.variables['LATITUDE'][ii]))
        lonmindiff2 = abs(MOB.variables['longitude'][:]-(profile.variables['LONGITUDE'][ii]))
        ilatmin2 = np.argmin(np.array(latmindiff2[:]))
        ilonmin2 = np.argmin(np.array(lonmindiff2[:]))
        ilat2 = MOB.variables['latitude'][ilatmin2].to_numpy()
        ilon2 = MOB.variables['longitude'][ilonmin2].to_numpy()
        distance2[ii]=gsw.distance([profile.variables['LONGITUDE'][ii], ilon2], [profile.variables['LATITUDE'][ii], ilat2], p=0)/1000
        MOB_mld_match[ii]=MOB.mlotst[0,ilatmin2,ilonmin2]
    
        # calculating observation MLD 
        subset=(pres[ii,:]*temp[ii,:]*sal[ii,:]).mask.squeeze() == False
        depthratio=np.diff(pres[ii,subset])[1::]/pres[ii,subset][1:-1]
        if len(pres[ii,subset])>5 and (max(depthratio)<2.5):
            if (min(pres[ii,subset]) <= presref and max(pres[ii,subset]) >= 2*presref) :
                try:
                    tpro=temp[ii,subset].squeeze()
                    spro=sal[ii,subset].squeeze()
                    dpro=dens[ii,subset].squeeze()
                    ppro=pres[ii,subset].squeeze()
                    only=np.where(ppro>presref)
                    dmld2,tmld2=mld.zmld_boyervar(spro[only], tpro[only], ppro[only],denscriteria=denscrit,pressref=presref)
                    fit_dmld2[ii]=dmld2
                    fit_tmld2[ii]=tmld2
                    h=HolteAndTalley(ppro[only], tpro[only], salinities=spro[only].tolist(), densities=dpro[only].tolist())
                    fit_dmld[ii]=h.densityMLD
                    fit_tmld[ii]=h.tempMLD
                    fit_smld[ii]=h.salinityMLD
                except:
                    print('An exception occured on profile ', ii)
                    problemprofile=np.concatenate((problemprofile,ii))
                    nbproblem=nbproblem+1

    C4obs_dmld.append(fit_dmld[:])
    C4obs_tmld.append(fit_tmld[:])
    C4obs_smld.append(fit_smld[:])       
    C4obs_dmld2.append(fit_dmld2[:])
    C4obs_tmld2.append(fit_tmld2[:])
    C4GLO12mld.append(GLO12_mld_match[:])   
    C4MOBmld.append(MOB_mld_match[:])
    C4lat.append(profile.variables['LATITUDE'][:])
    C4lon.append(profile.variables['LONGITUDE'][:])
    C4time.append(profile.variables['TIME'][:])
    MOBdistance.append(distance2[:])
    GLO12distance.append(distance[:])

    # build bin values in 5degree boxes time series 
    for count_ii, ii_lat in enumerate(latbin[0:-1]):
        #print(count_ii)    
        stii=count_ii*int(sizedeg/glo12res)
        enii=(count_ii+1)*int(sizedeg/glo12res)-1    
        stii2=count_ii*int(sizedeg/mobres)
        enii2=(count_ii+1)*int(sizedeg/mobres)-1     
        for count_jj, jj_lon in enumerate(lonbin[0:-1]):
            #print(count_jj)	
            stjj=count_jj*int(sizedeg/glo12res)
            enjj=(count_jj+1)*int(sizedeg/glo12res) -1
            stjj2=count_jj*int(sizedeg/mobres)
            enjj2=(count_jj+1)*int(sizedeg/mobres) -1	    
            take=np.where((np.array(profile.variables['LATITUDE'][:])>=ii_lat) & (np.array(profile.variables['LATITUDE'][:])<ii_lat+sizedeg) & (np.array(profile.variables['LONGITUDE'][:])>=jj_lon) & (np.array(profile.variables['LONGITUDE'][:])<jj_lon+sizedeg) )
            count[count_jj,count_ii,daycount]=np.shape(take)[1]
            C4GLO12mld_binned[count_jj,count_ii,daycount]=np.nanmean(GLO12.mlotst[0,stii:enii,stjj:enjj])
            C4MOBmld_binned[count_jj,count_ii,daycount]=np.nanmean(MOB.mlotst[0,stii2:enii2,stjj2:enjj2])
            C4GLO12mld_std[count_jj,count_ii,daycount]=np.nanstd(GLO12.mlotst[0,stii:enii,stjj:enjj])
            C4MOBmld_std[count_jj,count_ii,daycount]=np.nanstd(MOB.mlotst[0,stii2:enii2,stjj2:enjj2])	
            if np.shape(take)[1] != 0:
                C4obs_dmld_binned[count_jj,count_ii,daycount]=np.nanmean(fit_dmld[take])
                C4obs_tmld_binned[count_jj,count_ii,daycount]=np.nanmean(fit_tmld[take])
                C4obs_smld_binned[count_jj,count_ii,daycount]=np.nanmean(fit_smld[take])
                C4obs_dmld2_binned[count_jj,count_ii,daycount]=np.nanmean(fit_dmld2[take])
                C4obs_tmld2_binned[count_jj,count_ii,daycount]=np.nanmean(fit_tmld2[take])
                if np.shape(take)[1] != 1:
                    C4obs_dmld_std[count_jj,count_ii,daycount]=np.nanstd(fit_dmld[take])
                    C4obs_tmld_std[count_jj,count_ii,daycount]=np.nanstd(fit_tmld[take])
                    C4obs_smld_std[count_jj,count_ii,daycount]=np.nanstd(fit_smld[take])
                    C4obs_dmld2_std[count_jj,count_ii,daycount]=np.nanstd(fit_dmld2[take])
                    C4obs_tmld2_std[count_jj,count_ii,daycount]=np.nanstd(fit_tmld2[take])		

    daycount=daycount+1 
    loop_date=loop_datep1      ####################### end of daily loop


#save  in pickle files
match_obj = [C4obs_dmld, C4obs_tmld, C4obs_smld, C4GLO12mld, C4MOBmld, C4lat, C4lon, C4time, MOBdistance, GLO12distance]

with open('match_1month_sept2022.pickle', 'wb') as f:  
        pickle.dump(match_obj, f)
    
binned_obj = [C4GLO12mld_binned, C4MOBmld_binned, C4GLO12mld_std, C4MOBmld_std,C4obs_dmld_binned, C4obs_tmld_binned, C4obs_smld_binned, C4obs_dmld_std, C4obs_tmld_std, C4obs_smld_std]		

with open('binned_1month_sept2022.pickle', 'wb') as f:  
        pickle.dump(binned_obj, f)



#### if needed to retrieve output:

#with open('match_1month_sept2022.pickle', 'rb') as f:  
#    C4obs_dmld, C4obs_tmld, C4obs_smld, C4GLO12mld, C4MOBmld, C4lat, C4lon, C4time, MOBdistance, GLO12distance = pickle.load(f)
 
#with open('binned_1month_sept2022.pickle', 'rb') as f:  
#    C4GLO12mld_binned, C4MOBmld_binned, C4GLO12mld_std, C4MOBmld_std,C4obs_dmld_binned, C4obs_tmld_binned, C4obs_smld_binned, C4obs_dmld_std, C4obs_tmld_std, C4obs_smld_std = pickle.load(f)


### plot time series, at one 5x5 box, of binned MLD values and showing when there wasn't much data points 

count2p=np.copy(count)
count2p[(count<2)]=np.nan  # to mark data when more than one profile is averaged in a 5x5 degree bin

plt.figure()
plt.plot(tt,C4obs_tmld_binned[31,20,:]*count2p[31,20,:]/count2p[31,20,:],'bo',label='tMLD >2pts')
plt.plot(tt,C4obs_dmld_binned[31,20,:]*count2p[31,20,:]/count2p[31,20,:],'o',label='dMLD >2pts',color='orange')
plt.plot(tt,C4obs_tmld_binned[31,20,:],'*',label='tMLD Argo')       
plt.plot(tt,C4obs_dmld_binned[31,20,:],'.',label='dMLD Argo')
plt.plot(tt,C4GLO12mld_binned[31,20,:],label='GLO12')
plt.plot(tt,C4MOBmld_binned[31,20,:],'--',label='MOB')
plt.legend()


##map of number of points per bin for a given day

fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(2, 2)
s1=ax1.pcolormesh(lonbin,latbin,count[:,:,0].transpose(),cmap='viridis_r')
fig.colorbar(s1, ax=ax1)
plt.contour(bathylon,bathylat,bathybat.squeeze(),levels=[0,0.0001], colors='grey',linewidths=0.4)

ax1.pcolormesh(lonbin,latbin,count[:,:,0].transpose(),cmap='viridis_r')
cbar = plt.colorbar()
cbar.ax.set_ylabel('nb point')
plt.contour(bathylon,bathylat,bathybat.squeeze(),levels=[0,0.0001], colors='grey',linewidths=0.4)


### plot of normalized difference for a given day

plt.figure()
plt.scatter(profile.variables['LONGITUDE'][:],profile.variables['LATITUDE'][:],10,c=(C4obs_dmld[0]-C4GLO12mld[0])/C4obs_dmld[0],cmap='RdYlBu_r',vmin=-1, vmax=1)
plt.colorbar(extend='both')
plt.contour(bathylon,bathylat,bathybat.squeeze(),levels=[0,0.0001], colors='grey',linewidths=0.4)



#### same plot of normalized difference  points that are more than 100% off blacked out
### 
aa=(C4obs_dmld[0]-C4GLO12mld[0])/C4obs_dmld[0]
bb=aa[aa<-1]
#There are some points that are completely off, and it is misleading as it is often problematic profiles, so this is my way to remove them. 
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
plt.scatter(lon_all[:],lat_all[:],20,c=(C4obs_dmld[0]-C4GLO12mld[0])/C4obs_dmld[0],cmap='RdYlBu_r',vmin=-1, vmax=1)
plt.plot(lon_all[aa<-1],lat_all[aa<-1],'ok')

cbar=plt.colorbar()
#cbar.ax.set_ylabel('MLD (m)')
plt.ylim([-25,30])
plt.xlim([26,120])
plt.contour(bathylon,bathylat,bathybat.squeeze(),levels=[0,0.0001], colors='grey',linewidths=0.4)
major_ticks = np.arange(30, 120, 20)
minor_ticks = np.arange(30, 120, 5)
major_yticks = np.arange(-20, 30, 20)
minor_yticks = np.arange(-20, 30, 5)
ax.set_xticks(major_ticks)
ax.set_xticks(minor_ticks, minor=True)
ax.set_yticks(major_yticks)
ax.set_yticks(minor_yticks, minor=True)
# And a corresponding grid
ax.grid(which='both')

# Or if you want different settings for the grids:
ax.grid(which='minor', alpha=0.2)
ax.grid(which='major', alpha=0.5)
