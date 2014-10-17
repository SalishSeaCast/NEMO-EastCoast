% NAME: gen_obc.m
%
% AUTHOR: J.-P. Paquin 
%
% DATE: July 2013
%
% REVISIONS: Feb14 - J.-P. Paquin: 
%                    Add processing of T, S & SSH within this script
%                    Generalized to create all necessary OBC (N,S,E,W)
%            Jan14 - J.-P. Paquin: 
%                    Optimization of the code in f_opa_angle.m
%                    Optimization of the f_rot_rep
%                    Add loop over time
%         
%
% DESCRIPTION: Generate Open boundary conditions for Regional NEMO ("DESTINATION"). 
%              Developped to take data from GLORYS version 2v3 Data  ("SOURCE")
%              from ${datastorage}/DATA/GLORYS2V3
%           
%
% HYPOTHESES:  
%
%
% INPUTS:  INPUT dataset to extract the OBC from (SOURCE grid)
%            OPTION1: coordinates of SOURCE grid
%            OPTION2: files containing latitudes and longitudes for
%                     U, V, T and F on SOURCE grid 
%            DATA COVERAGE : 1993-2011 (monthly or annual)
%
%
%          DESTINATION 'NEMO regional' grid: coordinates, bathymetry, mesh_mask 
%
%
% OUTPUTS: Interpolated open boundary conditions for U and V on DESTINATION grid 
%
%
% CALLED PGM & SCRIPTS: f_open_netcdf
%                       f_readnetcdf
%                       f_writenetcdf
%
%                    PROCESSING T, S, and SSH
%                       f_convectTS (floodnan_opa3; rho_wright)
%
%                    PROCESSING U and V
%                       f_interpuv (->f_rot_rep->f_opa_angle)
%
%
% NOTES: 
%     1- REQUIRES M_MAP PACKAGE (version 1.4) for Matlab
%        THE VARIABLE time_counter is VERYIMPORTANT FOR NEMO TO KNOW WHEN 
%        TO READ NEW OBC RECORDS AND DEFINED AS THE NUMBER OF DAYS SINCE
%        1950-01-01 (see f_writenetcdf_obc.m)
%
%     2- North and East OBCs for CONCEPTS version 1.1.0 must be inverted to
%        have index 1 at the outer layer increasing inward
%
%   ***********************************************************************
%   *******  DURING INTERPOLATION, NaNs ARE REPLACED BY 0 !!!!!!!!  *******
%   ***********************************************************************
%
%   ***********************************************************************
%   ******* NETCDF VERSION SUPPORTED 3.6.2 : OLDER VERSION UNTESTED *******
%   *******             SEE KEY 'cdfversion' IN CODE                *******
%   ***********************************************************************
%
%--------------------------------------------------------------------------
clear all
display('***** CALL TO : gen_obc.m *****')
addpath('/media/Data/NEMO/JPP/0_common/m_map')

    
%##############################################################
%################## BEGIN USER MODIFICATIONS ##################
% - DESTINATION grid files required
dirbathy='/media/Data/NEMO/FC/NEMO_input/Pegasus/';
filebathy='Capesable_bathymetry_NWATL_ORCA05km.nc';
dircoord='/media/Data/NEMO/FC/NEMO_input/Pegasus/';
filecoord='Capesable_coordinates_FromORCA_05km.nc';
dirmesh ='/media/Data/NEMO/FC/NEMO_input/Pegasus/';
filemesh='CPSBL_mesh_mask.nc';
config_name='CPSBL';


% - SOURCE grid files/options required
inputdata='GLORYS2V3'; % GLORYS2 version 3
datadir=(['/media/Data/NEMO/Data/' inputdata]);

% OPTION 1 : coordinate file for SOURCE grid
dircoordG='';             
filecoordG='coord';

% OPTION 2 : load latitudes and longitudes form files of U,V (T and F) 
% TIME PERIOD TO EXTRACT
stayr=1993;
endyr=1993;
%list_months={'01','02','03','04','05','06','07','08','09','10','11','12'};
list_months={'01','02'};
[~,nbmonths]=size(list_months);


% - Input parameters
cdfversion='3.6.2'; %- Supported netcdf version: v3.6.2
                    %  Older versions untested

concepts_version='1.1.0'; % Switch for east and north boundary inversion
                          % Version 1.1.0           : inversion required
                          % Version 2.0.0 and above : not required 

list_obc={'north','south','east','west'}; % - List of open boundaries in the domain
                     
nbdeg=2; %- number of degrees around the DESTINATION grid required for 
         %  data extraction from the SOURCE grid to avoid (as much as
         %  possible) memory problems

nbptsbound=10; %- number of points in the boundary (standard is 10)
m_proj('mercator')
%##############################################################

time_genUV_beg=clock;                                                      %jpp


% -------------------------------------------------------------------------
% [1a] Bathymetry, mesh mask, coordinates
% -------------------------------------------------------------------------
display(' 1 - LOAD Bathy, coordinates, etc. ')
% use Bathy to get land-sea mask
if strcmp(cdfversion,'3.6.2')
  tmp_data=f_open_netcdf(dirbathy,filebathy,'Bathymetry');
  Bathy=permute(tmp_data,[2,1]); 
  clear tmp_data;
else
  display('**** OLDER NETCDF VERSION UNTESTED ****')  
  nc=netcdf('/home/zhaili/ESRF/NEMO/Nesting_tools/bin/1_bathy_OCRA12_GB.nc','nowrite');
  Bathy=nc{'Bathymetry'}(:); 
  close(nc);
end

% read information from the mesh mask (mdepth required)
[tmask,umask,vmask,fmask,mdepth,gdepw,dzz,~,~]...
                        =f_readmeshmask([dirmesh '/' filemesh],cdfversion);

% determine full cells (mlevel) by checking the bathy_meter file
[nz,ny,nx]=size(tmask);
mlevel=zeros(ny,nx);
for k=1:nz-1
    mlevel(Bathy>gdepw(k)) = k+1;
end
clear Bathy nx ny nz


% [1b] Create structure for coordinates (instead of path to coordinates files)
%      If no coordinate file exists for the SOURCE grid, use the latitudes 
%      and longitudes for T,U,V,F files 
%      *** N.B. Necessary to define all coordinates for use in function ***
%      ***      f_opa_angle.m                                           ***

% - load coordinates of DESTINATION grid
Rglamt=f_readnetcdf(dircoord,filecoord,'glamt',cdfversion);  glamt=squeeze(Rglamt);
Rglamu=f_readnetcdf(dircoord,filecoord,'glamu',cdfversion);  glamu=squeeze(Rglamu);
Rglamv=f_readnetcdf(dircoord,filecoord,'glamv',cdfversion);  glamv=squeeze(Rglamv);
Rglamf=f_readnetcdf(dircoord,filecoord,'glamf',cdfversion);  glamf=squeeze(Rglamf);
Rgphit=f_readnetcdf(dircoord,filecoord,'gphit',cdfversion);  gphit=squeeze(Rgphit);
Rgphiu=f_readnetcdf(dircoord,filecoord,'gphiu',cdfversion);  gphiu=squeeze(Rgphiu); 
Rgphiv=f_readnetcdf(dircoord,filecoord,'gphiv',cdfversion);  gphiv=squeeze(Rgphiv); 
Rgphif=f_readnetcdf(dircoord,filecoord,'gphif',cdfversion);  gphif=squeeze(Rgphif); 

%- Definition of a structure variable containing all DESTINATION grid lats/lons
coordREG=struct('glamt',glamt,'glamu',glamu,'glamv',glamv,'glamf',glamf, ...
                'gphit',gphit,'gphiu',gphiu,'gphiv',gphiv,'gphif',gphif  );
clear Rglamt Rglamu Rglamv Rglamf Rgphit Rgphiu Rgphiv Rgphif glamf gphif
            
       
% - load coordinates information of SOURCE grid
if ~exist([ dircoordG '/' filecoordG ],'file');
  display(' GLOBAL COORDINATE FILE NOT PROVIDED - READS T,U,V,F files')
  dirin=([datadir '/MONTHLY_' num2str(stayr)]);
  
  Tfile=([ inputdata '_ORCA025_' num2str(stayr) '01_gridT.nc']);
  Ufile=([ inputdata '_ORCA025_' num2str(stayr) '01_gridU.nc']);
  Vfile=([ inputdata '_ORCA025_' num2str(stayr) '01_gridV.nc']);
  Ffile=([ inputdata '_ORCA025_' num2str(stayr) '01_gridT.nc']); %- do not have f coord for GLORYS
  
  Gglamt=f_readnetcdf(dirin,Tfile,'nav_lon',cdfversion);  Gglamt=squeeze(Gglamt); 
  Gglamu=f_readnetcdf(dirin,Ufile,'nav_lon',cdfversion);  Gglamu=squeeze(Gglamu); 
  Gglamv=f_readnetcdf(dirin,Vfile,'nav_lon',cdfversion);  Gglamv=squeeze(Gglamv); 
  Gglamf=f_readnetcdf(dirin,Tfile,'nav_lon',cdfversion);  Gglamf=squeeze(Gglamf); 
  Ggphit=f_readnetcdf(dirin,Ffile,'nav_lat',cdfversion);  Ggphit=squeeze(Ggphit);   
  Ggphiu=f_readnetcdf(dirin,Ufile,'nav_lat',cdfversion);  Ggphiu=squeeze(Ggphiu);
  Ggphiv=f_readnetcdf(dirin,Vfile,'nav_lat',cdfversion);  Ggphiv=squeeze(Ggphiv);
  Ggphif=f_readnetcdf(dirin,Ffile,'nav_lat',cdfversion);  Ggphif=squeeze(Ggphif);
else
  display(' GLOBAL COORDINATE FILE PROVIDED')
  Gglamt=f_readnetcdf(dircoordG,filecoordG,'glamt',cdfversion);  Gglamt=squeeze(Gglamt); 
  Gglamu=f_readnetcdf(dircoordG,filecoordG,'glamu',cdfversion);  Gglamu=squeeze(Gglamu);
  Gglamv=f_readnetcdf(dircoordG,filecoordG,'glamv',cdfversion);  Gglamv=squeeze(Gglamv);
  Gglamf=f_readnetcdf(dircoordG,filecoordG,'glamf',cdfversion);  Gglamf=squeeze(Gglamf);
  Ggphit=f_readnetcdf(dircoordG,filecoordG,'gphit',cdfversion);  Ggphit=squeeze(Ggphit); 
  Ggphiu=f_readnetcdf(dircoordG,filecoordG,'gphiu',cdfversion);  Ggphiu=squeeze(Ggphiu); 
  Ggphiv=f_readnetcdf(dircoordG,filecoordG,'gphiv',cdfversion);  Ggphiv=squeeze(Ggphiv);
  Ggphif=f_readnetcdf(dircoordG,filecoordG,'gphif',cdfversion);  Ggphif=squeeze(Ggphif);
end
%- Definition of a structure variable containing all SOURCE grid lats/lons
coordGLOB=struct('glamt',Gglamt,'glamu',Gglamu,'glamv',Gglamv,'glamf',Gglamf, ...
                 'gphit',Ggphit,'gphiu',Ggphiu,'gphiv',Ggphiv,'gphif',Ggphif  );    
clear Gglamf Ggphif Gglamu Gglamv Gglamf Ggphit Ggphiu Ggphiv Ggphif
    

% -------------------------------------------------------------------------
% MAIN LOOP OVER THE USER DEFINED OPEN BOUNDARIES
global obc_name ;
for myobc = list_obc
  obc_name = myobc{1}; myobc;
  dirin=([datadir '/MONTHLY_' num2str(stayr)]);

% -------------------------------------------------------------------------
% [2] OPEN BOUNDARY LOCATION
%     The idea is to extract from the global files only the region
%     surrounding the open boundary (memory optimization).
% ** TEST FOR WESTRN BOUNDARY CONDITIONS ONLY :: PLAN TO GENERALIZE FOR **
% ** THE FOUR OBCS INCLUDING THE OPTION OF PERMUTING FOR NORTH AND EAST **
%
% N.B. The open boundary width is set by the parameter nbptsbound and must
%      later match the definition of the lateral boundary in NEMO code
% -------------------------------------------------------------------------
display(' 2 - DEFINE OPEN BOUNDARY')
[RdomY,RdomX]=size(coordREG.glamt);
%global obc_name ; obc_name='north'; 
[corners]=f_define_corners(RdomY,RdomX,obc_name,nbptsbound);

%- [2-2] Find minimum/maximum latitudes and longitudes of the open water
%        section on DESTINATION grid (based on the location of the 4
%        corners)
minlat=min([gphit(corners(1,1),corners(1,2)),gphit(corners(2,1),corners(2,2)), ...
            gphit(corners(4,1),corners(4,2)),gphit(corners(4,1),corners(4,2))]);
maxlat=max([gphit(corners(1,1),corners(1,2)),gphit(corners(2,1),corners(2,2)), ...
            gphit(corners(4,1),corners(4,2)),gphit(corners(4,1),corners(4,2))]);      
minlon=min([glamt(corners(1,1),corners(1,2)),glamt(corners(2,1),corners(2,2)), ...
            glamt(corners(4,1),corners(4,2)),glamt(corners(4,1),corners(4,2))]);
maxlon=max([glamt(corners(1,1),corners(1,2)),glamt(corners(2,1),corners(2,2)), ...
            glamt(corners(4,1),corners(4,2)),glamt(corners(4,1),corners(4,2))]);     

        
% -------------------------------------------------------------------------
% [3] Data preparation for call to f_interpuv (M. Dunphy) 
% -------------------------------------------------------------------------
% [3a] -- LOAD FIRST U FILE TO FIND THE MIN/MAX LATITUDES AND LONGITUDES
%         TO EXTRACT ON THE LARGE SCALE DATA
display(' 3a - IDENTIFY ZONE IN SOURCE GRID')
%- Read vertical levels of input file
globU    =f_readnetcdf(dirin,Ufile,'vozocrtx',cdfversion); globU=squeeze(globU);
dpthglobU=f_readnetcdf(dirin,Ufile,'deptht',cdfversion)  ; dpthglobU=squeeze(dpthglobU);
dpthglobV=f_readnetcdf(dirin,Vfile,'deptht',cdfversion)  ; dpthglobV=squeeze(dpthglobV);
[NZ,NY,NX]=size(globU);
 
%- Search for minimum and maximum indexes within the SOURCE data
%corresponding to the DESTINATION area of interest (high res domain +- nbdeg)
index_minii=NX;         index_maxii=-NX;
index_minjj=NY;         index_maxjj=-NY;
for globii=1:NX
  for globjj=1:NY
    % search if point into desired zone
    if ( (coordGLOB.gphiu(globjj,globii)>=(minlat-nbdeg)) && (coordGLOB.gphiu(globjj,globii)<=(maxlat+nbdeg)) )    
      if ( (coordGLOB.glamu(globjj,globii)>=(minlon-nbdeg)) && (coordGLOB.glamu(globjj,globii)<=(maxlon+nbdeg)) )    
         if ( globii <= index_minii ) ; index_minii=globii ; end  
         if ( globii >= index_maxii ) ; index_maxii=globii ; end  
         if ( globjj <= index_minjj ) ; index_minjj=globjj ; end  
         if ( globjj >= index_maxjj ) ; index_maxjj=globjj ; end  
      end
    end
  end
end
zoomdimj=index_maxjj-index_minjj+1;
zoomdimi=index_maxii-index_minii+1;

% -- Extract "regional" latitudes / longitudes of the SOURCE grid
reg_glob_latT=coordGLOB.gphit(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_lonT=coordGLOB.glamt(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_latU=coordGLOB.gphiu(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_lonU=coordGLOB.glamu(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_latV=coordGLOB.gphiv(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_lonV=coordGLOB.glamv(index_minjj:index_maxjj,index_minii:index_maxii);
clear globii globjj 


% [3b] -- LOOP OVER TIME TO EXTRACT OPEN BOUNDARY CONDITIONS 
display(' 3b - LOOP OVER MONTHS')
time_genUV_beg_load=clock;                                                 %jpp
NBREC=(endyr-stayr+1)*(nbmonths);
reg_glob_dataSSH=zeros(1,zoomdimj,zoomdimi,NBREC);
reg_glob_dataT  =zeros(NZ,zoomdimj,zoomdimi,NBREC);
reg_glob_dataS  =zeros(NZ,zoomdimj,zoomdimi,NBREC);
reg_glob_dataU  =zeros(NZ,zoomdimj,zoomdimi,NBREC);
reg_glob_dataV  =zeros(NZ,zoomdimj,zoomdimi,NBREC);
time_counter=zeros(NBREC,1);
year=stayr;  tt=1;
while year <= endyr
  dirin=([datadir '/MONTHLY_' num2str(year)]);

  for mymonth = list_months
    mm = mymonth{1} ; mymonth;
    display(['      LOADING DATA FOR: ' num2str(year) '/' mm])
    
    % -- EXTRACT SSH
    display([ '      LOAD SSH data for Y=' num2str(year) ' M=' mm ])
    SSHfile=([ inputdata '_ORCA025_' num2str(stayr) mm '_grid2D.nc']);
    globSSH=f_readnetcdf(dirin,SSHfile,'sossheig',cdfversion); 
    globSSH=squeeze(globSSH);
    reg_glob_dataSSH(1,:,:,tt)=squeeze(globSSH(index_minjj:index_maxjj,index_minii:index_maxii));
    clear globSSH
    
    % -- EXTRACT T 
    display([ '      LOAD   T data for Y=' num2str(year) ' M=' mm ])
    Tfile=([ inputdata '_ORCA025_' num2str(stayr) mm '_gridT.nc']);
    globT=f_readnetcdf(dirin,Tfile,'votemper',cdfversion); 
    globT=squeeze(globT);
    reg_glob_dataT(:,:,:,tt)=squeeze(globT(:,index_minjj:index_maxjj,index_minii:index_maxii));
    clear globT

    % -- EXTRACT S
    display([ '      LOAD   S data for Y=' num2str(year) ' M=' mm ])
    Sfile=([ inputdata '_ORCA025_' num2str(stayr) mm '_gridS.nc']);
    globS=f_readnetcdf(dirin,Sfile,'vosaline',cdfversion); 
    globS=squeeze(globS);
    reg_glob_dataS(:,:,:,tt)=squeeze(globS(:,index_minjj:index_maxjj,index_minii:index_maxii));
    clear globS
    
    % -- EXTRACT U 
    display([ '      LOAD   U data for Y=' num2str(year) ' M=' mm ])
    Ufile=([ inputdata '_ORCA025_' num2str(stayr) mm '_gridU.nc']);
    globU=f_readnetcdf(dirin,Ufile,'vozocrtx',cdfversion); 
    globU=squeeze(globU);
    reg_glob_dataU(:,:,:,tt)=squeeze(globU(:,index_minjj:index_maxjj,index_minii:index_maxii));
    clear globU

    % -- EXTRACT V
    display([ '      LOAD   V data for Y=' num2str(year) ' M=' mm ])
    Vfile=([ inputdata '_ORCA025_' num2str(stayr) mm '_gridV.nc']);
    globV=f_readnetcdf(dirin,Vfile,'vomecrty',cdfversion); 
    globV=squeeze(globV);
    reg_glob_dataV(:,:,:,tt)=squeeze(globV(:,index_minjj:index_maxjj,index_minii:index_maxii));
    clear globV
     
    % Calculate the variable time_counter defined as : days since 1950-01-01
    % VERY IMPORTANT FOR THE MODEL, THAT IS HOW NEMO KNOWS WHEN TO 
    % READ/CHANGE OBC CONDITIONS
    
    datestr='01/01/1950';
    dateend=[str2num(mm),'/01/',year];
    time_counter(tt)=daysact(datestr,dateend);
    tt=tt+1;
    
  end

year=year+1;
end % -- END OF TIME LOOP

time_genUV_end_load=clock;



% -------------------------------------------------------------------------
% [4] Processing scalar fields - Use griddata for coherence with U,V
% -------------------------------------------------------------------------
time_beg_scalar=clock;                                                     %jpp
% [4a] --  Interpolating and flooding from SOURCE to DESTINATION 
display(' 4a - Interpolating and flooding')
xi=1:zoomdimi;
yi=1:zoomdimj;
axi=corners(1,2):corners(3,2);
ayi=corners(1,1):corners(2,1);

% FCH: Put NAN where we have zero tracer
% Otherwise we get funny results from interpolation of zero especially in salinity
reg_glob_dataT(reg_glob_dataT==0)=NaN;
reg_glob_dataS(reg_glob_dataS==0)=NaN;

[outT]  =f_interp_scalar(reg_glob_dataT,reg_glob_lonT,reg_glob_latT,...
                         coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
                         dpthglobV, mdepth);
[outS]  =f_interp_scalar(reg_glob_dataS,reg_glob_lonT,reg_glob_latT,...
                         coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
                         dpthglobV, mdepth);
[outSSH]=f_interp_scalar(reg_glob_dataSSH,reg_glob_lonT,reg_glob_latT,...
                         coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
                         dpthglobV, mdepth);
time_end_scalar=clock;                                                     %jpp

% [4b] ---  Stratification verification and Convection
time_beg_convect=clock;                                                    %jpp
display(' 4b - Stratification and convection #####  REQUIRE OPTIMIZATION  ######')
[nzz,nyy,nxx,ntt]=size(outT);
for ti=1:ntt
  display(['    CONVECT FOR TIME: ' num2str(ti) ])
  for ii=1:nxx
    for jj=1:nyy
      kmtm=mlevel(jj,ii); % number of full cells at this point
      if (kmtm>1) % call TS_convect1 on FULL CELLS ONLY

        ttin=squeeze(outT(1:kmtm,jj,ii,ti));
        ssin=squeeze(outS(1:kmtm,jj,ii,ti));
        z1=mdepth(1:kmtm);
        dz1=dzz(1:kmtm);

        [tcc,scc]=TS_convect1_fast(ttin,ssin,z1,dz1);

        outT(1:kmtm,jj,ii,ti)=tcc;
        outS(1:kmtm,jj,ii,ti)=scc;
      end
    end
  end
end
time_end_convect=clock;                                                    %jpp
                 

% -------------------------------------------------------------------------
% [5] call to interpolation (and vector rotation) function
% -------------------------------------------------------------------------
time_genUV_beg_interp=clock;                                               %jpp
%  reg_glob_dataU : U velocities on SOURCE grid         (over extracted region)
%  reg_glob_dataV : V velocities on SOURCE grid        (over extracted region)
%  reg_glob_lonU  : U-grid longitudes from SOURCE file (over extracted region)
%  reg_glob_latU  : U-grid latitudes  from SOURCE file (over extracted region)
%  reg_glob_lonV  : V-grid longitudes from SOURCE file (over extracted region)
%  reg_glob_lonU  : V-grid latitudes  from SOURCE file (over extracted region)
%  tgtlonu      : U-grid longitudes for DESTINATION model grid
%  tgtlatu      : U-grid latitudes  for DESTINATION model grid
%  tgtlonv      : V-grid longitudes for DESTINATION model grid
%  tgtlatv      : v-grid latitudes  for DESTINATION model grid
%  dpthglob     : Depth levels of the SOURCE file 
%  mdepth       : Model depths (from mesh mask)
%  coordGLOB    : structure of coordinate files from SOURCE grid 
%  xi           : indices of the sub-section of SOURCE data
%  yi           : indices of the sub-section of SOURCE data
%  coordREG     : structure of coordinate file from DESTINATION grid  
%  axi          : array size of DESTINATION file 
%  ayi          : array size of DESTINATION file 
% xi=1:zoomdimi;
% yi=1:zoomdimj;
% axi=corners(1,2):corners(3,2);
% ayi=corners(1,1):corners(2,1);
display('       WARNING NaNs ARE REPLACED BY 0 FOR U and V INTERPOLATION!!!')
reg_glob_dataU(isnan(reg_glob_dataU))=0;
reg_glob_dataV(isnan(reg_glob_dataV))=0;

%dbstop in f_interpuv.m % debugging mode in f_interpuv
[outU outV]=f_interpuv(reg_glob_dataU,reg_glob_dataV,...
                       reg_glob_lonU ,reg_glob_latU,reg_glob_lonV,reg_glob_latV, ...
                       coordREG.glamu(ayi,axi),coordREG.gphiu(ayi,axi),...
                       coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),...
                       dpthglobV, mdepth, ...
                       coordGLOB, xi, yi, ...
                       coordREG ,axi,ayi);
time_genUV_end_interp=clock;                                               %jpp

% -----  PROCESS MASKS FOR OUTPUT
maskTout1=tmask(1,ayi,axi);
maskTout =tmask(:,ayi,axi); 
maskUout =umask(:,ayi,axi); 
maskVout =vmask(:,ayi,axi); 

%-------------------------------------------------------------------------
% [A] WRITE NETCDF FILES FOR OBC 
%-------------------------------------------------------------------------
display('WRITE NETCDF FILES')
time_genUV_beg_write=clock;                                               %jpp
%dbstop in f_writenetcdf.m % debugging mode in f_interpuv
fname = (['obc_' obc_name '_sossheig_' inputdata '_' config_name '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'sossheig',...
                       coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
                       1,time_counter,maskTout1,outSSH,obc_name,concepts_version);

fname = (['obc_' obc_name '_votemper_' inputdata '_' config_name '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'votemper',...
                       coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
                       mdepth,time_counter,maskTout,outT,obc_name,concepts_version);

fname = (['obc_' obc_name '_vosaline_' inputdata '_' config_name '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'vosaline',...
                       coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
                       time_counter,maskTout,outS,obc_name,concepts_version);

fname = (['obc_' obc_name '_vozocrtx_' inputdata '_' config_name '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'vozocrtx',...
                       coordREG.glamu(ayi,axi),coordREG.gphiu(ayi,axi),mdepth,...
                       time_counter,maskUout,outU,obc_name,concepts_version);

fname = (['obc_' obc_name '_vomecrty_' inputdata '_' config_name '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'vomecrty',...
                       coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
                       time_counter,maskVout,outV,obc_name,concepts_version);
time_genUV_end_write=clock;                                                %jpp

end % END LOOP ON LIST_OBC


time_genUV_end=clock;                                                      %jpp
display('*******************************************')                     %jpp
display('*****      SHOW EXECUTION TIMES       *****')                     %jpp
display([ 'OVERALL TIME BEG :' num2str(time_genUV_beg)       ])            %jpp
display([ 'OVERALL TIME END :' num2str(time_genUV_end)       ])            %jpp
display([ '  LOAD TIME BEG  :' num2str(time_genUV_beg_load)  ])            %jpp
display([ '  LOAD TIME END  :' num2str(time_genUV_end_load)  ])            %jpp
display([ ' SCALAR TIME BEG :' num2str(time_beg_scalar)      ])            %jpp
display([ ' SCALAR TIME END :' num2str(time_end_scalar)      ])            %jpp
display([ 'CONVECT TIME BEG :' num2str(time_beg_convect)     ])            %jpp
display([ 'CONVECT TIME END :' num2str(time_end_convect)     ])            %jpp
display([ ' INTERP TIME BEG :' num2str(time_genUV_beg_interp)])            %jpp
display([ ' INTERP TIME END :' num2str(time_genUV_end_interp)])            %jpp
display([ '  WRITE TIME BEG :' num2str(time_genUV_beg_write) ])            %jpp
display([ '  WRITE TIME END :' num2str(time_genUV_end_write) ])            %jpp
display('*******************************************')                                                            %jpp

fprintf('Finished\n');
