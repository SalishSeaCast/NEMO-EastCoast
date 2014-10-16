% NAME: gen_ini.m
%
% AUTHOR: J.-P. Paquin 
%
% DATE: Feb 2014 
%
% REVISIONS: 
%
% DESCRIPTION: Generate Initial conditions for Regional NEMO ("DESTINATION"). 
%              Developped to take data from GLORYS version 2v3 Data  ("SOURCE")
%              from ${datastorage}/DATA/GLORYS2V3
%              - One script for all variables T,S,U,V,SSH
%
% HYPOTHESES:  
%
%
% INPUTS:  INPUT dataset to extract the initial conditions (IC)
%            OPTION1: coordinates of SOURCE grid
%            OPTION2: files containing latitudes and longitudes for
%                     U, V, T and F on SOURCE grid 
%            DATA COVERAGE : 1993-2011 (monthly or annual)
%
%          DESTINATION 'NEMO regional' grid: coordinates, bathymetry, mesh_mask 
%
%
% OUTPUTS: Interpolated initial conditions for U and V on DESTINATION grid 
%
%
% CALLED PGM & SCRIPTS: f_open_netcdf
%                       f_readnetcdf
%                       f_writenetcdf
%                       f_readmeshmask
%
%                    PROCESSING T, S, and SSH
%                       f_create_mask (for floodnan_opa3)
%                       f_interp_scalar (floodnan_opa3; interp1q)
%
%                    PROCESSING U and V
%                       f_interpuv (->f_rot_rep->f_opa_angle)
%
%
% NOTES: REQUIRES M_MAP PACKAGE (version 1.4) for Matlab
%
%
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
display('***** CALL TO : gen_ini.m *****')
addpath('/media/Data/NEMO/JPP/0_common/')

    
%##############################################################
%################## BEGIN USER MODIFICATIONS ##################
% - DESTINATION grid files required
% dirbathy='/users/staff/jppaquin/scratch/WC3_PREP/3_mesh_mask';
% filebathy='bathy_meter.nc';
% dircoord='/users/staff/jppaquin/scratch/WC3_PREP/3_mesh_mask';
% filecoord='coordinates.nc';
% dirmesh ='/users/staff/jppaquin/scratch/WC3_PREP/3_mesh_mask';
% filemesh='mesh_mask.nc';
% config_name='WC3';
% - DESTINATION grid files required
dirbathy='/media/Data/NEMO/FC/NEMO_input/Pegasus/';
filebathy='Capesable_bathymetry_NWATL_R_smoothed_ChannelNorth.nc';
dircoord=dirbathy;
filecoord='Capesable_coordinates_R.nc';
dirmesh =dirbathy;
filemesh='mesh_mask_NWATL_R_smoothed_ChannelNorth.nc';
config_name='Capesable_FromAnna';
outputdir='/media/Data/NEMO/FC/NEMO_OUT_PREP/Initial';

% - SOURCE grid files/options required
%inputdata='GLORYS2V3'; % GLORYS2 version 3
inputdata='ANNA'; %FC
%datadir=(['/users/staff/jppaquin/scratch/DATA/' inputdata]);
datadir='/media/Data/NEMO/Anna/Anna_output/SCOTIAN_RUN2010_NOTIDES/';

%Sourcedir='/users/students/fateme/butterfly/FC/Scratch/DATA/CREG12';   %FC

% FC: For Anna's results
Tfile='SCOTIAN_T_2010_Bcl_notides_Extracted_CapeSable.nc';
Sfile=Tfile;
Ufile='SCOTIAN_U_2010_Bcl_notides_Extracted_CapeSable.nc';
Vfile='SCOTIAN_V_2010_Bcl_notides_Extracted_CapeSable.nc';
SSHfile='SCOTIAN_2D_2010_Bcl_notides_Extracted_CapeSable.nc';

% OPTION 1 : coordinate file for SOURCE grid
dircoordG='/media/Data/NEMO/Anna/Anna_output';             
filecoordG='coordinates_SCOTIAN_NEW_ANNA_CPSBL_Extracted.nc';
% OPTION 2 : load latitudes and longitudes form files of U,V (T and F) 
% TIME PERIOD TO EXTRACT
year='2010';
month='01';

% - Input parameters
cdfversion='3.6.2';  %- Supported netcdf version: v3.6.2
                     %  Older versions untested
          
nbdeg=0.1; %- number of degrees around the DESTINATION grid required for 
         %  data extraction from the SOURCE grid to avoid (as much as
         %  possible) memory problems
nbpts=10;%- "Real" search radius for land extrapolation (number of points)

m_proj('mercator')
%##############################################################

time_glob_beg=clock;                                                       %jpp


% -------------------------------------------------------------------------
% [1a] Bathymetry, mesh mask, coordinates
% -------------------------------------------------------------------------
display(' 1 - LOAD Bathy, coordinates, etc. ')
% use Bathy to get land-sea mask
if strcmp(cdfversion,'3.6.2')
  tmp_data=f_open_netcdf(dirbathy,filebathy,'Bathymetry');
  Bathy=permute(tmp_data,[2,1]); 
  clear tmp_data;
end

% read information from the mesh mask 
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
  %dir=([datadir '/MONTHLY_' num2str(year)]);
  
  %Tfile=([ inputdata '_ORCA025_' num2str(year) month '_gridT.nc']);
  %Ufile=([ inputdata '_ORCA025_' num2str(year) month '_gridU.nc']);
  %Vfile=([ inputdata '_ORCA025_' num2str(year) month '_gridV.nc']);
  %Ffile=([ inputdata '_ORCA025_' num2str(year) month '_gridT.nc']); %- do not have f coord for GLORYS
  
  % Added by FC for CREG12
%   dir=Sourcedir;
%   Tfile='CREG12_y2003m01d03h00m00_gridT.nc';
%   Sfile=Tfile;
%   Ufile='CREG12_y2003m01d03h00m00_gridU.nc';
%   Vfile='CREG12_y2003m01d03h00m00_gridV.nc';
%   Ffile=Tfile; %do not have fcoord for CREG12
%   SSHfile='CREG12_y2003m01d03h00m00_gridT2D.nc';
  
  
  Gglamt=f_readnetcdf(dir,Tfile,'nav_lon',cdfversion);  Gglamt=squeeze(Gglamt); 
  Gglamu=f_readnetcdf(dir,Ufile,'nav_lon',cdfversion);  Gglamu=squeeze(Gglamu); 
  Gglamv=f_readnetcdf(dir,Vfile,'nav_lon',cdfversion);  Gglamv=squeeze(Gglamv); 
  Gglamf=f_readnetcdf(dir,Tfile,'nav_lon',cdfversion);  Gglamf=squeeze(Gglamf); 
  Ggphit=f_readnetcdf(dir,Ffile,'nav_lat',cdfversion);  Ggphit=squeeze(Ggphit);   
  Ggphiu=f_readnetcdf(dir,Ufile,'nav_lat',cdfversion);  Ggphiu=squeeze(Ggphiu);
  Ggphiv=f_readnetcdf(dir,Vfile,'nav_lat',cdfversion);  Ggphiv=squeeze(Ggphiv);
  Ggphif=f_readnetcdf(dir,Ffile,'nav_lat',cdfversion);  Ggphif=squeeze(Ggphif);
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
clear Gglamf Ggphif
    

% -------------------------------------------------------------------------
% [2] Look for region to extract from the global files only the region
%     surrounding the model domain (memory optimization).
% -------------------------------------------------------------------------
display(' 2 - LOOK FOR REGION TO EXTRACT')
minlat=min(min(coordREG.gphit));
maxlat=max(max(coordREG.gphit));
minlon=min(min(coordREG.glamt));
maxlon=max(max(coordREG.glamt));

     
        
% -------------------------------------------------------------------------
% [3] Data preparation for call to f_interpuv (M. Dunphy) 
% -------------------------------------------------------------------------
% [3a] -- LOAD FIRST U FILE TO FIND THE MIN/MAX LATITUDES AND LONGITUDES
%         TO EXTRACT ON THE LARGE SCALE DATA
display(' 3a - IDENTIFY ZONE IN SOURCE GRID')
%- Read vertical levels of input file
globU    =f_readnetcdf(datadir,Ufile,'vozocrtx',cdfversion); globU=squeeze(globU);
%dpthglobU=f_readnetcdf(dir,Ufile,'deptht',cdfversion)  ; dpthglobU=squeeze(dpthglobU);
dpthglobU=f_readnetcdf(datadir,Ufile,'depthu',cdfversion)  ; dpthglobU=squeeze(dpthglobU); %FC for CREG12

%dpthglobV=f_readnetcdf(dir,Vfile,'deptht',cdfversion)  ; dpthglobV=squeeze(dpthglobV);
dpthglobV=f_readnetcdf(datadir,Vfile,'depthv',cdfversion)  ; dpthglobV=squeeze(dpthglobV); %FC for CREG12
[NZ,NY,NX,NT]=size(globU);
 
%- Search for minimum and maximum indexes within the SOURCE data
%corresponding to the DESTINATION area of interest (high res domain +- nbdeg)
index_minii=NX;         index_maxii=-NX;
index_minjj=NY;         index_maxjj=-NY;
for globii=1:NX
  for globjj=1:NY
    % search if point into desired zone
    if ( (Ggphiu(globjj,globii)>=(minlat-nbdeg)) && (Ggphiu(globjj,globii)<=(maxlat+nbdeg)) )    
      if ( (Gglamu(globjj,globii)>=(minlon-nbdeg)) && (Gglamu(globjj,globii)<=(maxlon+nbdeg)) )    
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
reg_glob_latT=Ggphit(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_lonT=Gglamt(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_latU=Ggphiu(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_lonU=Gglamu(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_latV=Ggphiv(index_minjj:index_maxjj,index_minii:index_maxii);
reg_glob_lonV=Gglamv(index_minjj:index_maxjj,index_minii:index_maxii);
clear globii globjj Ggphiu Gglamu Ggphiv Gglamv Ggphit Gglamt


% [3b] -- LOAD DATA 
display(' 3b - LOAD DATA')
time_beg_load=clock;                                                 %jpp
reg_glob_dataSSH=zeros(1,zoomdimj,zoomdimi);
reg_glob_dataT  =zeros(NZ,zoomdimj,zoomdimi);
reg_glob_dataS  =zeros(NZ,zoomdimj,zoomdimi);
reg_glob_dataU  =zeros(NZ,zoomdimj,zoomdimi);
reg_glob_dataV  =zeros(NZ,zoomdimj,zoomdimi);
time_counter=zeros(1,1);

% -- EXTRACT SSH
display([ '      LOAD SSH data for Y=' num2str(year) ' M=' month ])
%SSHfile=([ inputdata '_ORCA025_' num2str(year) month '_grid2D.nc']);
globSSH=f_readnetcdf(datadir,SSHfile,'sossheig',cdfversion);
globSSH=squeeze(globSSH);
% FCH: The first time step is considered for making the initial condition,
% i.e itime=1
% Due to the way JP permutes the variables in f_readnetcdf, the t has
% different index in 2D and 3D field

itime=1;
reg_glob_dataSSH(1,:,:,1)=squeeze(globSSH(itime,index_minjj:index_maxjj,index_minii:index_maxii));
clear globSSH

% -- EXTRACT T
display([ '      LOAD   T data for Y=' num2str(year) ' M=' month ])
%Tfile=([ inputdata '_ORCA025_' num2str(year) month '_gridT.nc']);
globT=f_readnetcdf(datadir,Tfile,'votemper',cdfversion);
globT=squeeze(globT);
reg_glob_dataT(:,:,:,1)=squeeze(globT(:,index_minjj:index_maxjj,index_minii:index_maxii,itime));
clear globT

% -- EXTRACT S
display([ '      LOAD   S data for Y=' num2str(year) ' M=' month ])
%Sfile=([ inputdata '_ORCA025_' num2str(year) month '_gridS.nc']);
globS=f_readnetcdf(datadir,Sfile,'vosaline',cdfversion);
globS=squeeze(globS);
reg_glob_dataS(:,:,:,1)=squeeze(globS(:,index_minjj:index_maxjj,index_minii:index_maxii,itime));
clear globS

% -- EXTRACT U
display([ '      LOAD   U data for Y=' num2str(year) ' M=' month ])
%Ufile=([ inputdata '_ORCA025_' num2str(year) month '_gridU.nc']);
globU=f_readnetcdf(datadir,Ufile,'vozocrtx',cdfversion);
globU=squeeze(globU);
reg_glob_dataU(:,:,:,1)=squeeze(globU(:,index_minjj:index_maxjj,index_minii:index_maxii,itime));
clear globU

% -- EXTRACT V
display([ '      LOAD   V data for Y=' num2str(year) ' M=' month ])
%Vfile=([ inputdata '_ORCA025_' num2str(year) month '_gridV.nc']);
globV=f_readnetcdf(datadir,Vfile,'vomecrty',cdfversion);
globV=squeeze(globV);
reg_glob_dataV(:,:,:,1)=squeeze(globV(:,index_minjj:index_maxjj,index_minii:index_maxii,itime));
clear globV

% Calculate the variable time_counter defined as : days since 1950-01-01
% VERY IMPORTANT FOR THE MODEL, THAT IS HOW NEMO KNOWS WHEN TO
% READ/CHANGE OBC CONDITIONS %% !!!!! NOT APPLIED HERE !!!!!

% FC: Better way for time counter
datestr='01/01/1950';
dateend=[month,'/01/',year];
time_counter(1,1)=daysact(datestr,dateend);

%  
% moduloyr4  =mod(year,4);
% moduloyr100=mod(year,100);
% moduloyr400=mod(year,400);
% if moduloyr400==0
%     julian_day=[1,32,61,92,122,153,183,214,245,275,306,336];
% elseif moduloyr4==0 && moduloyr100~=0
%     julian_day=[1,32,61,92,122,153,183,214,245,275,306,336];
% else
%     julian_day=[1,32,60,91,121,152,182,213,244,274,305,335];
% end
% time_counter(1,1)=(year-1950+1)*365.25+julian_day(str2num(month));
time_end_load=clock;


% -------------------------------------------------------------------------
% [4] Processing scalar fields - Use griddata for coherence with U,V
% -------------------------------------------------------------------------
time_beg_scalar=clock;                                                     %jpp
% [4a] -- Create modified mask for floodnan3_opa for T grid...
display(' 4a - Create mask for flooding')
[mod_mask,count]=f_create_mask(reg_glob_dataT,nbpts);

% [4b] --  Interpolating and flooding from SOURCE to DESTINATION 
display(' 4b - Interpolating and flooding')

% FCH: Assign NAN to overland values of T&S
% if I don't we get wierd values when interpolating horizontally
reg_glob_dataS(reg_glob_dataS==0)=NaN;
reg_glob_dataT(reg_glob_dataT==0)=NaN;
 
[outT]  =f_interp_scalar(reg_glob_dataT,reg_glob_lonT,reg_glob_latT,...
                         coordREG.glamt,coordREG.gphit,...
                         dpthglobV, mdepth,mod_mask);
 
[outS]  =f_interp_scalar(reg_glob_dataS,reg_glob_lonT,reg_glob_latT,...
                         coordREG.glamt,coordREG.gphit,...
                         dpthglobV, mdepth,mod_mask);

[outSSH]=f_interp_scalar(reg_glob_dataSSH,reg_glob_lonT,reg_glob_latT,...
                         coordREG.glamt,coordREG.gphit,...
                         dpthglobV, mdepth,mod_mask);
time_end_scalar=clock;                                                     %jpp

% [4b] ---  Stratification verification and Convection
time_beg_convect=clock;                                                    %jpp
display(' 4b - Stratification and convection #####  REQUIRE OPTIMIZATION  ######')
[nzz,nyy,nxx,ntt]=size(outT);
for ti=1:ntt
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
        
        % -- COPY BOTTOM CELLS DOWNWARD
        if (kmtm>1) && (kmtm<nzz)
          outT(kmtm+1:nzz,jj,ii,ti)=outT(kmtm,jj,ii,ti);
          outS(kmtm+1:nzz,jj,ii,ti)=outS(kmtm,jj,ii,ti);
        end
        
      end
    end
  end
end
time_end_convect=clock;                                                    %jpp
                 

% -------------------------------------------------------------------------
% [5] call to interpolation (and vector rotation) function
% -------------------------------------------------------------------------
time_beg_vector=clock;                                                     %jpp
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
xi=1:zoomdimi;
yi=1:zoomdimj;
axi=1:nxx;
ayi=1:nyy;
display('       WARNING NaNs ARE REPLACED BY 0 FOR U and V INTERPOLATION!!!')
reg_glob_dataU(isnan(reg_glob_dataU))=0;
reg_glob_dataV(isnan(reg_glob_dataV))=0;

%dbstop in f_interpuv.m % debugging mode in f_interpuv
[outU outV]=f_interpuv(reg_glob_dataU,reg_glob_dataV,...
                       reg_glob_lonU ,reg_glob_latU,reg_glob_lonV,reg_glob_latV, ...
                       coordREG.glamu,coordREG.gphiu(ayi,axi),...
                       coordREG.glamv,coordREG.gphiv(ayi,axi),...
                       dpthglobV, mdepth, ...
                       coordGLOB, xi, yi, ...
                       coordREG ,axi,ayi);
time_end_vector=clock;                                                     %jpp

% -----  PROCESS MASKS FOR OUTPUT
maskTout1=tmask(1,ayi,axi);
maskTout =tmask(:,ayi,axi); 
maskUout =umask(:,ayi,axi); 
maskVout =vmask(:,ayi,axi); 

%-------------------------------------------------------------------------
% [A] WRITE NETCDF FILES FOR OBC 
%-------------------------------------------------------------------------
display('WRITE NETCDF FILES')
time_beg_write=clock;                                                      %jpp
fname = ([outputdir,'/','IC_' config_name '_sossheig_' inputdata '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'sossheig',...
                       coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
                       1,time_counter,maskTout1,outSSH);

fname = ([outputdir,'/','IC_' config_name '_votemper_' inputdata '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'votemper',...
                       coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
                       mdepth,time_counter,maskTout,outT);

fname = ([outputdir,'/','IC_' config_name '_vosaline_' inputdata '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'vosaline',...
                       coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
                       time_counter,maskTout,outS);

fname = ([outputdir,'/','IC_' config_name '_vozocrtx_' inputdata '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'vozocrtx',...
                       coordREG.glamu(ayi,axi),coordREG.gphiu(ayi,axi),mdepth,...
                       time_counter,maskUout,outU);

fname = ([outputdir,'/','IC_' config_name '_vomecrty_' inputdata '.nc']);
[~] = f_writenetcdf(fname,cdfversion,'vomecrty',...
                       coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
                       time_counter,maskVout,outV);
                   
time_end_write=clock;                                                      %jpp
time_glob_end=clock;                                                       %jpp
display('*******************************************')                     %jpp
display('*****      SHOW EXECUTION TIMES       *****')                     %jpp
display([ 'OVERALL TIME BEG :' num2str(time_glob_beg)   ])                 %jpp
display([ 'OVERALL TIME END :' num2str(time_glob_end)   ])                 %jpp
display([ '   LOAD TIME BEG :' num2str(time_beg_load)   ])                 %jpp
display([ '   LOAD TIME END :' num2str(time_end_load)   ])                 %jpp
display([ ' SCALAR TIME BEG :' num2str(time_beg_scalar) ])                 %jpp
display([ ' SCALAR TIME END :' num2str(time_end_scalar) ])                 %jpp
display([ 'CONVECT TIME BEG :' num2str(time_beg_convect)])                 %jpp
display([ 'CONVECT TIME END :' num2str(time_end_convect)])                 %jpp
display([ ' VECTOR TIME BEG :' num2str(time_beg_vector) ])                 %jpp
display([ ' VECTOR TIME END :' num2str(time_end_vector) ])                 %jpp
display([ '  WRITE TIME BEG :' num2str(time_beg_write)  ])                 %jpp
display([ '  WRITE TIME END :' num2str(time_end_write)  ])                 %jpp
display('*******************************************')                                                            %jpp

fprintf('Finished\n');