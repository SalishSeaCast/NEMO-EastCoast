% NAME: gen_obc_3D.m
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
%            Sep14- Fatemeh Chegini (FCH)
%                   NEMO3.6 can read 2D and 3D boundaries separately
%                   This is useful for regional modeling where you want to
%                   give 2D B.C with a higher frequency compared to 3D
%                   Therefore, I've separated the creation of these B.C.s
%
%            Sep14- FCH: Generalized to read any nemo input file
%
%            Oct14- FCH: Add preparation of 3D for NEMO3.6
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
close all
display('***** CALL TO : gen_obc.m *****')
addpath('/media/Data/NEMO/JPP/0_common/m_map')


%##############################################################
%################## BEGIN USER MODIFICATIONS ##################
% - OUPUT directory
dirout='/media/Data/NEMO/FC/NEMO_OUT_PREP/OBC/ANNA/SCOTIAN_RUN2010_NOTIDES/3D_NEMO110/';

% - DESTINATION grid files required
dirbathy='/media/Data/NEMO/FC/NEMO_input/Pegasus/';
filebathy='Capesable_bathymetry_NWATL_R.nc';
dircoord='/media/Data/NEMO/FC/NEMO_input/Pegasus/';
filecoord='Capesable_coordinates_R.nc';
dirmesh ='/media/Data/NEMO/FC/NEMO_input/Pegasus/';
filemesh='CPSBL_mesh_mask_R.nc';
config_name='CPSBL_NEMO110';

% - SOURCE grid file names and related variables
inputdata='ANNA';
dirin=('/media/Data/NEMO/Anna/Anna_output/SCOTIAN_RUN2010_NOTIDES/');
Tfile='SCOTIAN_T_2010_Bcl_notides_Extracted_CapeSable.nc';
Ufile='SCOTIAN_U_2010_Bcl_notides_Extracted_CapeSable.nc';
Vfile='SCOTIAN_V_2010_Bcl_notides_Extracted_CapeSable.nc';
Uvar='vozocrtx';
Vvar='vomecrty';
Tvar='votemper';
Svar='vosaline';

% - SOURCE bathymetry file
% - FCH: required to make transport consistent
dirbathyG='/media/Data/NEMO/Anna/Anna_output';
filebathyG='bathy_meter_Scotian_NEW_ANNA_CPSBL_Extracted.nc';

% - SOURCE coordinate file
% OPTION 1 : coordinate file for SOURCE grid
dircoordG='/media/Data/NEMO/Anna/Anna_output/';
filecoordG='coordinates_SCOTIAN_NEW_ANNA_CPSBL_Extracted.nc';

% OPTION 2 : load latitudes and longitudes form files of U,V (T and F)

% - Input parameters
cdfversion='3.6.2'; %- Supported netcdf version: v3.6.2
%  Older versions untested

concepts_version='3.6'; % Switch for east and north boundary inversion
% Version 1.1.0           : inversion required
% Version 2.0.0 and above : not required

nemoversion='1.0';         % Switch for how to write boundary files
% Version 3.6 and 3.4 write fields as 1D arrays
% older versions: write fileds as 2D arrays


% FCH: OBC according to nambdy_index in namelist for NEMO3.6
% Note: the order is important, e.g obc.name(1),obc.nbdyind(1),... should
% correspond to first boundary.
if strcmp(nemoversion,'3.6')
    obc.name={'east','west','north','south'}; % - List of open boundaries in the domain
    obc.nbdyind=[384 2 249 2]; %-index of obc according to namelist
    obc.nbdybeg=[2, 2 , 2 ,2]; %- beg
    obc.nbdyend=[250,250,135,385]; %end
else
    obc.name={'east','west','north','south'}; % - List of open boundaries in the domain
end

%Start date of input file
if ~strcmp(nemoversion,'3.6')
    yearstr='2010';
    monthstr='01';
    daystr='01';
end

nbdeg=0.1; %- number of degrees around the DESTINATION grid required for
%  data extraction from the SOURCE grid to avoid (as much as
%  possible) memory problems

nbptsbound=10; %- number of points in the boundary


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
end

% read information from the mesh mask (mdepth required)
[tmask,umask,vmask,fmask,mdepth,gdepw,dzz,~,~] ...
    =f_readmeshmask([dirmesh '/' filemesh],cdfversion,nemoversion);

% determine full cells (mlevel) by checking the bathy_meter file
[nz,ny,nx]=size(tmask);
mlevel=zeros(ny,nx);
for k=1:nz-1
    mlevel(Bathy>gdepw(k)) = k+1;
end
clear nx ny nz


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

%% -------------------------------------------------------------------------------
% SOURCE DATA
% - load coordinates information of SOURCE grid

% - load coordinates information of SOURCE grid
if ~exist([ dircoordG '/' filecoordG ],'file');
    display(' GLOBAL COORDINATE FILE NOT PROVIDED - READS T,U,V,F files')
    
    Ffile=Tfile; %- do not have f coord for ANNA's output
    
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

% - load bathymetry file
BathyG=f_readnetcdf(dirbathyG,filebathyG,'Bathymetry',cdfversion); BathyG=squeeze(BathyG);

%- Definition of a structure variable containing all SOURCE grid lats/lons
coordGLOB=struct('glamt',Gglamt,'glamu',Gglamu,'glamv',Gglamv,'glamf',Gglamf, ...
    'gphit',Ggphit,'gphiu',Ggphiu,'gphiv',Ggphiv,'gphif',Ggphif  );
clear Gglamf Ggphif Gglamu Gglamv Gglamf Ggphit Ggphiu Ggphiv Ggphif

%% READ SOURCE DATA

% -- READ U3D
display('LOAD U data ')
globU=f_readnetcdf(dirin,Ufile,Uvar,cdfversion);

% -- READ V3D
display('LOAD V data ')
globV=f_readnetcdf(dirin,Vfile,Vvar,cdfversion);

% -- READ T & S
display('LOAD T&S data ')
globT=f_readnetcdf(dirin,Tfile,Tvar,cdfversion);
globS=f_readnetcdf(dirin,Tfile,Svar,cdfversion);

% -- READ depthu,depthv,deptht
dpthglobU=f_readnetcdf(dirin,Ufile,'depthu',cdfversion)  ; dpthglobU=squeeze(dpthglobU);
dpthglobV=f_readnetcdf(dirin,Vfile,'depthv',cdfversion)  ; dpthglobV=squeeze(dpthglobV);
dpthglobT=f_readnetcdf(dirin,Tfile,'deptht',cdfversion)  ; dpthglobT=squeeze(dpthglobT);

[NZ,NY,NX,NTIME]=size(globT);

% -- READ TIME_COUNTER
% Calculate the variable time_counter defined as : days since 1950-01-01
% VERY IMPORTANT FOR THE MODEL, THAT IS HOW NEMO KNOWS WHEN TO
% READ/CHANGE OBC CONDITIONS
% FCH NOTE: In NEMO3.6 the time_counter unit is seconds since stdate
% where stdate is indicated in files name
% Therefore no need to change time_counter to days since 1950-01-01
timecounterseconds=f_readnetcdf(dirin,Tfile,'time_counter',cdfversion);
if strcmp(nemoversion,'3.6')
    time_units=ncreadatt([dirin,'/',Tfile],'time_counter','units');
    time_counter=timecounterseconds;
else
    dateref='01/01/1950';
    datestr=[monthstr,'/',daystr,'/',yearstr];
    datediff=daysact(dateref,datestr);
    timecounterdays=timecounterseconds/24/3600;
    time_counter=timecounterdays+datediff;
end

%% -------------------------------------------------------------------------
% MAIN LOOP OVER THE USER DEFINED OPEN BOUNDARIES
%global obcname nbdyind nbdybeg nbdyend ;
nobc=length(obc.name);
for iobc=1:nobc
    obcname = obc.name{iobc};
    display( ['Processing boundary: ' ,obcname]) 
    if strcmp(nemoversion,'3.6')
        nbdyind=obc.nbdyind(iobc);
        nbdybeg=obc.nbdybeg(iobc);
        nbdyend=obc.nbdyend(iobc);
    end
    
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
    
    if strcmp(nemoversion,'3.6')
        [corners] = f_define_corners_NEMO36(nbdyind,nbdybeg,nbdyend,obcname,nbptsbound);
    else
        [RdomY,RdomX]=size(coordREG.glamt);
        [corners]=f_define_corners(RdomY,RdomX,obcname,nbptsbound);
    end
    
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
    % [3a] -- FIND THE MIN/MAX LATITUDES AND LONGITUDES
    %         TO EXTRACT ON THE LARGE SCALE DATA
    display(' 3a - IDENTIFY ZONE IN SOURCE GRID')
    
    %- Search for minimum and maximum indexes within the SOURCE data
    
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
    
    reg_glob_latT=coordGLOB.gphit(index_minjj:index_maxjj,index_minii:index_maxii);
    reg_glob_lonT=coordGLOB.glamt(index_minjj:index_maxjj,index_minii:index_maxii);
    reg_glob_latU=coordGLOB.gphiu(index_minjj:index_maxjj,index_minii:index_maxii);
    reg_glob_lonU=coordGLOB.glamu(index_minjj:index_maxjj,index_minii:index_maxii);
    reg_glob_latV=coordGLOB.gphiv(index_minjj:index_maxjj,index_minii:index_maxii);
    reg_glob_lonV=coordGLOB.glamv(index_minjj:index_maxjj,index_minii:index_maxii);
    
    % Plot to make check extracted region
    figure()
    plot([minlon,maxlon,maxlon,minlon,minlon],[minlat,minlat,maxlat,maxlat,minlat],'ro-')
    hold on; plot(reg_glob_lonT(:),reg_glob_latT(:),'.')
    
    % [3b] --  EXTRACT SOURCE DATA
    display(' 3b - EXTRACTING SOURCE DATA')
    time_genUV_beg_load=clock; %jpp
    
    reg_glob_dataT=squeeze(globT(:,index_minjj:index_maxjj,index_minii:index_maxii,:));
    reg_glob_dataS=squeeze(globS(:,index_minjj:index_maxjj,index_minii:index_maxii,:));
    reg_glob_dataU=squeeze(globU(:,index_minjj:index_maxjj,index_minii:index_maxii,:));
    reg_glob_dataV=squeeze(globV(:,index_minjj:index_maxjj,index_minii:index_maxii,:));
    
    reg_glob_dataT(reg_glob_dataT==0)=NaN;
    reg_glob_dataS(reg_glob_dataS==0)=NaN;
    
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
    
    [outT]  =f_interp_scalar(reg_glob_dataT,reg_glob_lonT,reg_glob_latT,...
        coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
        dpthglobT, mdepth);
    [outS]  =f_interp_scalar(reg_glob_dataS,reg_glob_lonT,reg_glob_latT,...
        coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
        dpthglobT, mdepth);
    
    time_end_scalar=clock;                                                     %jpp
    
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
    maskUout =umask(:,ayi,axi); %**
    maskVout =vmask(:,ayi,axi); %**
    
    %-------------------------------------------------------------------------
    % [A] WRITE NETCDF FILES FOR OBC
    %-------------------------------------------------------------------------
    display('WRITE NETCDF FILES')
    time_genUV_beg_write=clock;                                               %jpp
    %dbstop in f_writenetcdf.m % debugging mode in f_interpuv
    if strcmp(nemoversion,'3.6')
        fname = ([dirout,'/','obc_' obcname '_votemper_' config_name '.nc']);
        [~] = f_writenetcdf_3D(fname,cdfversion,'votemper',...
            coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
            mdepth,time_counter,maskTout,outT,obcname);
        
        fname = ([dirout,'/','obc_' obcname '_vosaline_' config_name '.nc']);
        [~] = f_writenetcdf_3D(fname,cdfversion,'vosaline',...
            coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
            time_counter,maskTout,outS,obcname);
        
        fname = ([dirout,'/','obc_' obcname '_vozocrtx_' config_name '.nc']);
        [~] = f_writenetcdf_3D(fname,cdfversion,'vozocrtx',...
            coordREG.glamu(ayi,axi),coordREG.gphiu(ayi,axi),mdepth,...
            time_counter,maskUout,outU,obcname);
        
        fname = ([dirout,'/','obc_' obcname '_vomecrty_' config_name '.nc']);
        [~] = f_writenetcdf_3D(fname,cdfversion,'vomecrty',...
            coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
            time_counter,maskVout,outV,obcname);
    else
        fname = ([dirout,'/','obc_' obcname '_votemper_' config_name '.nc']);
        [~] = f_writenetcdf(fname,cdfversion,'votemper',...
            coordREG.glamt(ayi,axi),coordREG.gphit(ayi,axi),...
            mdepth,time_counter,maskTout,outT,obcname,concepts_version);
        
        fname = ([dirout,'/','obc_' obcname '_vosaline_' config_name '.nc']);
        [~] = f_writenetcdf(fname,cdfversion,'vosaline',...
            coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
            time_counter,maskTout,outS,obcname,concepts_version);
        
        fname = ([dirout,'/','obc_' obcname '_vozocrtx_'  config_name '.nc']);
        [~] = f_writenetcdf(fname,cdfversion,'vozocrtx',...
            coordREG.glamu(ayi,axi),coordREG.gphiu(ayi,axi),mdepth,...
            time_counter,maskUout,outU,obcname,concepts_version);
        
        fname = ([dirout,'/','obc_' obcname '_vomecrty_' config_name '.nc']);
        [~] = f_writenetcdf(fname,cdfversion,'vomecrty',...
            coordREG.glamv(ayi,axi),coordREG.gphiv(ayi,axi),mdepth,...
            time_counter,maskVout,outV,obcname,concepts_version);
    end
    
    
     time_genUV_end_write=clock;  %jpp
    
    clear reg_glob_dataSSH reg_glob_dataU reg_glob_dataV
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
display([ ' INTERP TIME BEG :' num2str(time_genUV_beg_interp)])            %jpp
display([ ' INTERP TIME END :' num2str(time_genUV_end_interp)])            %jpp
display([ '  WRITE TIME BEG :' num2str(time_genUV_beg_write) ])            %jpp
display([ '  WRITE TIME END :' num2str(time_genUV_end_write) ])            %jpp
display('*******************************************')                                                            %jpp

fprintf('Finished\n');
