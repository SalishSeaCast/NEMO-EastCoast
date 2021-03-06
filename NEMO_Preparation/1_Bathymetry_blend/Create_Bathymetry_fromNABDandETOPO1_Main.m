% ===================================================================
%
% NAME: Create_Bathymetry_fromNABDandEtopo1_Main.m
%
% Description:
% Script for Creating Bathymetry dataset by combining
% North Atlantic Bathymetry in-situ data and ETOPO1 background data
%
% Namelist Input:
%   *Selected Domain Area: minimum and maximum longitude and latitude
%   *Name of Input NABD (Herman) Data file
%   *Webtide input files including path to files & name of .nod file
%   *Tidal constituents to consider for datum conversion
%   *Etopo file (xyz format)
%   *Destination grid file (default is netcdf file with NEMO format)
%   *Parameters for Optimal interpolation method
%   *Outputdir : Directory where you want to store your output files
%
% Processing Steps:
%   1. Preprocess data
%     1.1. Extract data that are in Selected Domain
%     1.2. Eliminate bad and interpolated data these include:
%          Eliminate lines in file that do not have 14 columns
%          Eliminate lines that do not have proper reference level
%          Eliminate lines that have suspisous unit for depth
%          i.e. depth is greater than 10,000 m
%          or given depth is not negative
%          (default NABD depth is negative down)
%          Eliminate all interpolated data
%     1.3. Separate data that should be corrected with WebTide
%       if datum is MSL or depth>200m do not need correction
%       if datum is LLWLT or depth < 200m data should be corrected
%
%   2. Unify vertical datum using webtide
%
%   3. Grid insitu Bathymetry to destination grid using median value
%   4. Blend Gridded insitu data
%   with Background bathymetry such as ETOPO1
%
% ====================================================================
%
% % Author: Fatemeh Chegini (fatemeh.chegini@dal.ca)
% Date: June 2014
% ====================================================================

%================ BEGIN USER INPUT PARAMETERS ========================

% Select destination domain area
% by giving maximum and minimum longitude and latititude

lonmin=-90;
lonmax=-75;
latmin=73;
latmax=75;

% nput file for NABD Bathymetry dataset
Hermandata='/media/Data/NEMO/Data/nwat21.txt';

% Output directory where output files will be stored
outputdir='/media/Data/NEMO/FC/Bathymetry/output/';

% Path to WebTide data files including .s2c and .nod files
pathWeb='/home/fateme/NEMO/Data/WebTide/WebTide-install/data/arctic9/';

% Name of Webtide .nod file
nodeWeb='arctic9.nod'; %WebTide node file name

% List of tidal constituents to be used
list_compo={'M2' , 'S2', 'N2','K1','O1'};

% Grid destination file
NemoGridFile='/media/Data/NEMO/FC/Capesable_coordinates.nc';

% Etopo file
EtopoFile='/media/Data/NEMO/Data/ETOPO1/Etopo.xyz';

% parameters for OI method
Ls=0.02;
sigobs=1;
dlon=0.2;
dlat=0.2;

% Steps to be processed
Step1=true;
Step2=false;
Step3=false;
Step4=false;

%================ End USER INPUT PARAMETERS==========================

% Other defined parameters , Do not change
fillnumber=999999;
MSLfile=[outputdir,'/MSLData.lld'];
LLWLTfile=[outputdir,'/LLWLTData.lld'];

%% Read input data

% Read destination grid longitude and latitude (Only if Step3 is true)
if Step3
   % DestinationGrid.lon=double(ncread(NemoGridFile,'nav_lon'));
   % DestinationGrid.lat=double(ncread(NemoGridFile,'nav_lat'));
   x=[lonmin:0.1:lonmax];
   y=[latmin:0.1:latmax];
   [X,Y]=meshgrid(x,y);
   DestinationGrid.lon=X;
   DestinationGrid.lat=Y;
end

% Read Etopo data (Only if Step4 is true)
if Step4
    tmp=load(EtopoFile);
    Etopo.lon=tmp(:,1);
    Etopo.lat=tmp(:,2);
    Etopo.depth=tmp(:,3);
    clear tmp
end

%% STEP 1 : PREPROCESS IN-SITU DATA
% This step is carried out by running
% the bash script 1_Preprocess_NABD_data.sh

% First write some input data to a file
if Step1
    system('rm 0_Input.sh')
    fileID=fopen('0_Input.sh','w');
    fprintf(fileID,'%10s\n','#!/bin/sh')
    fprintf(fileID,'%10s\n',['lonmin=',num2str(lonmin)])
    fprintf(fileID,'%10s\n',['latmin=',num2str(latmin)])
    fprintf(fileID,'%10s\n',['lonmax=',num2str(lonmax)])
    fprintf(fileID,'%10s\n',['latmax=',num2str(latmax)])
    fprintf(fileID,'%50s\n',['Hermandata=',Hermandata])
    fprintf(fileID,'%50s\n',['outputdir=',outputdir])
    
%   Run the script
    system('. ./1_Preprocess_NABD_data.sh')
    
end
%% STEP 2: UNIFY VERTICAL DATUM
if Step2
    
    %   Read pre-processed data
    MSLData=load(MSLfile);
    LLWLTData=load(LLWLTfile);
    
    % Convert LLWLT to MSL
    LLWLT2MSL=f_WebTide_datum_conversion ...
        (pathWeb,nodeWeb,LLWLTData,list_compo);
    
    % Add all data together
    tmp=[LLWLT2MSL;MSLData];
    InsituData.lon=tmp(:,1);
    InsituData.lat=tmp(:,2);
    InsituData.depth=tmp(:,3);
    clear tmp
end

%% STEP 3: GRID PROCESSED INSITU DATA TO MODEL GRIDS
if Step3
    GriddedInsituData = f_grid_insitu_bathymetry_2modelgrid ...
        (InsituData,DestinationGrid);
end
%% STEP 4 : COMBINE ETOPO1 and GRIDDED INSITU DATA

% First interploate ETOP1 to model grids
if Step4
    InterpEtopo = griddata(Etopo.lon,Etopo.lat,Etopo.depth, ...
        DestinationGrid.lon(:),DestinationGrid.lat(:));
    InterpEtopo(InterpEtopo>=0)=NaN;
    
    InterpEtopo=reshape(InterpEtopo,size(DestinationGrid.lon));
    
    % Blend ETOPO1 and in-situ data
    
    BlendedData.depth = f_Optimal_interpolation ...
        (DestinationGrid.lon,DestinationGrid.lat, ...
        InterpEtopo,GriddedInsituData,dlon,dlat,Ls,sigobs);
    
    %% WRITE OUTPUTS TO FILES
    BlendedData.lon=DestinationGrid.lon;
    BlendedData.lat=DestinationGrid.lat;
    
    save ([outputdir,'BlendedBathymetry.mat'], 'BlendedData')
end


