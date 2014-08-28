% NAME: interp_tides_Webtide.m
%
% AUTHOR: J.-P. Paquin
%
% DATE: July 2013
%
% REVISIONS: 
%    July13: JP Paquin: - Adapt to netcdf 3.6.2 read, permute and write
%                       - Option for extrapolation on land of input data
%                         to avoid discontinuities in the input fields
%                       - Add comments and details
%
%
% DESCRIPTION: Extract and interpolate tidal components from Webtide 
%
% HYPOTHESES:  
%
% INPUTS : DESTINATION 'NEMO regional' grid: coordinates, bathymetry
%          ASCII tidal components from Webtide: 
%               'K1','K2','M2','N2','O1','P1','Q1','S2'
%          www.bio.gc.ca/science/research-recherche/ocean/webtide/index-eng.php
%
% OUTPUTS: All Tidal components grouped in files :
%               tide_elev.nc, tide_ubar.nc, tide_vbar.nc
%
%
% CALLED PGM & SCRIPTS: f_open_netcdf
%                       f_readnetcdf
%                       f_writenetcdf_tides
%                       f_opa_angle
%                       floodnan3_opa
%			            floodnan4_opa
%
% NOTES:
%
%   ***********************************************************************
%   ******* NETCDF VERSION SUPPORTED 3.6.2 : OLDER VERSION UNTESTED *******
%   *******             SEE KEY 'cdfversion' IN CODE                *******
%   ***********************************************************************
%
%--------------------------------------------------------------------------
clear all
display('***** CALL TO : interp_tides_Webtide.m *****')

%##############################################################
%################## BEGIN USER MODIFICATIONS ##################
%- Specify the requested tidal components:

list_compo={'M2' , 'S2', 'N2' };

% - Specify input directory for NEMO coordinates and Bathymetry
% - Coordinate file must be in the execution directory...
pathcoord='/media/Data/NEMO/FC/NEMO_input'; 
coordfile='Capesable_coordinates_ORCA05km.nc';

%- New bathymetry with j>545 = 0
pathbathy ='/media/Data/NEMO/FC/NEMO_input';
filebathy = 'Capesable_bathymetry_ORCA05km.nc';


% - Specify input directory for WebTide
inputdata='Webtide_sshelf';
pathWeb ='/media/Data/NEMO/Data/WebTide/WebTide-install/data/sshelf';
fileWebnode ='sshelf_ll.nod';
fileWebBathy='sshelf.bat';


% - Input parameters
landextrapolation='yes'; %- Option if data need to be extrapolated on land
                         %  use of floodnan3_opa 
cdfversion='3.6.2'; %- Supported netcdf version: v3.6.2
                    %  Older versions untested

addpath('/media/Data/NEMO/JPP/0_common');
addpath('/media/Data/NEMO/JPP/5_tide');
%##############################################################


%JPP - for showing execution times
[~,ncompo]=size(list_compo);
time_T_beg_interp =zeros(1,6,ncompo);                                      %jpp
time_T_end_interp =zeros(1,6,ncompo);                                      %jpp
time_T_beg_angles =zeros(1,6,ncompo);                                      %jpp
time_T_end_angles =zeros(1,6,ncompo);                                      %jpp
time_T_beg_reste  =zeros(1,6,ncompo);                                      %jpp
time_T_end_reste  =zeros(1,6,ncompo);                                      %jpp


time_T_beg_total=clock;                                                    %jpp

% ------------------------------------------------------------------------
% [1] load coordinates and bathymetry fileng data and necessary constants
% ------------------------------------------------------------------------
% - Load bathymetry
display('1- LOAD BATHYMETRY AND COORDINATES')
if strcmp(cdfversion,'3.6.2')
  tmp_data=f_open_netcdf(pathbathy,filebathy,'Bathymetry');
  bathy=permute(tmp_data,[2,1]); 
  clear tmp_data;
end
mask=bathy;
mask(mask>1)=1;
[NY,NX]=size(bathy);

% - Load coordinates
if strcmp(cdfversion,'3.6.2')
  tmp_data=f_open_netcdf(pathcoord,coordfile,'glamt');   Alont=permute(tmp_data,[2,1]); clear tmp_data;
  tmp_data=f_open_netcdf(pathcoord,coordfile,'gphit');   Alatt=permute(tmp_data,[2,1]); clear tmp_data;
  tmp_data=f_open_netcdf(pathcoord,coordfile,'glamu');   Alonu=permute(tmp_data,[2,1]); clear tmp_data;
  tmp_data=f_open_netcdf(pathcoord,coordfile,'gphiu');   Alatu=permute(tmp_data,[2,1]); clear tmp_data;
  tmp_data=f_open_netcdf(pathcoord,coordfile,'glamv');   Alonv=permute(tmp_data,[2,1]); clear tmp_data;
  tmp_data=f_open_netcdf(pathcoord,coordfile,'gphiv');   Alatv=permute(tmp_data,[2,1]); clear tmp_data;

  % - load coordinates of DESTINATION grid
  Rglamt=f_readnetcdf(pathcoord,coordfile,'glamt',cdfversion);  glamt=squeeze(Rglamt);
  Rglamu=f_readnetcdf(pathcoord,coordfile,'glamu',cdfversion);  glamu=squeeze(Rglamu);
  Rglamv=f_readnetcdf(pathcoord,coordfile,'glamv',cdfversion);  glamv=squeeze(Rglamv);
  Rglamf=f_readnetcdf(pathcoord,coordfile,'glamf',cdfversion);  glamf=squeeze(Rglamf);
  Rgphit=f_readnetcdf(pathcoord,coordfile,'gphit',cdfversion);  gphit=squeeze(Rgphit);
  Rgphiu=f_readnetcdf(pathcoord,coordfile,'gphiu',cdfversion);  gphiu=squeeze(Rgphiu); 
  Rgphiv=f_readnetcdf(pathcoord,coordfile,'gphiv',cdfversion);  gphiv=squeeze(Rgphiv);   
  Rgphif=f_readnetcdf(pathcoord,coordfile,'gphif',cdfversion);  gphif=squeeze(Rgphif); 

  %- Definition of a structure variable containing all DESTINATION grid lats/lons
  coordREG=struct('glamt',glamt,'glamu',glamu,'glamv',glamv,'glamf',glamf, ...
                  'gphit',gphit,'gphiu',gphiu,'gphiv',gphiv,'gphif',gphif  );
  clear Rglamt Rglamu Rglamv Rglamf Rgphit Rgphiu Rgphiv Rgphif glamt glamf gphit gphif

end



% ------------------------------------------------------------------------
% [2] Create land-sea mask for flooding
% ------------------------------------------------------------------------
display('2- CREATE LAND-SEA MASK')
tmpmask=(['flood_mask_' inputdata '.mat']); %- Check if temporary file exist 
time_T_beg_floodM =clock;                                                 %jpp
if ~exist(tmpmask,'file'); 
  nanmask=mask;
  nanmask(nanmask==0)=nan;

  mask0=nanmask.*0;
  mask1=zeros(size(mask));

  mask2=floodnan4_opa(mask0,mask1,3);

  save (tmpmask,'mask2')  % save temporary file
  time_T_end_floodM =clock;                                                %jpp
else
  display(['   LOAD INTERMEDIATE FILE: ' tmpmask])
  load (tmpmask,'mask2') 
  time_T_end_floodM =clock;                                                %jpp
end

mask2(isnan(mask2))=0; %FCH: mask doesn't work for me ?!

% ------------------------------------------------------------------------
% [3] Loop over all required tidal components
% ------------------------------------------------------------------------
display('3- INTERPOLATE TIDAL COMPONENT(S)')

% - [3-1] Load latitudes/ longitudes and bathymetry
display('  3-1 LOAD WEBTIDE LAT/LON & BATHY')
% - Read WebTide Latitudes and Longitudes
latlonweb  =([pathWeb ,'/', fileWebnode]);           % Webtide Lat-Lon
HEADERLINES=0; DELIMITER=' ';
tmpdata = importdata(latlonweb, DELIMITER, HEADERLINES);
lonin=tmpdata(:,2)  ; % longitude
latin=tmpdata(:,3)  ; % latitude
ii=find(lonin>0)  ;lonin(ii) = lonin(ii)-360 ; % Remove lines where long
clear tmpdata                                  % are superior to 0

% - Read Webtide Bathymetry
bathyweb   =([pathWeb ,'/',fileWebBathy]);              % WebTide bathymetry
HEADERLINES=0; DELIMITER=' ';
tmpdata = importdata(bathyweb, DELIMITER, HEADERLINES);
bathyin=tmpdata(:,2)  ; % bathymetry   
clear tmpdata


% - PROCESS TIDAL COMPONENTS - 
dataH(ncompo*2,NY,NX)=NaN;
dataU(ncompo*2,NY,NX)=NaN;
dataV(ncompo*2,NY,NX)=NaN;
allvarsH{ncompo*2,1}=' ';
allvarsU{ncompo*2,1}=' ';
allvarsV{ncompo*2,1}=' ';
nbcompo=1; 
for mycomponents = list_compo
  compo = mycomponents{1}; mycomponents; 

  time_T_beg_interp(:,:,nbcompo)=clock;                                  %jpp
  
  display([ '  3-2 INTERPOLATE COMPONENT : ' compo ])
  tidefile   =([pathWeb '/' compo '.barotropic.v2c']); % velocity
  tideHfile  =([pathWeb '/' compo '.barotropic.s2c']); % elevation
  
  % - Load vlocities amplitude and phase
  HEADERLINES=3; DELIMITER=' ';                       
  tmp= importdata(tidefile, DELIMITER, HEADERLINES);
  tmpdata=tmp.data ; % clear tmp 
  Ua=tmpdata(:,2)  ; % u amplitude
  up=tmpdata(:,3)  ; % u phase
  Va=tmpdata(:,4)  ; % v amplitude
  vp=tmpdata(:,5)  ; % v phase
  clear tmpdata
  
  % - Load amplitude and phase
%  ha=[];hp=[];
  HEADERLINES=3; DELIMITER=' ';      
  tmp= importdata(tideHfile, DELIMITER, HEADERLINES);
  tmpdata=tmp.data ; % clear tmp 
  ha=tmpdata(:,2)  ; % H amplitude
  hp=tmpdata(:,3)  ; % H phase
  clear tmpdata
  
  % ORIGINAL CALCULATION BY Li Zhai
  % - Multiply Depth average velocity by bathymetry
  Ua=Ua.*bathyin  ; % u**2/s
  Va=Va.*bathyin  ; % v**2/s

  % - Transfer to sin and cos part of U and V
  u1=Ua.*cosd(up);
  u2=Ua.*sind(up);
  v1=Va.*cosd(vp);
  v2=Va.*sind(vp);
  h1=ha.*cosd(hp);
  h2=ha.*sind(hp);

  clear Ua Va

  % Interpolate to DESTINATION grid
  display(['    INTERPOLATING COMPONENT : ' compo ])
  
  %- Interpolate Webtide bathymetry to model grid...
  dumba=bathyin;         %dum(dum==0)=nan;
  dumba1 = griddata(lonin,latin,dumba,Alont,Alatt,'linear');
  bathyWT1=dumba1.*mask;
  bathyWT2=floodnan3_opa(bathyWT1,mask2,3);
  clear dum dum1

  dum=u1;         dum(dum==0)=nan;
  dum1 = griddata(lonin,latin,dum,Alonv,Alatv,'linear');
  utt1=dum1.*mask;
  utt1=floodnan3_opa(utt1,mask2,3);
  clear dum dum1 u1

  dum=u2;         dum(dum==0)=nan;
  dum1 = griddata(lonin,latin,dum,Alonv,Alatv,'linear');
  utt2=dum1.*mask;
  utt2=floodnan3_opa(utt2,mask2,3);
  clear dum dum1 u2

  dum=v1;         dum(dum==0)=nan;
  dum1 = griddata(lonin,latin,dum,Alonv,Alatv,'linear');
  vtt1=dum1.*mask;
  vtt1=floodnan3_opa(vtt1,mask2,3);
  clear dum dum1 v1

  dum=v2;         dum(dum==0)=nan;
  dum1 = griddata(lonin,latin,dum,Alonv,Alatv,'linear');
  vtt2=dum1.*mask;
  vtt2=floodnan3_opa(vtt2,mask2,3);
  clear dum dum1 v2

  dum=h1;         dum(dum==0)=nan;
  dum1 = griddata(lonin,latin,dum,Alont,Alatt,'linear');
  htt1=dum1.*mask;
  htt1=floodnan3_opa(htt1,mask2,3);
  clear dum dum1 h1

  dum=h2;         dum(dum==0)=nan;
  dum1 = griddata(lonin,latin,dum,Alont,Alatt,'linear');
  htt2=dum1.*mask;
  htt2=floodnan3_opa(htt2,mask2,3);
  clear dum dum1 h2

  time_T_end_interp(:,:,nbcompo)=clock;                                  %jpp

  
  % get the angle
  display(['    COMPUTING ANGLES : ' compo ])
  time_T_beg_angles(:,:,nbcompo)=clock;                                  %jpp
  %[gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf] = opa_angle(coordfile);
  xi=1:size(Alatt,2);
  yi=1:size(Alatt,1);
  [gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf] = f_opa_angle(coordREG,xi,yi);
  time_T_end_angles(:,:,nbcompo)=clock;                                  %jpp

  
  time_T_beg_reste(:,:,nbcompo)=clock;                                   %jpp

  % rotate the velocity relative to the ji-axis in NEMO
  u1n=utt1  .* gcosu+vtt1 .* gsinu;
  u2n=utt2  .* gcosu+vtt2 .* gsinu;
  v1n=-utt1 .* gsinv+vtt1 .* gcosv;
  v2n=-utt2 .* gsinv+vtt2 .* gcosv; 
  clear utt1 utt2 vtt1 vtt2
  
  
  % Interpolate u back to u point 
  dum=u1n;
  dum1 = griddata(Alonv,Alatv,dum,Alonu,Alatu,'linear');
  u1n=dum1; clear dum dum1
  dum=u2n;
  dum1 = griddata(Alonv,Alatv,dum,Alonu,Alatu,'linear');
  u2n=dum1; clear dum dum1

 
  ha=abs(htt1-htt2.*sqrt(-1));
  hg=-angle(htt1-htt2.*sqrt(-1))*180/pi;
  ua=abs(u1n-u2n.*sqrt(-1));
  ug=-angle(u1n-u2n.*sqrt(-1))*180/pi;
  va=abs(v1n-v2n.*sqrt(-1));
  vg=-angle(v1n-v2n.*sqrt(-1))*180/pi;

  dataH(nbcompo*2-1,:,:)=ha;
  dataH(nbcompo*2  ,:,:)=hg;
  allvarsH{nbcompo*2-1,1}=[ compo '_amp_elev' ];
  allvarsH{nbcompo*2,1}  =[ compo '_phi_elev' ];

  dataU(nbcompo*2-1,:,:)=ua;
  dataU(nbcompo*2  ,:,:)=ug;
  allvarsU{nbcompo*2-1}=[ compo '_amp_u' ];
  allvarsU{nbcompo*2}  =[ compo '_phi_u' ];

  dataV(nbcompo*2-1,:,:)=va;
  dataV(nbcompo*2  ,:,:)=vg;
  allvarsV{nbcompo*2-1}=[ compo '_amp_v' ];
  allvarsV{nbcompo*2}  =[ compo '_phi_v' ];
  
  time_T_end_reste(:,:,nbcompo)=clock;                                  %jpp

  nbcompo=nbcompo+1;
end
% replicate nan points 
dataH(isnan(dataH))=0;
dataU(isnan(dataU))=0;
dataV(isnan(dataV))=0;



% - Write Netcdf file
display('*****  REMOVE EXISTING FILES  *****')
if exist('tide_elev.nc','file') ; delete('tide_elev.nc') ; end
if exist('tide_ubar.nc','file') ; delete('tide_ubar.nc') ; end
if exist('tide_vbar.nc','file') ; delete('tide_vbar.nc') ; end

% - Write tide ssh file
[~] = f_writenetcdf_tides('tide_elev.nc',cdfversion,allvarsH,Alont,Alatt,dataH);
% - Write tide u file
[~] = f_writenetcdf_tides('tide_ubar.nc',cdfversion,allvarsU,Alonu,Alatu,dataU);
% - Write tide v file
[~] = f_writenetcdf_tides('tide_vbar.nc',cdfversion,allvarsV,Alonv,Alatv,dataV);
  
time_T_end_total=clock;                                                    %jpp

  
  % JPP display execution times
display('#####################################################')
display('###  DISPLAY EXECUTION TIMES')
display([ 'Total beg:' num2str(time_T_beg_total) ])                              %jpp
display([ 'Total end:' num2str(time_T_end_total) ])                              %jpp
display([ 'Flood beg:' num2str(time_T_beg_floodM) ])                             %jpp
display([ 'Flood beg:' num2str(time_T_end_floodM) ])                             %jpp
for ii=1:ncompo                                                                  %jpp
display([ 'Interp ' list_compo{ii} ' beg:' num2str(time_T_beg_interp(:,:,ii)) ]) %jpp
display([ 'Interp ' list_compo{ii} ' end:' num2str(time_T_end_interp(:,:,ii)) ]) %jpp
display([ 'Angles ' list_compo{ii} ' beg:' num2str(time_T_beg_angles(:,:,ii)) ]) %jpp
display([ 'Angles ' list_compo{ii} ' end:' num2str(time_T_end_angles(:,:,ii)) ]) %jpp
display([ 'Reste  ' list_compo{ii} ' beg:' num2str(time_T_beg_reste(:,:,ii)) ])  %jpp
display([ 'Reste  ' list_compo{ii} ' end:' num2str(time_T_end_reste(:,:,ii)) ])  %jpp
end                                                                              %jpp
display('#####################################################')
