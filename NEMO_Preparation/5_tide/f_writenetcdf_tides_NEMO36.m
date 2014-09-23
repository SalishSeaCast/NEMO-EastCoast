%JP Paquin - Jul2013 : Allows all components to be appended in 1 file
%          - Jun2013 : Adapt for netcdf 3.6.2
%                      Older MATLAB UNTESTED
function[isok]=f_writenetcdf_tides_NEMO36(namefile,cdfversion,allvarNC,...
                             Alon,Alat,datain)
[nbvar,NY,NX]=size(datain);
idVAR=zeros(nbvar);

if strcmp(cdfversion,'3.6.2')

ncout=netcdf.create(namefile,'WRITE');
  
% DEFINE DIMENSIONS AND ATTRIBUTES
dimidx = netcdf.defDim(ncout,'x',NX);
dimidy = netcdf.defDim(ncout,'y',NY); 

idnavlon = netcdf.defVar(ncout,'nav_lon'     ,'float' ,[dimidx,dimidy]);
  netcdf.putAtt(ncout,idnavlon,'units','degrees_east');
idnavlat = netcdf.defVar(ncout,'nav_lat'     ,'float' ,[dimidx,dimidy]);
  netcdf.putAtt(ncout,idnavlat,'units','degrees_north');

for vv=1:nbvar
switch allvarNC{vv}
% - ELEVATION 
case {'K1_z1','K2_z1','M2_z1','N2_z1',...
      'O1_z1','P1_z1','Q1_z1','S2_z1','L2_z1'}
  idVAR(vv) = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'longname','tidal elevation: cosine');
    netcdf.putAtt(ncout,idVAR(vv),'units','m');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);

case {'K1_z2','K2_z2','M2_z2','N2_z2',...
      'O1_z2','P1_z2','Q1_z2','S2_z2','L2_z2'}
  idVAR(vv) = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'longname','tidal elevation: sin');
    netcdf.putAtt(ncout,idVAR(vv),'units','m');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);
   
    
% - "U" VELOCITY    
case {'K1_u1','K2_u1','M2_u1','N2_u1',...
      'O1_u1','P1_u1','Q1_u1','S2_u1','L2_u1'}
  idVAR(vv)  = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'longname','tidal x-velocity: cosine');
    netcdf.putAtt(ncout,idVAR(vv),'units','m/s');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);

case {'K1_u2','K2_u2','M2_u2','N2_u2',...
      'O1_u2','P1_u2','Q1_u2','S2_u2','L2_u2'}
  idVAR(vv)  = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'longname','tidal x-velocity: sine');
    netcdf.putAtt(ncout,idVAR(vv),'units','m/s');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);    
    
    
 % - "V" VELOCITY    
case {'K1_v1','K2_v1','M2_v1','N2_v1',...
      'O1_v1','P1_v1','Q1_v1','S2_v1','L2_v1'}
  idVAR(vv)  = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'longname','tidal y-velocity: cosine');
    netcdf.putAtt(ncout,idVAR(vv),'units','m/s');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);

case {'K1_v2','K2_v2','M2_v2','N2_v2',...
      'O1_v2','P1_v2','Q1_v2','S2_v2','L2_v2'}
  idVAR(vv) = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'longname','tidal y-velocity: sine');
    netcdf.putAtt(ncout,idVAR(vv),'units','m/s');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);         
end
end
netcdf.endDef(ncout);

% PUT VARIABLES
netcdf.putVar(ncout,idnavlon,permute(Alon,[2,1]))
netcdf.putVar(ncout,idnavlat,permute(Alat,[2,1]))
for vv=1:nbvar
  tmpdata=squeeze(datain(vv,:,:));
  netcdf.putVar(ncout,idVAR(vv),permute(tmpdata,[2,1]))
end

netcdf.close(ncout)

%% Write 1D tide files specific for NEMO3.6




else 
  display('***** OLDER VERSIONS NETCDF NOT TESTED *****')
  isok=0;
  return
end

isok=1;
end