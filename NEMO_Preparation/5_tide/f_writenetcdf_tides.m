%JP Paquin - Jul2013 : Allows all components to be appended in 1 file
%          - Jun2013 : Adapt for netcdf 3.6.2
%                      Older MATLAB UNTESTED
function[isok]=f_writenetcdf_tides(namefile,cdfversion,allvarNC,...
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
case {'K1_amp_elev','K2_amp_elev','M2_amp_elev','N2_amp_elev',...
      'O1_amp_elev','P1_amp_elev','Q1_amp_elev','S2_amp_elev'}
  idVAR(vv) = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'units','m');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);

case {'K1_phi_elev','K2_phi_elev','M2_phi_elev','N2_phi_elev',...
      'O1_phi_elev','P1_phi_elev','Q1_phi_elev','S2_phi_elev'}
  idVAR(vv) = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'units','degree|GMT');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);
   
    
% - "U" VELOCITY    
case {'K1_amp_u','K2_amp_u','M2_amp_u','N2_amp_u',...
      'O1_amp_u','P1_amp_u','Q1_amp_u','S2_amp_u'}
  idVAR(vv)  = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'units','m2/s');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);

case {'K1_phi_u','K2_phi_u','M2_phi_u','N2_phi_u',...
      'O1_phi_u','P1_phi_u','Q1_phi_u','S2_phi_u'}
  idVAR(vv)  = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'units','degree|GMT');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);    
    
    
 % - "V" VELOCITY    
case {'K1_amp_v','K2_amp_v','M2_amp_v','N2_amp_v',...
      'O1_amp_v','P1_amp_v','Q1_amp_v','S2_amp_v'}
  idVAR(vv)  = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'units','m2/s');
    netcdf.putAtt(ncout,idVAR(vv),'missing_value',0);

case {'K1_phi_v','K2_phi_v','M2_phi_v','N2_phi_v',...
      'O1_phi_v','P1_phi_v','Q1_phi_v','S2_phi_v'}
  idVAR(vv) = netcdf.defVar(ncout,allvarNC{vv}      ,'float', [dimidx,dimidy] );
    netcdf.putAtt(ncout,idVAR(vv),'units','degree|GMT');
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

else 
  display('***** OLDER VERSIONS NETCDF NOT TESTED *****')
  isok=0;
  return
end

isok=1;
end