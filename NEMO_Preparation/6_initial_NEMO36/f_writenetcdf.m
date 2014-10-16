%JP Paquin - Jun2013 : Adapt for netcdf 3.6.2
%                      Older MATLAB UNTESTED
function[isok]=f_writenetcdf(namefile,cdfversion,varNC,...
                             Alon,Alat,mdepth,time_counter,mask,datain)
isok=0;
[NZ,NY,NX,NT]=size(datain);
if strcmp(cdfversion,'3.6.2')
    
    
ncout=netcdf.create(namefile,'NC_write');
% DEFINE DIMENSIONS AND ATTRIBUTES
  dimidx = netcdf.defDim(ncout,'x',NX);
  dimidy = netcdf.defDim(ncout,'y',NY); 
  dimidz = netcdf.defDim(ncout,'z',NZ);
  dimidt = netcdf.defDim(ncout,'time_counter',netcdf.getConstant('NC_UNLIMITED'));

  % Global attributes
  globatt=netcdf.getConstant('NC_GLOBAL');
  netcdf.putAtt(ncout,globatt,'Description','Initial Condition From ANNA')
  netcdf.putAtt(ncout,globatt,'Author','Fatemeh Chegini')
  netcdf.putAtt(ncout,globatt,'Date',date)
  netcdf.putAtt(ncout,globatt,'Convention','GDT 1.2')
  netcdf.putAtt(ncout,globatt,'file_name',namefile)
  
  idnavlon = netcdf.defVar(ncout,'nav_lon'     ,'float' ,[dimidx,dimidy]);
    netcdf.putAtt(ncout,idnavlon,'units','degrees_east');
    netcdf.putAtt(ncout,idnavlon,'valid_min', min(min(Alon)) );
    netcdf.putAtt(ncout,idnavlon,'valid_max', max(max(Alon)) );
    netcdf.putAtt(ncout,idnavlon,'long_name','Longitude at t-point');
    
  idnavlat = netcdf.defVar(ncout,'nav_lat'     ,'float' ,[dimidx,dimidy]);
    netcdf.putAtt(ncout,idnavlat,'units','degrees_north');
    netcdf.putAtt(ncout,idnavlat,'valid_min', min(min(Alat)) );
    netcdf.putAtt(ncout,idnavlat,'valid_max', max(max(Alat)) );
    netcdf.putAtt(ncout,idnavlat,'long_name','Latitude at t-point');

switch varNC 
case {'votemper','vosaline'}    
  iddepth  = netcdf.defVar(ncout,'deptht'      ,'double', dimidz );
    netcdf.putAtt(ncout,iddepth,'units','model_levels');
    netcdf.putAtt(ncout,iddepth,'valid_min',1);
    netcdf.putAtt(ncout,iddepth,'valid_max',6000);
    netcdf.putAtt(ncout,iddepth,'long_name','Model levels');
case {'vozocrtx'}
  iddepth  = netcdf.defVar(ncout,'depthu'      ,'double', dimidz );
    netcdf.putAtt(ncout,iddepth,'units','model_levels');
    netcdf.putAtt(ncout,iddepth,'valid_min',1);
    netcdf.putAtt(ncout,iddepth,'valid_max',6000);
    netcdf.putAtt(ncout,iddepth,'long_name','Model levels');
case {'vomecrty'}
  iddepth  = netcdf.defVar(ncout,'depthv'      ,'double', dimidz );
    netcdf.putAtt(ncout,iddepth,'units','model_levels');
    netcdf.putAtt(ncout,iddepth,'valid_min',1);
    netcdf.putAtt(ncout,iddepth,'valid_max',6000);
    netcdf.putAtt(ncout,iddepth,'long_name','Model levels');
case {'sossheig'}
  iddepth  = netcdf.defVar(ncout,'deptht'      ,'double', dimidz );
    netcdf.putAtt(ncout,iddepth,'units','model_levels');
    netcdf.putAtt(ncout,iddepth,'valid_min',1);
    netcdf.putAtt(ncout,iddepth,'valid_max',1);
    netcdf.putAtt(ncout,iddepth,'long_name','Model levels');
end
        
  idtime   = netcdf.defVar(ncout,'time_counter','double', dimidt );
    netcdf.putAtt(ncout,idtime,'units','days since 1950-01-01 00:00:00');
    netcdf.putAtt(ncout,idtime,'calendar','gregorian');
    netcdf.putAtt(ncout,idtime,'title','Time');
    netcdf.putAtt(ncout,idtime,'long_name','Time axis');
    netcdf.putAtt(ncout,idtime,'time_origin','1950-01-15 00:00:00');

switch varNC
case {'votemper'}
  idmask  = netcdf.defVar(ncout,'tmask'       ,'float' ,[dimidx,dimidy,dimidz]);
  idtemp = netcdf.defVar(ncout,'votemper'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','deg_C');
    netcdf.putAtt(ncout,idtemp,'valid_min',-2.00);
    netcdf.putAtt(ncout,idtemp,'valid_max',40.0);
    netcdf.putAtt(ncout,idtemp,'long_name','votemper');
case {'vosaline'}
  idmask  = netcdf.defVar(ncout,'tmask'       ,'float' ,[dimidx,dimidy,dimidz]);
  idtemp = netcdf.defVar(ncout,'vosaline'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','PSU');
    netcdf.putAtt(ncout,idtemp,'valid_min',0. );
    netcdf.putAtt(ncout,idtemp,'valid_max',50.);
    netcdf.putAtt(ncout,idtemp,'long_name','vosaline');  
case {'sossheig'}
  idmask  = netcdf.defVar(ncout,'tmask'       ,'float' ,[dimidx,dimidy,dimidz]);
  idtemp = netcdf.defVar(ncout,'sossheig'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','m');
    netcdf.putAtt(ncout,idtemp,'valid_min',-100);
    netcdf.putAtt(ncout,idtemp,'valid_max',100);
    netcdf.putAtt(ncout,idtemp,'long_name','Sea surface height');
case {'vozocrtx'} % U
  idmask  = netcdf.defVar(ncout,'umask'       ,'float' ,[dimidx,dimidy,dimidz]);
  idtemp = netcdf.defVar(ncout,'vozocrtx'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','m s-1');
    netcdf.putAtt(ncout,idtemp,'valid_min',-100);
    netcdf.putAtt(ncout,idtemp,'valid_max',100);
    netcdf.putAtt(ncout,idtemp,'long_name','Zonal velocity');
case {'vomecrty'} % V
  idmask  = netcdf.defVar(ncout,'vmask'       ,'float' ,[dimidx,dimidy,dimidz]);
  idtemp = netcdf.defVar(ncout,'vomecrty'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','m s-1');
    netcdf.putAtt(ncout,idtemp,'valid_min',-100);
    netcdf.putAtt(ncout,idtemp,'valid_max',100);
    netcdf.putAtt(ncout,idtemp,'long_name','Meridional velocity');
end
    
netcdf.endDef(ncout);
% PUT VARIABLES
netcdf.putVar(ncout,idnavlon,permute(Alon,[2,1]))
netcdf.putVar(ncout,idnavlat,permute(Alat,[2,1]))
netcdf.putVar(ncout,iddepth ,mdepth)
netcdf.putVar(ncout,idtime  ,0,NT,time_counter)
netcdf.putVar(ncout,idmask ,permute(mask,[3,2,1]))
%netcdf.putVar(ncout,idtemp,permute(datain,[4,3,2,1]))
netcdf.putVar(ncout,idtemp,permute(datain,[3,2,1,4]))

netcdf.close(ncout)

else 
  display('***** OLDER VERSIONS NETCDF NOT TESTED *****')
  isok=0;
  return
end

isok=1;
end