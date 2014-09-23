% Program for writing OBC files
% JP Paquin - Jun2013 : Adapt for netcdf 3.6.2
%                      Older MATLAB UNTESTED
%
% Fatemeh Chegini - Sep2014 : Adapt for NEMO 3.6

function f_writenetcdf_2D(namefile,varNC,...
                             Alon,Alat,time_counter,time_units,datain,...
                             obc_name)

datain=squeeze(datain);
[NX,NT]=size(datain);
Alon=squeeze(Alon);
Alat=squeeze(Alat);

    
ncout=netcdf.create(namefile,'NC_write');
% DEFINE DIMENSIONS AND ATTRIBUTES
  dimidx = netcdf.defDim(ncout,'x',NX);
  dimidy = netcdf.defDim(ncout,'y',1); 
  dimidt = netcdf.defDim(ncout,'time_counter',netcdf.getConstant('NC_UNLIMITED'));

  % Global attributes
  globatt=netcdf.getConstant('NC_GLOBAL');
  netcdf.putAtt(ncout,globatt,'Description','OBC Source from ANNA')
  netcdf.putAtt(ncout,globatt,'Author','Fatemeh Chegini')
  netcdf.putAtt(ncout,globatt,'Date',date)
  netcdf.putAtt(ncout,globatt,'Convention','GDT 1.2')
  netcdf.putAtt(ncout,globatt,'file_name',namefile)
  netcdf.putAtt(ncout,globatt,'OBC', obc_name)


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
     
  idtime   = netcdf.defVar(ncout,'time_counter','double', dimidt );
    netcdf.putAtt(ncout,idtime,'units',time_units);
    netcdf.putAtt(ncout,idtime,'calendar','gregorian');
    netcdf.putAtt(ncout,idtime,'title','Time');
    netcdf.putAtt(ncout,idtime,'long_name','Time axis');

switch varNC
case {'sossheig'}
  idtemp = netcdf.defVar(ncout,'sossheig'    ,'float' ,[dimidx,dimidy,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','m');
    netcdf.putAtt(ncout,idtemp,'valid_min',-100);
    netcdf.putAtt(ncout,idtemp,'valid_max',100);
    netcdf.putAtt(ncout,idtemp,'long_name','Sea surface height');
case {'vozocrtx'} % U
  idtemp = netcdf.defVar(ncout,'vozocrtx'    ,'float' ,[dimidx,dimidy,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','m s-1');
    netcdf.putAtt(ncout,idtemp,'valid_min',-100);
    netcdf.putAtt(ncout,idtemp,'valid_max',100);
    netcdf.putAtt(ncout,idtemp,'long_name','Zonal velocity');
case {'vomecrty'} % V
  idtemp = netcdf.defVar(ncout,'vomecrty'    ,'float' ,[dimidx,dimidy,dimidt]);
    netcdf.putAtt(ncout,idtemp,'units','m s-1');
    netcdf.putAtt(ncout,idtemp,'valid_min',-100);
    netcdf.putAtt(ncout,idtemp,'valid_max',100);
    netcdf.putAtt(ncout,idtemp,'long_name','Meridional velocity');
end
    
netcdf.endDef(ncout);

% Add y dimension to data
tmpdatain(:,1,:)=datain; datain=tmpdatain;
tmpAlon(:,1)=Alon; Alon=tmpAlon;
tmpAlat(:,1)=Alat; Alat=tmpAlat;
clear tmpAlon tmpAlat tmpdatain


% PUT VARIABLES
netcdf.putVar(ncout,idnavlon,Alon)
netcdf.putVar(ncout,idnavlat,Alat)
netcdf.putVar(ncout,idtime  ,0,NT,time_counter)
netcdf.putVar(ncout,idtemp,datain)

netcdf.close(ncout)

end