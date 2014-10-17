% Program for writing OBC files
% JP Paquin - Jun2013 : Adapt for netcdf 3.6.2
%                      Older MATLAB UNTESTED
%
% Fatemeh Chegini - Oct2014 : Adapt for NEMO 3.6 for 3D fields
%
% ---- Copied from NEMO book version 3.4:
% Each open boundary set is defined as a list of points. The information is stored
% in the arrays nbi, nbj, and nbr in the idx bdy structure. The nbi and nbj arrays
% define the local (i, j) indices of each point in the boundary zone and the nbr array
% defines the discrete distance from the boundary with nbr = 1 meaning that the
% point is next to the edge of the model domain and nbr > 1 showing that the point
% is increasingly further away from the edge of the model domain
%
% The data files contain the data arrays in the order in which the points are defined
% in the nbi and nbj arrays. The data arrays are dimensioned on : a time dimension ;
% xb which is the index of the boundary data point in the horizontal ; and yb which
% is a degenerate dimension of 1 to enable the file to be read by the standard NEMO
% I/O routines. The 3D fields also have a depth dimension.
% 
% 1. The data points must be in order of increasing nbr, ie. all the nbr = 1 points,
% then all the nbr = 2 points etc.
% 2. All the data for a particular boundary set must be in the same order. (Prior
% to 3.4 it was possible to define barotropic data in a different order to the data
% for tracers and baroclinic velocities).


function f_writenetcdf_3D(namefile,varNC,...
                             Alon,Alat,time_counter,time_units,datain,...
                             obc_name)
datain=squeeze(datain);
Alon=squeeze(Alon);
Alat=squeeze(Alat);

[NZ,NY,NX,NT]=size(datain);
nxb=NY*NX;

switch obc_name
    case 'north'
        Alon=Alon(NY:-1:1,:);
        Alat=Alat(NY:-1:1,:);
        datain=datain(:,NY:-1:1,:,:);
    case 'east'
        Alon=Alon(:,NX:-1:1);
        Alat=Alat(:,NX:-1:1);
        datain=datain(:,:,NX:-1:1,:);
end
        
dataline=reshape(datain,NZ,1,nxb,NT);
lonline=reshape(Alon,1,nxb);
latline=reshape(Alat,1,nxb);
    
ncout=netcdf.create(namefile,'NC_write');
% DEFINE DIMENSIONS AND ATTRIBUTES
  dimidx = netcdf.defDim(ncout,'x',nxb);
  dimidy = netcdf.defDim(ncout,'y',1); 
  dimidz = netcdf.defDim(ncout,'z',NZ);
  dimidt = netcdf.defDim(ncout,'time_counter',netcdf.getConstant('NC_UNLIMITED'));

  % Global attributes
  globatt=netcdf.getConstant('NC_GLOBAL');
  netcdf.putAtt(ncout,globatt,'Description','OBC Source from ANNA')
  netcdf.putAtt(ncout,globatt,'Author','Fatemeh Chegini')
  netcdf.putAtt(ncout,globatt,'Date',date)
  netcdf.putAtt(ncout,globatt,'NEMO_version','3.6')
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
    case {'votemper'}
        idvar = netcdf.defVar(ncout,'votemper'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
        netcdf.putAtt(ncout,idvar,'units','deg_C');
        netcdf.putAtt(ncout,idvar,'valid_min',-2.00);
        netcdf.putAtt(ncout,idvar,'valid_max',40.0);
        netcdf.putAtt(ncout,idvar,'long_name','votemper');
    case {'vosaline'}
        idvar = netcdf.defVar(ncout,'vosaline'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
        netcdf.putAtt(ncout,idvar,'units','PSU');
        netcdf.putAtt(ncout,idvar,'valid_min',0. );
        netcdf.putAtt(ncout,idvar,'valid_max',50.);
        netcdf.putAtt(ncout,idvar,'long_name','vosaline');
        
    case {'vozocrtx'} % U
        idvar = netcdf.defVar(ncout,'vozocrtx'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
        netcdf.putAtt(ncout,idvar,'units','m s-1');
        netcdf.putAtt(ncout,idvar,'valid_min',-100);
        netcdf.putAtt(ncout,idvar,'valid_max',100);
        netcdf.putAtt(ncout,idvar,'long_name','Zonal velocity');
    case {'vomecrty'} % V
        idvar = netcdf.defVar(ncout,'vomecrty'    ,'float' ,[dimidx,dimidy,dimidz,dimidt]);
        netcdf.putAtt(ncout,idvar,'units','m s-1');
        netcdf.putAtt(ncout,idvar,'valid_min',-100);
        netcdf.putAtt(ncout,idvar,'valid_max',100);
        netcdf.putAtt(ncout,idvar,'long_name','Meridional velocity');
end

netcdf.endDef(ncout);


% PUT VARIABLES
netcdf.putVar(ncout,idnavlon,lonline)
netcdf.putVar(ncout,idnavlat,latline)
netcdf.putVar(ncout,idtime  ,0,NT,time_counter)
netcdf.putVar(ncout,idvar,dataline)

netcdf.close(ncout)

end