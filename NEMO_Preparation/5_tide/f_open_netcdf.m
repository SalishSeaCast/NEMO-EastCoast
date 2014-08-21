%
% NAME: f_open_netcdf.m
%
% AUTHOR: Jean-Philippe PAQUIN
%
% DATE: Novembre 2012
%
% REVISIONS: 
%
%
% DESCRIPTION: Get netCDF Variable Name from file name prefix
%
% INPUTS:  File Prefix and missing value
%
% OUTPUTS: netCDF Varname
%
%-----------------------------------------------------------------------
function [data]=f_open_netcdf(path,file,varNC)

  filein=([path '/' file]);
  
  ncid=netcdf.open(filein,'NC_NOWRITE');
  varid = netcdf.inqVarID(ncid,varNC);
  if isempty(ncid), error('## Bad netcdf operation.'), end
  [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,varid);
  factor=1.; offset=0;
  for att=1:natts-1
     attname = netcdf.inqAttName(ncid,varid,att);
     if     strcmp(attname,'add_offset')
        offset = netcdf.getAtt(ncid,varid,attname);
     elseif strcmp(attname,'scale_factor')
        factor = netcdf.getAtt(ncid,varid,attname);
     end
  end
%  tmp_data=netcdf.getvar(ncid,varid,'double');
  tmp_data=netcdf.getVar(ncid,varid,'double');
  data=tmp_data*factor+offset;
  netcdf.close(ncid)
  clear filein path file ncid varid tmp_data 
  clear varname xtype dimids natts factor offset

end
