%JP Paquin - Jun2013 : Adapt for netcdf 3.6.2
%                      Older MATLAB UNTESTED
function[tmask,umask,vmask,fmask,mdepth,gdepw,dzz,Alon,Alat]=f_readmeshmask(meshfile,cdfversion)

if strcmp(cdfversion,'3.6.2')
  ncid=netcdf.open(meshfile,'NC_nowrite');
    % --- READ ALL MASKS ---
    varid = netcdf.inqVarID(ncid,'tmask');
    tmp_tmask = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'umask');
    tmp_umask = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'vmask');
    tmp_vmask = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'fmask');
    tmp_fmask = netcdf.getVar(ncid,varid); 
    % --- READ OTHER FIELDS ---
    varid = netcdf.inqVarID(ncid,'gdept_0');
    mdepth = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'gdepw_0');
    gdepw = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'e3w_0');
    dzz = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'nav_lat');
    tmp_lat = netcdf.getVar(ncid,varid);
    varid = netcdf.inqVarID(ncid,'nav_lon');
    tmp_lon = netcdf.getVar(ncid,varid);
  netcdf.close(ncid)

  %- Permute for compatibility with older code.
  tmask=permute(tmp_tmask,[3,2,1,4]);
  umask=permute(tmp_umask,[3,2,1,4]);
  vmask=permute(tmp_vmask,[3,2,1,4]);
  fmask=permute(tmp_fmask,[3,2,1,4]);

  Alon=permute(tmp_lon,[2,1]);
  Alat=permute(tmp_lat,[2,1]);    
  
end
end