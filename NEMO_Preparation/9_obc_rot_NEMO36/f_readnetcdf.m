%JP Paquin - Jun2013 : Adapt for netcdf 3.6.2
%                      Checks for offsets and scaling factor
%                      Older MATLAB UNTESTED
function[data]=f_readnetcdf(dirin,file,varNC,cdfversion)

fillvalue=NaN;
filein=([dirin '/' file]);
if strcmp(cdfversion,'3.6.2')
    
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
     elseif strcmp(attname,'_FillValue')
        fillvalue = netcdf.getAtt(ncid,varid,attname);
     elseif strcmp(attname,'missing_value')
        fillvalue = netcdf.getAtt(ncid,varid,attname);
     end
     
  end
  tmp_data=netcdf.getVar(ncid,varid,'double');
  if ~isnan(fillvalue)
    tmp_data(tmp_data==fillvalue)=NaN;
  end
  tmp_data=tmp_data.*factor+offset;
  data=permute(tmp_data,[3,2,1,4]);
  netcdf.close(ncid)

else % Older netcdf code : N.B. UNTESTED
  display('**** OLDER CODE UNTESTED ****')
  nc=netcdf(filein);
  data = nc{varNC}(1,:,:,:)  ;
  close(nc);
    
end
end