function [outV] = f_interp_scalar(varin, lon, lat, tgtlon, tgtlat, ...
                                  z, tgtz)
% NAME: f_interp_scalar
%
% AUTHOR: J.-P. Paquin 
%
% DATE: Feb14
%
% REVISIONS: 
%
% DESCRIPTION: interpolate scalar field from SOURCE data to a DESTINATION grid 
%              *** Hypothesis: DESTINATION grid might not be a subset
%              *** of the SOURCE grid so the following steps are necessary:
%              1- Extrapolate ocean data on land (floodnan) 
%              2- Interpolate from SOURCE to DESTINATION
%
% NOTE : REQUIRES M_MAP PACKAGE
%
% CALLED PGM & SCRIPTS: 
%              floodnan3_opa
%              interp1q
%--------------------------------------------------------------------------
m_proj=('mercator');

% prepare projected coodinates for interpolation in x-y space
[OXU,OYU]=m_ll2xy(double(lon)   ,double(lat));    % SOURCE grid
[AXU,AYU]=m_ll2xy(double(tgtlon),double(tgtlat)); % DESTINATION (target) grid

[NZ,NY,NX,NT]=size(varin);
[ny,nx]      =size(tgtlon);
[nz,~]       =size(tgtz);


% --- Land extrapolation of SOURCE data (floodnan3_opa)
mod_mask=zeros(NY,NX);
tmpvar=zeros(NZ,ny,nx,NT);
fprintf('%s: horiz interp to DESTINATION and flooding: ');
for zi=1:NZ
  %fprintf('%s: horiz interp to DESTINATION and flooding: zi=%d ...\n',mfilename,zi);
  for ti=1:NT    
    datain = squeeze(varin(zi,:,:,ti)) ;
    %figure ; pcolor(datain) ; shading flat ; colorbar
    if isnan(datain)
      tmpvar(zi,:,:,ti)=NaN;
    else
      tmp=floodnan3_opa(datain,mod_mask,1);
      tmp2=griddata(OXU,OYU,tmp,AXU,AYU);
      tmp3=1e20*ones(size(tmp)); tmp3=tmp2;
      tmpvar(zi,:,:,ti) = tmp3;    
      clear tmp tmp2 tmp3
    end
  end
end

fprintf('%s: vertical interpolation ');
% Now interpolate vertically from SOURCE levels to DESTINATION levels
if ( NZ>1 && nz>1)
  tmp_moddpth=zeros(nz,1);
  for zz=1:nz
    if (  tgtz(zz) < z(1)) 
       tmp_moddpth(zz) = z(1);
    else
       tmp_moddpth(zz) = tgtz(zz);
    end    
  end

  outV=zeros(nz,ny,nx,NT);
  fprintf('%s: vertical interpolation from SOURCE to DESTINATION levels: ...\n',mfilename);
  for ti=1:NT
    for ii=1:nx
    for jj=1:ny
      tmp=(squeeze(tmpvar(:,jj,ii,ti)));
      outV(:,jj,ii,ti) = interp1q(z,tmp,tmp_moddpth);
    end
    end
  end
else % 2D fields no vertical interpolation...
  for ti=1:NT  
    outV(1,:,:,ti)= tmpvar(1,:,:,ti);   
  end
end 
end % end of function