function [outV] = f_interp_scalar_2D(varin, lon, lat, tgtlon, tgtlat)
                                  
% NAME: f_interp_scalar_2D
%
% AUTHOR: J.-P. Paquin 
%
% DATE: Feb14
%
% REVISIONS: 
%              FCH : modified for 2D variables only such as bathymetry
%
% DESCRIPTION: interpolate scalar field from SOURCE data to a DESTINATION grid 
%              *** Hypothesis: DESTINATION grid might not be a subset
%              *** of the SOURCE grid so the following steps are necessary:
%              1- Extrapolate ocean data on land (floodnan) 
%              2- Interpolate from SOURCE to DESTINATION
% 
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

[NY,NX]=size(varin);
[ny,nx]=size(tgtlon);

% --- Land extrapolation of SOURCE data (floodnan3_opa)
mod_mask=zeros(NY,NX);
tmpvar=zeros(ny,nx);
 
datain = varin ;
%figure ; pcolor(datain) ; shading flat ; colorbar
if isnan(datain)
    tmpvar(:,:)=NaN;
else
    tmp=floodnan3_opa(datain,mod_mask,1);
    %      tmp=datain;
    size(tmp)
    tmp2=griddata(OXU,OYU,tmp,AXU,AYU);
    tmp3=1e20*ones(size(tmp)); tmp3=tmp2;
    tmpvar(:,:) = tmp3;
    clear tmp tmp2 tmp3
end
outV=tmpvar;

end % end of function