function [uout vout] = f_interpuv(uin, vin, lonu, latu, lonv, latv, ...
                                  tgtlonu, tgtlatu, tgtlonv, tgtlatv, ...
                                  z, tgtz, ...
                                  coordGLOB, xi, yi, ...
                                  coordREG , axi, ayi)
% NAME: gen_UVini.m
%
% AUTHOR: Michael Dunphy 
%
% DATE: ???
%
% REVISIONS: 
%    Jun2013: JP Paquin: - Change coordinates files/paths to "structure"
%                          variable containing informations (for function
%                          f_opa_angle.m)
%                        - Update comments
%
%
% DESCRIPTION: obcinterpuv: interpolate u,v from SOURCE data to a
%              DESTINATION grid 
%              *** Hypothesis: DESTINATION grid might not be a subset
%              *** of the SOURCE grid so the following steps are necessary:
%              1: Move SOURCE V from V points to U points
%              2: Rotate this (U,V) at U points into east-north
%              3: Interpolate to DESTINATION points
%              4: Rotate from east-north into DESTINATION i,j coordinates
%              5: Vertical interpolation from SOURCE to DESTINATION levels
%
% NOTE : REQUIRES M_MAP PACKAGE
%
% CALLED PGM & SCRIPTS: f_rot_rep.m (f_opa_angle.m)
%--------------------------------------------------------------------------
m_proj=('mercator');

% prepare projected coodinates for interpolation in x-y space
[OXU,OYU]=m_ll2xy(double(lonu),double(latu)); % orca u points in x-y space
[OXV,OYV]=m_ll2xy(double(lonv),double(latv)); % orca v points in x-y space
[AXU,AYU]=m_ll2xy(double(tgtlonu),double(tgtlatu)); % arc016 u points in x-y space
[AXV,AYV]=m_ll2xy(double(tgtlonv),double(tgtlatv)); % arc016 v points in x-y space


[NZ,NY,NX,NT]=size(uin);
[ny,nx]      =size(tgtlonu);
[nz,~]       =size(tgtz);


% move v-at-v points to v-at-u points for SOURCE data
for zi=1:NZ
  fprintf('%s: SOURCE griddata v-at-v -> v-at-u: zi=%d ...\n',mfilename,zi);
  for ti=1:NT
    tmp=squeeze(vin(zi,:,:,ti)); 
    tmp2=griddata(OXV,OYV,tmp,OXU,OYU);
    tmp3=1e20*ones(size(tmp)); tmp3=tmp2;
    vin(zi,:,:,ti) = tmp3;
  end
end
clear tmp tmp2 tmp3

% now rotate (uin,vin) at the SOURCE u points into east-north
for zi=1:NZ
  fprintf('%s: SOURCE rotate i,j -> east,north: zi=%d ...\n',mfilename,zi);
  for ti=1:NT
    tmpu=squeeze(uin(zi,:,:,ti));
    tmpv=squeeze(vin(zi,:,:,ti));
    uin(zi,:,:,ti) = f_rot_rep(tmpu,tmpv, 'U','ij->e',coordGLOB,xi,yi);
    vin(zi,:,:,ti) = f_rot_rep(tmpu,tmpv, 'U','ij->n',coordGLOB,xi,yi);
  end
end
clear tmpu tmpv

% Interpolate east-north velocity from SOURCE u points
% onto DESTINATION u points and v points.  Then rotate this
% east north velocity into DESTINATION i-j coordinates
u50=zeros(NZ,ny,nx,NT);
v50=zeros(NZ,ny,nx,NT);
for zi=1:NZ
    fprintf('%s: horiz interp to DESTINATION and rotate east,north -> i,j: zi=%d ...\n',mfilename,zi);
    for ti=1:NT
    % east velocity to u,v points
    tmp=squeeze(uin(zi,:,:,ti));
    uup=griddata(OXU,OYU,tmp,AXU,AYU);
    uvp=griddata(OXU,OYU,tmp,AXV,AYV);

    % north velocity to u,v points
    tmp=squeeze(vin(zi,:,:,ti));
    vup=griddata(OXU,OYU,tmp,AXU,AYU);
    vvp=griddata(OXU,OYU,tmp,AXV,AYV);

    % rotate into REGIONAL coordinates
    u50(zi,:,:,ti) = f_rot_rep(uup, vup, 'U','en->i',coordREG,axi,ayi);
    v50(zi,:,:,ti) = f_rot_rep(uvp, vvp, 'V','en->j',coordREG,axi,ayi);
    end
end
clear uup uvp vup vvp


% Now interpolate vertically from SOURCE levels to DESTINATION levels
tmp_moddpth=zeros(nz,1);
for zz=1:nz
  if (  tgtz(zz) < z(1)) 
     tmp_moddpth(zz) = z(1);
  else
     tmp_moddpth(zz) = tgtz(zz);
  end    
end
uout=zeros(nz,ny,nx,NT);
vout=zeros(nz,ny,nx,NT);

fprintf('%s: vertical interpolation from SOURCE to DESTINATION levels: ...\n',mfilename);
for ti=1:NT
  for ii=1:nx
  for jj=1:ny
    tmp=(squeeze(u50(:,jj,ii,ti)));
    uout(:,jj,ii,ti) = interp1q(z,tmp,tmp_moddpth);
    tmp=squeeze(v50(:,jj,ii,ti));
    vout(:,jj,ii,ti) = interp1q(z,tmp,tmp_moddpth);

  end
  end
end

end % end of function