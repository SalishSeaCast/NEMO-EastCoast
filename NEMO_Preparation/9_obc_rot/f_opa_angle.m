%   SUBROUTINE angle
function [gsint,gcost,gsinu,gcosu,gsinv,gcosv,gsinf,gcosf]=f_opa_angle(coordinates,xi,yi)
% !!----------------------------------------------------------------------
% !!                  ***  ROUTINE opa_angle  ***
% !! 
% !! ** Purpose :   Compute angles between model grid lines and the North direction
% !!
% !! ** Method  :
% !!
% !! ** Action  :   Compute (gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf) arrays:
% !!      sinus and cosinus of the angle between the north-south axe and the 
% !!      j-direction at t, u, v and f-points
% !!
% !! History :
% !!   7.0  !  96-07  (O. Marti )  Original code
% !!   8.0  !  98-06  (G. Madec )
% !!   8.5  !  98-06  (G. Madec )  Free form, F90 + opt.
% !!   9.2  !  07-04  (S. Masson)  Add T, F points and bugfix in cos lateral boundary
% !!----------------------------------------------------------------------
% !! * local declarations
%    INTEGER ::   ji, jj      ! dummy loop indices
%    REAL(wp) ::   &
%    zlam, zphi,            &  ! temporary scalars
%    zlan, zphh,            &  !    "         "
%    zxnpt, zynpt, znnpt,   &  ! x,y components and norm of the vector: T point to North Pole
%    zxnpu, zynpu, znnpu,   &  ! x,y components and norm of the vector: U point to North Pole
%    zxnpv, zynpv, znnpv,   &  ! x,y components and norm of the vector: V point to North Pole
%    zxnpf, zynpf, znnpf,   &  ! x,y components and norm of the vector: F point to North Pole
%    zxvvt, zyvvt, znvvt,   &  ! x,y components and norm of the vector: between V points below and above a T point
%    zxffu, zyffu, znffu,   &  ! x,y components and norm of the vector: between F points below and above a U point
%    zxffv, zyffv, znffv,   &  ! x,y components and norm of the vector: between F points left  and right a V point
%    zxuuf, zyuuf, znuuf       ! x,y components and norm of the vector: between U points below and above a F point
% !!----------------------------------------------------------------------
% 
% ! ============================= !
% ! Compute the cosinus and sinus !
% ! ============================= !
% ! (computation done on the north stereographic polar plane)
% 
%    DO jj = 2, jpjm1
%      DO ji = fs_2, jpi   ! vector opt.

%- Read Structure of coordinates
glamu=coordinates.glamu;
gphiu=coordinates.gphiu;
glamv=coordinates.glamv;
gphiv=coordinates.gphiv;
glamt=coordinates.glamt;
gphit=coordinates.gphit;
glamf=coordinates.glamf;
gphif=coordinates.gphif;

gcost=zeros(size(glamt));
gcosu=zeros(size(glamt));
gcosv=zeros(size(glamt));
gcosf=zeros(size(glamt));
gsint=zeros(size(glamt));
gsinu=zeros(size(glamt));
gsinv=zeros(size(glamt));
gsinf=zeros(size(glamt));

[jpj jpi] = size(glamt);

rpi=pi;
rad=pi/180;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   JPPAQUIN : OPTIMIZED CODE FOR ANGLE CALCULATION ...
%             ! north pole direction & modulous (at t-point)
zlam = rad .* glamt(2:jpj-1,2:jpi);
zphi = rpi/4. - (rad/2.) .*gphit(2:jpj-1,2:jpi);
zxnpt = - 2. .* cos(zlam) .* tan(zphi);
zynpt = - 2. .* sin(zlam) .* tan(zphi);
znnpt = zxnpt.*zxnpt + zynpt.*zynpt;
clear zlam zphi            

%             ! north pole direction & modulous (at u-point)
zlam = rad .* glamu(2:jpj-1,2:jpi);
zphi = rpi/4. - (rad/2.) .*gphiu(2:jpj-1,2:jpi);
zxnpu = - 2. .* cos(zlam) .* tan(zphi);
zynpu = - 2. .* sin(zlam) .* tan(zphi);
znnpu = zxnpu.*zxnpu + zynpu.*zynpu;
clear zlam zphi

%             ! north pole direction & modulous (at v-point)
zlam = rad .* glamv(2:jpj-1,2:jpi);
zphi = rpi/4. - (rad/2.) .* gphiv(2:jpj-1,2:jpi);
zxnpv = - 2. .* cos(zlam) .* tan(zphi);
zynpv = - 2. .* sin(zlam) .* tan(zphi);
znnpv = zxnpv.*zxnpv + zynpv.*zynpv;
clear zlam zphi

%             ! north pole direction & modulous (at f-point)
zlam = rad .* glamf(2:jpj-1,2:jpi);
zphi = rpi/4. - (rad/2.) .* gphif(2:jpj-1,2:jpi);
zxnpf = - 2. .* cos(zlam) .* tan(zphi);
zynpf = - 2. .* sin(zlam) .* tan(zphi);
znnpf = zxnpf.*zxnpf + zynpf.*zynpf;
clear zlam zphi

%             ! j-direction: v-point segment direction (around t-point)
zlam = rad .* glamv(2:jpj-1,2:jpi);
zphi = rpi/4. - (rad/2.) .* gphiv(2:jpj-1,2:jpi);
zlan = rad .* glamv(1:jpj-2,2:jpi);
zphh = rpi/4. - (rad/2.) .* gphiv(1:jpj-2,2:jpi);
zxvvt =  2. .* cos(zlam) .* tan(zphi)   ...
      -  2. .* cos(zlan) .* tan(zphh);
zyvvt =  2. .* sin(zlam) .* tan(zphi)   ...
      -  2. .* sin(zlan) .* tan(zphh);
znvvt = sqrt( znnpt .* ( zxvvt.*zxvvt + zyvvt.*zyvvt )  );
znvvt = max( znvvt, 1.e-14 );
clear zlam zphi zlan zphh

%             ! j-direction: f-point segment direction (around u-point)
zlam = rad .* glamf(2:jpj-1,2:jpi);
zphi = rpi/4. - (rad/2.) .* gphif(2:jpj-1,2:jpi);
zlan = rad .* glamf(1:jpj-2,2:jpi);
zphh = rpi/4. - (rad/2.) .* gphif(1:jpj-2,2:jpi);
zxffu =  2. .* cos(zlam) .* tan(zphi)   ...
      -  2. .* cos(zlan) .* tan(zphh);
zyffu =  2. .* sin(zlam) .* tan(zphi)   ...
      -  2. .* sin(zlan) .* tan(zphh);
znffu = sqrt( znnpu .* ( zxffu.*zxffu + zyffu.*zyffu )  );
znffu = max( znffu, 1.e-14 );
clear zlam zphi zlan zphh

%             ! i-direction: f-point segment direction (around v-point)
zlam = rad.*glamf(2:jpj-1,2:jpi);
zphi = rpi/4. - (rad/2.) .*gphif(2:jpj-1,2:jpi);
zlan = rad.*glamf(2:jpj-1,1:jpi-1);
zphh = rpi/4. - (rad/2.) .*gphif(2:jpj-1,1:jpi-1);
zxffv =  2. .* cos(zlam) .* tan(zphi)   ...
      -  2. .* cos(zlan) .* tan(zphh);
zyffv =  2. .* sin(zlam) .* tan(zphi)   ...
      -  2. .* sin(zlan) .* tan(zphh);
znffv = sqrt( znnpv .* ( zxffv.*zxffv + zyffv.*zyffv )  );
znffv = max( znffv, 1.e-14 );
clear zlam zphi zlan zphh

%             ! j-direction: u-point segment direction (around f-point)
zlam = rad.*glamu(3:jpj  ,2:jpi);
zphi = pi/4. - (rad/2.) .* gphiu(3:jpj  ,2:jpi);
zlan = rad.*glamu(2:jpj-1,2:jpi);
zphh = pi/4. - (rad/2.) .* gphiu(2:jpj-1,2:jpi);
zxuuf =  2. .* cos(zlam) .* tan(zphi )   ...
      -  2. .* cos(zlan) .* tan(zphh);
zyuuf =  2. .* sin(zlam) .* tan(zphi)   ...
      -  2. .* sin(zlan) .* tan(zphh);
znuuf = sqrt( znnpf .* ( zxuuf.*zxuuf + zyuuf.*zyuuf )  );
znuuf = max( znuuf, 1.e-14 );
clear zlam zphi zlan zphh            

%             ! cosinus and sinus using scalar and vectorial products
gsint(2:jpj-1,2:jpi) = ( zxnpt.*zyvvt - zynpt.*zxvvt ) ./ znvvt;
gcost(2:jpj-1,2:jpi) = ( zxnpt.*zxvvt + zynpt.*zyvvt ) ./ znvvt;

gsinu(2:jpj-1,2:jpi) = ( zxnpu.*zyffu - zynpu.*zxffu ) ./ znffu;
gcosu(2:jpj-1,2:jpi) = ( zxnpu.*zxffu + zynpu.*zyffu ) ./ znffu;

gsinf(2:jpj-1,2:jpi) = ( zxnpf.*zyuuf - zynpf.*zxuuf ) ./ znuuf;
gcosf(2:jpj-1,2:jpi) = ( zxnpf.*zxuuf + zynpf.*zyuuf ) ./ znuuf;

%             ! (caution, rotation of 90 degres)
gsinv(2:jpj-1,2:jpi) = ( zxnpv.*zxffv + zynpv.*zyffv ) ./ znffv;
gcosv(2:jpj-1,2:jpi) =-( zxnpv.*zyffv - zynpv.*zxffv ) ./ znffv;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% hopefully the following lines do the same as lbc_lnk in opa does. 
% if not, the error is pretty insignificant anyway

% 1:  first row is copied from second row
      gsint(:,1)=gsint(:,2);
      gcost(:,1)=gcost(:,2);
      gsinu(:,1)=gsinu(:,2);
      gcosu(:,1)=gcosu(:,2);
      gsinv(:,1)=gsinv(:,2);
      gcosv(:,1)=gcosv(:,2);
      gsinf(:,1)=gsinf(:,2);
      gcosf(:,1)=gcosf(:,2);

% 2: last row copied from second last row
      gsint(:,end)=gsint(:,end-1);
      gcost(:,end)=gcost(:,end-1);
      gsinu(:,end)=gsinu(:,end-1);
      gcosu(:,end)=gcosu(:,end-1);
      gsinv(:,end)=gsinv(:,end-1);
      gcosv(:,end)=gcosv(:,end-1);
      gsinf(:,end)=gsinf(:,end-1);
      gcosf(:,end)=gcosf(:,end-1);
     
% first column copied from second column
      gsint(1,:)=gsint(2,:);
      gcost(1,:)=gcost(2,:);
      gsinu(1,:)=gsinu(2,:);
      gcosu(1,:)=gcosu(2,:);
      gsinv(1,:)=gsinv(2,:);
      gcosv(1,:)=gcosv(2,:);
      gsinf(1,:)=gsinf(2,:);
      gcosf(1,:)=gcosf(2,:);
      
% 4: last column copied from second last column
      gsint(end,:)=gsint(end-1,:);
      gcost(end,:)=gcost(end-1,:);
      gsinu(end,:)=gsinu(end-1,:);
      gcosu(end,:)=gcosu(end-1,:);
      gsinv(end,:)=gsinv(end-1,:);
      gcosv(end,:)=gcosv(end-1,:);
      gsinf(end,:)=gsinf(end-1,:);
      gcosf(end,:)=gcosf(end-1,:);

% 5: opa_angle_obc:  modified to cut out a subsection from input (yi,xi)
      gsint=gsint(yi,xi);
      gsinu=gsinu(yi,xi);
      gsinv=gsinv(yi,xi);
      gsinf=gsinf(yi,xi);

      gcost=gcost(yi,xi);
      gcosu=gcosu(yi,xi);
      gcosv=gcosv(yi,xi);
      gcosf=gcosf(yi,xi);

end
