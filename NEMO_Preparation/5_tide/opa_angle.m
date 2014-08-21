%   SUBROUTINE angle
function [gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf] = opa_angle(ncfile)
%       !!----------------------------------------------------------------------
%       !!                  ***  ROUTINE opa_angle  ***
%       !!
%       !! ** Purpose :   Compute angles between model grid lines and the North direction
%       !!
%       !! ** Method  :
%       !!
%       !! ** Action  :   Compute (gsint, gcost, gsinu, gcosu, gsinv, gcosv, gsinf, gcosf) arrays:
%       !!      sinus and cosinus of the angle between the north-south axe and the
%       !!      j-direction at t, u, v and f-points
%       !!
%       !! History :
%       !!   7.0  !  96-07  (O. Marti )  Original code
%       !!   8.0  !  98-06  (G. Madec )
%       !!   8.5  !  98-06  (G. Madec )  Free form, F90 + opt.
%       !!   9.2  !  07-04  (S. Masson)  Add T, F points and bugfix in cos lateral boundary
%       !!----------------------------------------------------------------------
%       !! * local declarations
%       INTEGER ::   ji, jj      ! dummy loop indices
%
%       REAL(wp) ::   &
%          zlam, zphi,            &  ! temporary scalars
%          zlan, zphh,            &  !    "         "
%          zxnpt, zynpt, znnpt,   &  ! x,y components and norm of the vector: T point to North Pole
%          zxnpu, zynpu, znnpu,   &  ! x,y components and norm of the vector: U point to North Pole
%          zxnpv, zynpv, znnpv,   &  ! x,y components and norm of the vector: V point to North Pole
%          zxnpf, zynpf, znnpf,   &  ! x,y components and norm of the vector: F point to North Pole
%          zxvvt, zyvvt, znvvt,   &  ! x,y components and norm of the vector: between V points below and above a T point
%          zxffu, zyffu, znffu,   &  ! x,y components and norm of the vector: between F points below and above a U point
%          zxffv, zyffv, znffv,   &  ! x,y components and norm of the vector: between F points left  and right a V point
%          zxuuf, zyuuf, znuuf       ! x,y components and norm of the vector: between U points below and above a F point
%       !!----------------------------------------------------------------------
%
%       ! ============================= !
%       ! Compute the cosinus and sinus !
%       ! ============================= !
%       ! (computation done on the north stereographic polar plane)
%
%       DO jj = 2, jpjm1
%          DO ji = fs_2, jpi   ! vector opt.
ncload(ncfile);

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

for jj=2:jpj-1
    for ji=2:jpi
        % ! north pole direction & modulous (at t-point)
        zlam = glamt(jj,ji);
        zphi = gphit(jj,ji);
        zxnpt = 0. - 2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        zynpt = 0. - 2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        znnpt = zxnpt*zxnpt + zynpt*zynpt;

        % ! north pole direction & modulous (at u-point)
        zlam = glamu(jj,ji);
        zphi = gphiu(jj,ji);
        zxnpu = 0. - 2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        zynpu = 0. - 2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        znnpu = zxnpu*zxnpu + zynpu*zynpu;

        % ! north pole direction & modulous (at v-point)
        zlam = glamv(jj,ji);
        zphi = gphiv(jj,ji);
        zxnpv = 0. - 2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        zynpv = 0. - 2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        znnpv = zxnpv*zxnpv + zynpv*zynpv;

        % ! north pole direction & modulous (at f-point)
        zlam = glamf(jj,ji);
        zphi = gphif(jj,ji);
        zxnpf = 0. - 2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        zynpf = 0. - 2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. );
        znnpf = zxnpf*zxnpf + zynpf*zynpf;

        % ! j-direction: v-point segment direction (around t-point)
        zlam = glamv(jj,  ji);
        zphi = gphiv(jj,  ji);
        zlan = glamv(jj-1,ji);
        zphh = gphiv(jj-1,ji);
        zxvvt =  2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * cos( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        zyvvt =  2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * sin( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        znvvt = sqrt( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  );
        znvvt = max( znvvt, 1.e-14 );

        % ! j-direction: f-point segment direction (around u-point)
        zlam = glamf(jj  ,ji );
        zphi = gphif(jj  ,ji);
        zlan = glamf(jj-1,ji);
        zphh = gphif(jj-1,ji);
        zxffu =  2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * cos( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        zyffu =  2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * sin( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        znffu = sqrt( znnpu * ( zxffu*zxffu + zyffu*zyffu )  );
        znffu = max( znffu, 1.e-14 );

        % ! i-direction: f-point segment direction (around v-point)
        zlam = glamf(jj,ji  );
        zphi = gphif(jj,ji  );
        zlan = glamf(jj,ji-1);
        zphh = gphif(jj,ji-1);
        zxffv =  2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * cos( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        zyffv =  2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * sin( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        znffv = sqrt( znnpv * ( zxffv*zxffv + zyffv*zyffv )  );
        znffv = max( znffv, 1.e-14 );

        % ! j-direction: u-point segment direction (around f-point)
        zlam = glamu(jj+1,ji);
        zphi = gphiu(jj+1,ji);
        zlan = glamu(jj  ,ji);
        zphh = gphiu(jj  ,ji);
        zxuuf =  2. * cos( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * cos( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        zyuuf =  2. * sin( rad*zlam ) * tan( rpi/4. - rad*zphi/2. )   ...
              -  2. * sin( rad*zlan ) * tan( rpi/4. - rad*zphh/2. );
        znuuf = sqrt( znnpf * ( zxuuf*zxuuf + zyuuf*zyuuf )  );
        znuuf = max( znuuf, 1.e-14 );

        % ! cosinus and sinus using scalar and vectorial products
        gsint(jj,ji) = ( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt;
        gcost(jj,ji) = ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt;

        gsinu(jj,ji) = ( zxnpu*zyffu - zynpu*zxffu ) / znffu;
        gcosu(jj,ji) = ( zxnpu*zxffu + zynpu*zyffu ) / znffu;

        gsinf(jj,ji) = ( zxnpf*zyuuf - zynpf*zxuuf ) / znuuf;
        gcosf(jj,ji) = ( zxnpf*zxuuf + zynpf*zyuuf ) / znuuf;

        % ! (caution, rotation of 90 degres)
        gsinv(jj,ji) = ( zxnpv*zxffv + zynpv*zyffv ) / znffv;
        gcosv(jj,ji) =-( zxnpv*zyffv - zynpv*zxffv ) / znffv;

    end
end
%          END DO
%       END DO

%       ! =============== !
%       ! Geographic mesh !
%       ! =============== !
%

% MD:
% COMMENT THIS OUT - IT SEEMS TO CAUSE PROBLEMS
% for jj=2:jpj-1
%     for ji=2:jpi
%             if( mod( abs( glamv(jj,ji) - glamv(jj-1,ji) ), 360. ) < 1.e-5 )
%                gsint(jj,ji) = 0.;
%                gcost(jj,ji) = 1.;
%             end
%             if( mod( abs( glamf(jj,ji) - glamf(jj-1,ji) ), 360. ) < 1.e-5 )
%                gsinu(jj,ji) = 0.;
%                gcosu(jj,ji) = 1.;
%             end
%             if(      abs( gphif(jj,ji) - gphif(jj,ji-1) )         < 1.e-5 )
%                gsinv(jj,ji) = 0.;
%                gcosv(jj,ji) = 1.;
%             end
%             if( mod( abs( glamu(jj,ji) - glamu(jj+1,ji) ), 360. ) < 1.e-5 )
%                gsinf(jj,ji) = 0.;
%                gcosf(jj,ji) = 1.;
%             end
%     end
% end




%       ! =========================== !
%       ! Lateral boundary conditions !
%       ! =========================== !
%
%       ! lateral boundary cond.: T-, U-, V-, F-pts, sgn
%       CALL lbc_lnk ( gcost, 'T', 1. )   ;   CALL lbc_lnk( gsint, 'T', -1. )
%       CALL lbc_lnk ( gcosu, 'U', 1. )   ;   CALL lbc_lnk( gsinu, 'U', -1. )
%       CALL lbc_lnk ( gcosv, 'V', 1. )   ;   CALL lbc_lnk( gsinv, 'V', -1. )
%       CALL lbc_lnk ( gcosf, 'F', 1. )   ;   CALL lbc_lnk( gsinf, 'F', -1. )


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

end
