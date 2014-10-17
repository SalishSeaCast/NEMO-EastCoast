function [corners] = f_define_corners(NY,NX,obc,nbptsbound)
% NAME: f_define_corners
%
% AUTHOR: J.-P. Paquin 
%
% DATE: Feb14
%
% REVISIONS: 
%
% DESCRIPTION: Define the four corners of the selected obc 
%
% NOTE : EXTRACTION OF THE NEMO ROUTINE obc_par_mer.F90 
%        for the NFS036 experiments (Li Zhai)
%***********************************************************************
%    !! * EAST open boundary
%       (...)
%       jpieob = (/jpiglo-2/), & !: i-localization of the East open boundary
%       jpjedt = (/       2/), & !: j-starting indice of the East open boundary
%       jpjeft = (/jpjglo-1  /)    !: j-ending   indice of the East open boundary
% 
%    !! * WEST open boundary
%       (...)
%       jpiwob = (/    2  /), &  !: i-localization of the West open boundary
%       jpjwdt = (/   2   /), &  !: j-starting indice of the West open boundary
%       jpjwft = (/jpjglo-1/)     !: j-ending   indice of the West open boundary
%
%    !! * NORTH open boundary
%       (...)
%       jpjnob = (/jpjglo-2/), &  !: j-localization of the North open boundary
%       jpindt = (/       2/), &  !: i-starting indice of the North open boundary
%       jpinft = (/jpiglo-1/)     !: i-ending   indice of the North open boundary
% 
%    !! * SOUTH open boundary
%       (...)
%       jpjsob = (/       2/), &  !: j-localization of the South open boundary
%       jpisdt = (/       2/), &  !: i-starting indice of the South open boundary
%       jpisft = (/jpiglo-1  /)     !: i-ending   indice of the South open boundary
%***********************************************************************
% 
% CALLED PGM & SCRIPTS: 
%--------------------------------------------------------------------------


%%%%% IMPORTANT NOTE: THE INVERSION REQUIRED FOR THE      %%%%% 
%%%%% NORTHERN AND EASTERN BOUNDARIES IS NOT HANDLED HERE %%%%%
%%%%% BUT JUST BEFORE WRITING THE NETCDF FILES            %%%%%
switch obc
  case 'east'  %       Y                 X   
    corners(1,:)=[    2,                (NX-1)-nbptsbound+1];  % lower left
    corners(2,:)=[ NY-1,                (NX-1)-nbptsbound+1];  % lower right 
    corners(3,:)=[    2,                 NX-1              ];  % upper left
    corners(4,:)=[ NY-1,                 NX-1              ];  % upper right   
    
  case 'west' %       Y                  X
    corners(1,:)=[    2,                 2                 ];  % lower left
    corners(2,:)=[ NY-1,                 2                 ];  % lower right 
    corners(3,:)=[    2,                 2+nbptsbound-1    ];  % upper left
    corners(4,:)=[ NY-1,                 2+nbptsbound-1    ];  % upper right
 
  case 'north' %       Y                 X
    corners(1,:)=[ (NY-1)-nbptsbound+1,  2                 ];  % lower left
    corners(2,:)=[  NY-1              ,  2                 ];  % lower right
    corners(3,:)=[ (NY-1)-nbptsbound+1,  NX-1              ];  % upper left
    corners(4,:)=[  NY-1              ,  NX-1              ];  % upper right
     
  case 'south' %       Y                 X    
    corners(1,:)=[    2,                 2                 ];  % lower left
    corners(2,:)=[    2+nbptsbound-1,    2                 ];  % lower right
    corners(3,:)=[    2,                 NX-1              ];  % upper left
    corners(4,:)=[    2+nbptsbound-1,    NX-1              ];  % upper right             
  otherwise
    display( 'ERROR !')
end
display([ num2str(corners) ])
end

