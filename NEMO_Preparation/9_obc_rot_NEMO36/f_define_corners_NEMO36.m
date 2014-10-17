function [corners] = f_define_corners_NEMO36(nbdyind,nbdybeg,nbdyend,obcname,nbptsbound)
% NAME: f_define_corners_NEMO36
%
% AUTHOR: Fatemeh Chegini
% Adapted from JP's f_define_corners script
%
% DATE: Oct 2014
%
% REVISIONS: 
%
% DESCRIPTION: Define the four corners of the selected obc 
%
%
%--------------------------------------------------------------------------

switch obcname
  case 'east'  %       Y                 X   
    corners(1,:)=[ nbdybeg,             nbdyind-nbptsbound+1];  % lower left
    corners(2,:)=[ nbdyend,             nbdyind-nbptsbound+1];  % lower right 
    corners(3,:)=[ nbdybeg,             nbdyind              ];  % upper left
    corners(4,:)=[ nbdyend,             nbdyind              ];  % upper right   
    
  case 'west' %       Y                  X
    corners(1,:)=[ nbdybeg,             nbdyind                 ];  % lower left
    corners(2,:)=[ nbdyend,             nbdyind                 ];  % lower right 
    corners(3,:)=[ nbdybeg,             nbdyind+nbptsbound-1    ];  % upper left
    corners(4,:)=[ nbdyend,             nbdyind+nbptsbound-1    ];  % upper right
 
  case 'north' %       Y                 X
    corners(1,:)=[ nbdyind-nbptsbound+1,  nbdybeg             ];  % lower left
    corners(2,:)=[  nbdyind              ,  nbdybeg           ];  % lower right
    corners(3,:)=[ nbdyind-nbptsbound+1,  nbdyend             ];  % upper left
    corners(4,:)=[  nbdyind              ,  nbdyend           ];  % upper right
     
  case 'south' %       Y                 X    
    corners(1,:)=[    nbdyind,                 nbdybeg        ];  % lower left
    corners(2,:)=[    nbdyind+nbptsbound-1,    nbdybeg        ];  % lower right
    corners(3,:)=[    nbdyind,                 nbdyend        ];  % upper left
    corners(4,:)=[    nbdyind+nbptsbound-1,    nbdyend        ];  % upper right             
  otherwise
    display( 'ERROR !')
end
display([ num2str(corners) ])
end

