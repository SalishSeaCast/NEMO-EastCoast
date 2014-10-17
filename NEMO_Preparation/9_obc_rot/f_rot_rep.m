function [prot] = f_rot_rep (pxin,pyin,cd_type,cdtodo,coordinates,xi,yi)
% SUBROUTINE rot_rep ( pxin, pyin, cd_type, cdtodo, prot )
% !!----------------------------------------------------------------------
% !!                  ***  ROUTINE rot_rep  ***
% !!
% !! ** Purpose :   Rotate the Repere: Change vector componantes between
% !!                geographic grid <--> stretched coordinates grid.
% !!
% !! History :
% !!   9.2  !  07-04  (S. Masson)  
% !!                  (O. Marti ) Original code (repere and repcmo)
% !!----------------------------------------------------------------------
% !! * Arguments 
% REAL(wp), DIMENSION(jpi,jpj), INTENT( IN ) ::   pxin, pyin   ! vector componantes
% CHARACTER(len=1),             INTENT( IN ) ::   cd_type      ! define the nature of pt2d array grid-points
% CHARACTER(len=5),             INTENT( IN ) ::   cdtodo       ! specify the work to do:
% !!                                                           ! 'en->i' east-north componantes to model i componante
% !!                                                           ! 'en->j' east-north componantes to model j componante
% !!                                                           ! 'ij->e' model i-j componantes to east componante
% !!                                                           ! 'ij->n' model i-j componantes to east componante
% REAL(wp), DIMENSION(jpi,jpj), INTENT(out) ::   prot      
%
% !!----------------------------------------------------------------------
% ! Initialization of gsin* and gcos* at first call
% ! -----------------------------------------------
% CALL angle       ! initialization of the transformation
persistent gsint gcost gsinu gcosu gsinv gcosv gsinf gcosf
if isempty(gsint)
   display('COMPUTE ANGLES')
   [gsint,gcost,gsinu,gcosu,gsinv,gcosv,gsinf,gcosf]=f_opa_angle(coordinates,xi,yi);
else
   ssxi=size(xi);
   ssyi=size(yi);
   ssgsint=size(gsint);
   if (ssxi(2)~=ssgsint(2)) || (ssyi(2)~=ssgsint(1))
     display('RECOMPUTE ANGLES')
     clear gsint gcost gsinu gcosu gsinv gcosv gsinf gcosf
     [gsint,gcost,gsinu,gcosu,gsinv,gcosv,gsinf,gcosf]=f_opa_angle(coordinates,xi,yi);
   end
end


switch cdtodo
  case 'en->i'
    switch cd_type
      case 'T' 
        prot(:,:) = pxin(:,:) .* gcost(:,:) + pyin(:,:) .* gsint(:,:);
      case 'U'
        prot(:,:) = pxin(:,:) .* gcosu(:,:) + pyin(:,:) .* gsinu(:,:);
      case 'V'
        prot(:,:) = pxin(:,:) .* gcosv(:,:) + pyin(:,:) .* gsinv(:,:);
      case 'F'
        prot(:,:) = pxin(:,:) .* gcosf(:,:) + pyin(:,:) .* gsinf(:,:);
    end
    
  case 'en->j'
    switch cd_type
      case 'T'
        prot(:,:) = pyin(:,:) .* gcost(:,:) - pxin(:,:) .* gsint(:,:);
      case 'U'
        prot(:,:) = pyin(:,:) .* gcosu(:,:) - pxin(:,:) .* gsinu(:,:);
      case 'V'
        prot(:,:) = pyin(:,:) .* gcosv(:,:) - pxin(:,:) .* gsinv(:,:);  
      case 'F'
        prot(:,:) = pyin(:,:) .* gcosf(:,:) - pxin(:,:) .* gsinf(:,:);
    end
    
  case 'ij->e'
    switch cd_type
      case 'T'
        prot(:,:) = pxin(:,:) .* gcost(:,:) - pyin(:,:) .* gsint(:,:);
      case 'U'
        prot(:,:) = pxin(:,:) .* gcosu(:,:) - pyin(:,:) .* gsinu(:,:);
      case 'V'
        prot(:,:) = pxin(:,:) .* gcosv(:,:) - pyin(:,:) .* gsinv(:,:);
      case 'F'
        prot(:,:) = pxin(:,:) .* gcosf(:,:) - pyin(:,:) .* gsinf(:,:);
    end
    
  case 'ij->n'
    switch cd_type
      case 'T'
        prot(:,:) = pyin(:,:) .* gcost(:,:) + pxin(:,:) .* gsint(:,:);
      case 'U'
        prot(:,:) = pyin(:,:) .* gcosu(:,:) + pxin(:,:) .* gsinu(:,:);
      case 'V'
        prot(:,:) = pyin(:,:) .* gcosv(:,:) + pxin(:,:) .* gsinv(:,:);
      case 'F'
        prot(:,:) = pyin(:,:) .* gcosf(:,:) + pxin(:,:) .* gsinf(:,:);
    end
end
end
%!    END SUBROUTINE rot_rep

