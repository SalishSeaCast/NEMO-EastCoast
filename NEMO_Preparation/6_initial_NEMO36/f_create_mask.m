%JP Paquin - Jun2013 : Computes 2D mask for land extrapolation
%
function[mask,count]=f_create_mask(data,nbpts)

[inz,iny,inx]=size(data);
mask=data ;
mask(~isnan(mask)) = 0; 
count(inz,1)=0;
for zz=1:inz
for ii=1:inx
for jj=1:iny
    if isnan(mask(zz,jj,ii)) % if land point
        
      minii=max(  1,ii-nbpts);
      maxii=min(inx,ii+nbpts);
      minjj=max(  1,jj-nbpts);
      maxjj=min(iny,jj+nbpts);        
        
      for sii=minii:maxii
      for sjj=minjj:maxjj
          if ( ~isnan(data(zz,sjj,sii)) ) %- if ocean point present within 10 grid points
             mask(zz,jj,ii)=0;
             count(zz)=count(zz)+1;
          end
      end
      end
    end
end
end
end
    
end