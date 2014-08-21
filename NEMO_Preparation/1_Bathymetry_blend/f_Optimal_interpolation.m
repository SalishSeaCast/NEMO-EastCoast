
function [Blend] = f_Optimal_interpolation ...
    (lon,lat,Background,Insitu,delx,dely,Ls,sigobs)

% First check if background data has missing points
% and replace with in-situ data

I=find(isnan(Background)==1 & isnan(Insitu)==0);
Backgroundfill=Background; Backgroundfill(I)=Insitu(I);

increment=Insitu-Backgroundfill;

%% Map increment to grid using OI 
%  with Gaspari-Cohn (compact support) corrlation function

increment_pred=nan*increment; 
for lat1=min(lat(:)):dely:max(lat(:))
    for lon1=min(lon(:)):delx:max(lon(:))
        
        lat2=lat1+dely; lon2=lon1+delx;
        
        Ig=find(lat>=lat1&lat<lat2&lon>=lon1&lon<lon2);
        x1=lon(Ig); x2=lat(Ig);
        
        I=find(lat>=lat1&lat<lat2&lon>=lon1&lon<lon2&isnan(increment)==0);
        x1obs=lon(I); x2obs=lat(I); yobs=increment(I);
        
        ypred=f_OI_gaspari_cohn(x1,x2,Ls,x1obs,x2obs,yobs,sigobs);
        increment_pred(Ig)=ypred;
        
    end
end

Blend=Backgroundfill+increment_pred;

