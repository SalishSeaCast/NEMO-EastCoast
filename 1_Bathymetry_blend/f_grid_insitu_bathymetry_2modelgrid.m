function [GriddedInsituDepth] = f_grid_insitu_bathymetry_2modelgrid ...
    (InsituData,DestinationGrid)
% =====================================================================
% 
% NAME: Grid_Hermanbathy_2modelgrid.m
% 
% Description:
% Script for gridding Processed Herman Bathymetry data to NEMO grids
% using median value 
%
%======================================================================
% Input:
% 1. InsituData: Insitu Processed data file
%  InsituData.lon : longitude (range: -180:180)
%  InsituData.lat : latitude  (range: -90:90)
%  Insitude.depth : depth
%
% 2. Destination grid data
%  DestinationGrid.lon : longitude
%  DestinationGrid.lat : latitude
% =====================================================================
% Output:
% Gridded insitu bathymetry to destination grids
% 
% =====================================================================
% Process Steps (Adam's idea): 
% 1. Find the nearest NEMO grid for each Herman data point,
% 2. If a grid has points around it, take the median number of them;
% 3. If not, put 999999 as flag of no data.
% 4. Multiply with the NEMO mask to get final result.
% 
% =====================================================================
% Author: Ji Lei (Ji.Lei@dfo-mpo.gc.ca) 
% Modified: Fatemeh Chegini (fatemeh.chegini@dal.ca) Feb 2014


%% Start gridding process =============================================

%Find nearest destinate NEMO cell for each data point 

sizeDestination=size(DestinationGrid.lon);
ndata=length(DestinationGrid.lon(:)); 
surf=griddata(DestinationGrid.lon(:),DestinationGrid.lat(:), ...
[1:ndata],InsituData.lon(:),InsituData.lat(:),'nearest');

%Searching the points that shares the same nearest box

GriddedInsituDepth=DestinationGrid.lon(:);
for ind=1:ndata
    dum=find(surf==ind);
    if (length(dum) >=1)
        GriddedInsituDepth(ind)=median(InsituData.depth(dum));
    else
        GriddedInsituDepth(ind)=NaN;
    end
end

GriddedInsituDepth=reshape(GriddedInsituDepth,sizeDestination);

GriddedInsituDepth(:,1)=NaN;
GriddedInsituDepth(1,:)=NaN;






    
