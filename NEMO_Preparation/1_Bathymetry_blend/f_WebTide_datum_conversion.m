
function [CorrectedData] = f_WebTide_datum_conversion ...
    (pathWeb,nodeWeb,inputdata,list_compo)

% NAME: f_WebTide_datum_conversion.m
%
% AUTHOR: Fatemeh Chegini
%
% DATE: Feb 2014
%
%
% DESCRIPTION: Convert datum for bathymetry 
%              from LLWLT (low level water, large tides) 
%              to MSL (mean sea level)
%              using tidal amplitude from Webtide
%
%
% INPUTS : 
% 1. File containing destination grid (inputfile)
%    file format: ASCII
%    file information: Longitude, Latitude , LLWLT depth (neagitve down)
%          
% 2. Webtide files 
%    including tidal amplitudes (.s2c) and node information (.nod) files
%    file format: ASCII
%    Download and install from:
%    www.bio.gc.ca/science/research-recherche/ocean/webtide/index-eng.php
%
% OUTPUTS: 
% 1. Converted depths 
%    file format: ASCII
%    file information: Longitude, Latitude, MSL depth (negative down)
%    
%
% DETAILS:
% This script converts datum for input depths with LLWLT reference to MSL. 
% This is done by adding the tidal amplitude from WebTide to the given 
% input depths. Before this, the tidal amplitudes are interpolated 
% from WebTide unstructured grids to desnitaion grids. 
% The program also extrapolates for those nodes that are not in 
% WebTide domain. To avoid errors specially when extrapolating, 
% inter/extrapolation is done using the nearest neighbour method. 
% The tidal constituents to be included is defined by the
% user. The program only sums up the amplitudes for defined constituents.
%
% NOTES: Parts of this script were taken from interp_tides_Webtide 
%        by Jean-Philippe Paquin
%        and also the some of the calculation is from the bathycor.c
%        program by David Greenberg. I have compared the results with
%        bathycor output and they are similar.
% 
%======================================================================


%% 1.1. Read input longitude and latitude

loninput=inputdata(:,1) ;  % longitude
latinput=inputdata(:,2);   % latitude
depthinput=inputdata(:,3); % LLWLT depth

nnodesinput=length(loninput); %number of given points


%% 1.2. Read WebTide longitudes and latitude

latlonweb  =([pathWeb '/' nodeWeb]);           % Webtide Lat-Lon
HEADERLINES=0; DELIMITER=' ';
tmpdata = importdata(latlonweb, DELIMITER, HEADERLINES);
lonWeb=tmpdata(:,2)  ; % longitude
latWeb=tmpdata(:,3)  ; % latitude
ii=find(lonWeb>0)  ;lonWeb(ii) = lonWeb(ii)-360 ; % Remove lines where long
clear tmpdata                                  % are superior to 0

nnodesWeb=length(lonWeb); % number of webtide nodes

WebTotalAmplitude=zeros(nnodesWeb,1); 

%% 1.3. Read tidal Amplitudes for all constituents 
%  and sum them up to get total value

for mycomponents = list_compo
  compo = mycomponents{1}; mycomponents; 
  
  display([ 'Reading data for tidal component : ' compo ])
  tidefile   =([pathWeb '/' compo '.barotropic.s2c']); % elevation
    
  HEADERLINES=3; DELIMITER=' ';      
  tmp= importdata(tidefile, DELIMITER, HEADERLINES);
  tmpdata=tmp.data ; % clear tmp 
  WebCompAmplitude=tmpdata(:,2)  ; % H amplitude
  clear tmpdata
  
  WebTotalAmplitude=WebTotalAmplitude+WebCompAmplitude; %sum up all amplitudes
  
end

%% Interpolate 

display('Interpolating data on given grids')

InterpAmplitude= griddata(lonWeb,latWeb,WebTotalAmplitude,loninput,latinput,'nearest');

% Although we're using nearest, still check for wierd results
for i=1:nnodesinput
    if InterpAmplitude(i)<0 
        display([ ' Wierd negative amplitude in point ', i , '.... setting to zero '])
        InterpAmplitude(i)=0;
    end
end

%% Correct datum from LLWLT to MSL
CorrectedDepth = zeros(nnodesinput,1);
for i=1:nnodesinput
    
    % Calculation of Ffact was copied from bathycor.c program written by
    % David Greenberg
    Amp=InterpAmplitude(i) ;
    if(Amp<=2.)
        Ffact=0.85;
    elseif(Amp>=5)
        Ffact=1.0;
    elseif(Amp<5 && Amp>2)
        Ffact = .85 + .15*(Amp-2.)/3.;
    end
    
    depthpositivedown = -depthinput(i);
    CorrectedDepth(i) = -(depthpositivedown+Ffact*Amp) ; %negative down
  
end

CorrectedData=[loninput,latinput,CorrectedDepth];




