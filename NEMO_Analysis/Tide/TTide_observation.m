% =========================================================================
% TTide_observation
%  
% DESCRIPTION:
% Script for Tidal harmonic analysis of observed stations using TTide
% Download TTide from:
% http://www.isdm-gdsi.gc.ca/isdm-gdsi/twl-mne/index-eng.htm
%
% INPUT: sea level time series 
%  format: read from Environment Canada, cvs file, hourly data
%          8 lines header, data: date, sea level
%  format example: 
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  Station_Name,YARMOUTH
%  Station_Number,365
%  Latitude_Decimal_Degrees,43.833333
%  Longitude_Decimal_Degrees,66.116667
%  Datum,CD
%  Time_Zone,AST
%  SLEV=Observed Water Level
%  Obs_date,SLEV(metres)
%  2014/05/01 00:00,4.62,
%  2014/05/01 01:00,4.22,
%  ...
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% OUTPUT: Amplitude and phase for each constituent and errors
% output is written to file
% number of constituents depend on the length of your time series
%
% Note: Environment Canada time series have occasional gap
%       in this script these gaps are detected
%==========================================================================
% Author: Fatemeh Chegini
% Date: 30/06/2014
%==========================================================================

% ============= User defined parameters ==================================
inputfile='/media/Data/NEMO/Data/TideGauge/365-01-JAN-2013_slev.csv';
outputfile='/media/Data/NEMO/FC/NEMO_Validate/Tides/Observation/Yarmouth_consituents.txt';
% ============= End of user defined parameters ===========================


HEADERLINES=8; DELIMITER=',';
tmpdata = importdata(inputfile, DELIMITER, HEADERLINES);

data.date=tmpdata.textdata(9:end);
data.elev=tmpdata.data(:);
clear tmpdata

%First check to see if data has any gaps (default is hourly data)
ndata=length(data.date);
for i=2:ndata
    if (datenum(data.date(i))-datenum(data.date(i-1))-0.0417)>1e-4
        display(['There is a gap in:',num2str(i)])
        
    end
end

% Call TTide to do harmonic analysis

[TIDESTRUC,XOUT]=t_tide(data.elev(1:96),'output',outputfile);
