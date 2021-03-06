function [lat,lon,glon,glat] = convert_hdf_to_lat(filename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[Time,locn,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5_new(filename);

% read the basic parameters into some usefully named arrays
[y,m,d,h,mn,s]=datevec(Time(1,:));
ut_time=h+mn/60.;
alt=par2D(:,:,2);
ne=par2D(:,:,3);
te=par2D(:,:,4);%.*par2D(:,:,5);
ti=par2D(:,:,5);
azm=par1D(:,1);
eln=par1D(:,2);
rng=par2D(:,:,1);

% add the software path for the loc2gg program to convert from elevation, range and azimuth to geopgraphic latitude and longitude 
addpath 'W:\COURSE MTR & DATA StudentsReadOnly\AGF\AGF_304_software\project_scripts'
 for j= 1:length(rng(:,1));
     for jj=1:length(rng(1,:)); 
   dd=loc2gg(locn,[eln(jj), azm(jj), rng(j,jj)]);
   lat(j,jj)=dd(1);
   lon(j,jj)=dd(2);
   end;
 end;
 
 
  % convert to geomagnetic co-ords using IGRF model - input the year, month
  % and day of the data set *******************************
addpath 'W:\COURSE MTR & DATA StudentsReadOnly\AGF\AGF_304_software\project_scripts\m_map\'
 m_coord('IGRF-geomagnetic',datenum(2016,2,16));  % define the coord system and input the correct date here!!!!
 for j= 1:length(rng(:,1));
     for jj=1:length(rng(1,:));
[glon(j,jj),glat(j,jj)]=m_geo2mag(lon(j,jj),lat(j,jj));   % Convert geographic to geomagnetic
%[LON,LAT]=m_mag2geo(gLON,gLAT);   % This is the inverse
     end
 end
end

