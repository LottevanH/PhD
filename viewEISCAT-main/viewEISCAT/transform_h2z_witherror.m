function [ vz_N, vz_E, vz_z ] = transform_h2z_witherror( vh_E, vh_N, vh_h, ...
    inc)
%TRANSFORM_H2Z Summary of this function goes here
%   Detailed explanation goes here
  %addpath /eiscat/MODELS/IGRF/geomag_v60

% [yy, mm, dd] = datevec(dn);
% decimalYear = decyear(yy, mm, dd);
% [magFieldVector,horIntensity,declination,inclination,totalIntensity,    ...
%     magFieldSecVariation,secVariationHorizontal,secVariationDeclination, ...
%     secVariationInclination,secVariationTotal] =    ...
%     igrfmagm(height,lat,lon,decimalYear,12); 
inc=inc/180*pi;
R=[  0  sin(inc) cos(inc);                                %  (R) matrix
     1     0        0    ;
     0  cos(inc) -sin(inc)];
% R = [0          1       0;
%     sin(inc)    0       cos(inc);
%     cos(inc)    0       sin(ind)];
err_R=[];

[m, n] = size(vh_N.val);
vv_h = [vh_E.val(:)'; vh_N.val(:)'; vh_h.val(:)'];
err_vv_h = [vh_E.err(:)'; vh_N.err(:)'; vh_h.err(:)'];


vv_z = R * vv_h;
err_vv_z = sqrt(R.^2 * err_vv_h.^2);

vz_N.val = reshape(vv_z(1, :), [m, n]);
vz_E.val = reshape(vv_z(2, :), [m, n]);
vz_z.val = reshape(vv_z(3, :), [m, n]);

vz_N.err = reshape(err_vv_z(1, :), [m, n]);
vz_E.err = reshape(err_vv_z(2, :), [m, n]);
vz_z.err = reshape(err_vv_z(3, :), [m, n]);
end
