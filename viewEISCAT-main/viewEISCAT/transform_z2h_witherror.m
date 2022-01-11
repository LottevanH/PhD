function [ vh_E, vh_N, vh_h ] = transform_z2h_witherror( vz_N, vz_E, vz_z, ...
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
% R=[  0  sin(inc) cos(inc);                                %  (R) matrix
%      1     0        0    ;
%      0  cos(inc) -sin(inc)];
R = [0          1       0;
    sin(inc)    0       cos(inc);
    cos(inc)    0       sin(inc)];
err_R=[];

[m, n] = size(vz_N.val);
vv_z = [vz_N.val(:)'; vz_E.val(:)'; vz_z.val(:)'];
err_vv_z = [vz_N.err(:)'; vz_E.err(:)'; vz_z.err(:)'];


vv_h = R * vv_z;
err_vv_h = sqrt(R.^2 * err_vv_z.^2);

vh_E.val = reshape(vv_h(1, :), [m, n]);
vh_N.val = reshape(vv_h(2, :), [m, n]);
vh_h.val = reshape(vv_h(3, :), [m, n]);

vh_E.err = reshape(err_vv_h(1, :), [m, n]);
vh_N.err = reshape(err_vv_h(2, :), [m, n]);
vh_h.err = reshape(err_vv_h(3, :), [m, n]);

end
