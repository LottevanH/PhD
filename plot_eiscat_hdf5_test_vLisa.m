clear all
% program to use the function in guisdap load_param_hdf5.m to load up data
% from a hdf5 file which has been downloaded from madrigal

% Instructions:
% 1. Download the hdf5 data file from madrigal
% 2. start guisdap in matlab and then run the following command (e.g):
% [Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('D:\EISCAT_data\analyzed\eiscat_website_hdf5\EISCAT_2019-11-20_folke_64@42mb.hdf5')
% Variable names: 
% par2D: [Ran,Alt,Ne,Te,Ti,Vi,Coll,Comp,Res]
% par1D: [Az,El,Pt,Tsys]
% rpar2D: [ppRan,ppAlt,ppNe]
% err2d_id  = {'var_Ne' 'var_Tr' 'var_Ti' 'var_Vi' 'var_Collf'};
%   Note 1: 'Tr' as input will generate 'Te'. In this case 'Ti' is needed as input as well.
%   Note 2: 'Te' will generate NaN...
% Input:
%   pars1d  = names of parameters wanted in par1D
%   pars2d  = names of parameters wanted in par2D
%   do_rpar = true/false (extract pp-data or not)
%   errs2d  = generated in function: those parameters in pars2d with corresponding variances in the
%     HDF5 file will be
% Default (load_param_hdf5(hdf5file)) sets
%   pars1d  = {'az' 'el' 'Pt' 'Tsys'};
%   par2d_id  = {'range' 'h' 'Ne' 'Tr' 'Ti' 'Vi' 'Collf' 'po+' 'res'};
%   do_rpars = true
%   errs2d_id  = {'var_Ne' 'var_Tr' 'var_Ti' 'var_Vi' 'var_Collf'}

% [Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('D:\EISCAT_data\analyzed\eiscat_website_hdf5\EISCAT_2019-11-20_folke_64@42mb.hdf5');
[Time,par2D,par1D,rpar2D,err2D,errs2d]=load_param_hdf5('C:\Svalbard\PhD\Data\Test_data\EISCAT_2019-11-20_folke_64@42mb.hdf5');

[y,m,d,h,mn,s]=datevec(Time(1,:));
ut_time=h+mn/60.;
alt=par2D(:,:,2);
ne=par2D(:,:,3);
te=par2D(:,:,4).*par2D(:,:,5);
ti=par2D(:,:,5);

n1 = datenum([2019 11 20 16 00 00]);
n2 = datenum([2019 11 20 20 00 00]);
xLimits=[n1 n2];
yLimits=[100 600];
f1=figure;
colormap jet;

subplot(3,1,1);
pcolor(datenum([y(:),m(:),d(:),h(:),mn(:),s(:)]),alt,log10(ne)),shading flat;
% give the limits of the colourbar
caxis([10 11]);
% % Add colorbar to the right of the plot
colorbar;
set(get(colorbar,'label'),'FontSize', 8,'string','Electron density (m^{-3})')
xlim(xLimits);
 ylim(yLimits) ;
datetick('x',13,'keeplimits')
% % Add labels to axes and colorbar
ylabel('Altitude [km]','FontSize', 8,'FontName','Arial')
 xlabel('Time, UT','FontSize', 8,'FontName','Arial')
 
