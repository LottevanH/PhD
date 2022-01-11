function [cmap]=vis_colormaps(go)
  addpath([pwd '/colormaps/']);
  if ~go
    cmap{1}='jet';
    %cmap{1}=pmkmp(256,'CubicYF');
    %cmap{1}=colormap_mypalette(1,'whole',256,1);
%     fp=[pwd '/colormaps/extended_kindlmann/'];
%     fn='extended-kindlmann-table-float-0512.csv';
%     dat=csvread([fp fn],1,0);
%     cmap{1}=dat(:,2:4);
    
%     fp=[pwd '/colormaps/extended_blackbody/'];
%     fn='extended-black-body-table-float-0512.csv';
%     dat=csvread([fp fn],1,0);
%     cmap{1}=dat(:,2:4);
    
    
  else
    cmap{1}=myb(1,'whole',80,1);
    cmap{2}=myb(7,'whole',120,1); % wind magnitude
    cmap{3}=myb(4,'whole',100,0); % neutral winds
    cmap{4}=myb(6,'whole',180,0); % wind clock angles?
  end
end