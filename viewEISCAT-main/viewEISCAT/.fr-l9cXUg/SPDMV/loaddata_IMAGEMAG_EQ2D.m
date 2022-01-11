function [EQ1D]=IMAGEMAG_loaddata_EQ2D(fp_dat , fn_dat)

  global datasetinfo
  
  eval(['load ' [fp_dat fn_dat]]);
  
  lon0=20;
  if isfield(datasetinfo,'Jeq')
    if isfield(datasetinfo.Jeq,'lon')
      lon0=datasetinfo.Jeq.lon;
    end
  end
  datasetinfo.Jeq.lon=lon0;
  lat0=69.58;
  if isfield(datasetinfo,'Jeq')
    if isfield(datasetinfo.Jeq,'lat')
      lat0=datasetinfo.Jeq.lat;
    end
  end
  datasetinfo.Jeq.lat=lat0;
  
  lat_vec=Jlat(1,:);
  lon_vec=Jlon(:,1)';
  
  [minval, indlon]=min(abs(lon_vec-lon0));
  [minval, indlat]=min(abs(lat_vec-lat0));
  
  EQ1EArr=[];
  for i=1:length(dn)
    JY=reshape(JYext(i,:),size(Jlat));
    EQ1EArr=[EQ1EArr; JY(indlon,:)];
  end
  EQ1E.val=EQ1EArr';
  EQ1E.lat=lat_vec';
  EQ1E.lon=lon0;
  EQ1E.monosite.val=EQ1EArr(:,indlat)';
  EQ1E.monosite.lat=lat0;
  EQ1E.monosite.lon=lon0;
  EQ1E.tl=dn;
  
  EQ1NArr=[];
  for i=1:length(dn)
    JX=reshape(JXext(i,:),size(Jlat));
    EQ1NArr=[EQ1NArr; JX(indlon,:)];
  end
  EQ1N.val=EQ1NArr';
  EQ1N.lat=lat_vec';
  EQ1N.lon=lon0;
  EQ1N.monosite.val=EQ1EArr(:,indlat)';
  EQ1N.monosite.lat=lat0;
  EQ1N.monosite.lon=lon0;
  EQ1N.tl=dn;
  
  EQ1D.Elat=EQ1E;
  EQ1D.Nlat=EQ1N;
end
