%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% data manager for all types of data %%%
function [NaN]=SPDManager(varargin)
global datasetinfo drawopt dataset

if isempty(varargin)
  dataset=[];
  datasetinfo=[];
  % Project name which is shown in the data file name
  datasetinfo.projname='IMF/AE/IE';
  
  % Set parameters to be fetched (refer to initial files)
  datasetinfo.paralist={'Bz_gsm','By_gsm','AE','AL','AU','IE','IL','IU'};
  datasetinfo.save=0; % 0 - no save; 1 - save
  
  % Set plotting pannels
  drawopt.order={[1,-1;2,-2],[3;4;5],[6;7;8]};
  drawopt.plottype=nan(size(drawopt.order));
  drawopt.visual='on';
  drawopt.save=1;
  drawopt.figure.mode=1;
  
elseif length(varargin)==1
  datasetinfo=varargin{1};
  drawopt.visual='off';
elseif length(varargin)==1
  datasetinfo=varargin{1};
  drawopt=varargin{2};
end

% Set dates
datelist=getdatelist(1);

% Initialize database
para_init; 

% Fetch and plot
for i=1:size(datelist,1)
 
  
  datasetinfo.dateran=datelist(i,:);
  
  %% Read data
  for j=1:length(datasetinfo.filereader)
    datasetinfo.filereaderrec=j;
    [pathstr, name, ext] = fileparts(datasetinfo.filereader(j).name);
    eval(name)
  end

  %% Start plotting
  if strcmp(drawopt.visual,'off')
    continue
  else
    vis;
  end  
  
  %% Savedata
  if datasetinfo.save
    savedata;
  end
end
end

function [datelist]=getdatelist(id)

  if id==0
    prompt='Input starting date and time, e.g., [2005 09 01 00 00 00]:\n';
    dv=input(prompt);
    st_dn=datenum(dv);
    prompt='Input ending date and time, e.g., [2005 09 01 24 00 00]:\n';
    dv=input(prompt);
    ed_dn=datenum(dv);
  elseif id==1
    st_dn=datenum([2015 02 15 16 00 00]);
    ed_dn=datenum([2015 02 15 24 00 00]);
  elseif id==2
    dns=getdnlist(1); % dns in column
    st_time=[16 00 00]; % [HH MM SS];
    ed_time=[24 00 00];
    st_dn=dns+datenum([0 0 0 st_time]);
    ed_dn=dns+datenum([0 0 0 ed_time]);
  end

  datelist=[st_dn ed_dn];
end

function [dns]=getdnlist(id)
  if id==1
    dns=[datenum(2003,11,11):datenum(2003,11,19),        ...
      datenum(2005,09,01):datenum(2005,09,30)];
    dns=dns';
    dns=datenum([2011 12 17;
      2012 10 17;
      2012 11 25;
      2012 11 30;
      2012 12 01;
      2013 02 10;
      2013 03 08;
      2013 09 28;
      2015 02 15;
      2015 10 22;
      2016 10 20;
      2016 10 22;
      2016 10 23;
      2016 11 04;
      2016 11 06;
      2016 11 08;
      2016 11 16;
      2013 02 08;]);
  end
  
end
