%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% data manager for all types of data %%%
function [NaN]=viewEISCAT(varargin)
if isempty(varargin)
  clear global
end
global datasetinfo drawopt dataset


addpath_before_start;

if isempty(varargin)
  dataset=[];
  datasetinfo=[];
  
  % Project name which is shown in the data file name
  datasetinfo.projname='EISCAT/Swarm';
  
  % Set EISCAT experiment time range
  dt_fr = datenum([2020 11 21 22 00 00]);
  dt_to = datenum([2020 11 22 01 00 00]);
  datelist=getdatelist(1, dt_fr, dt_to);
  
  % Set Swarm trajectory parameters
  swarm_id = "A";
  dt_swarm = datenum([2020 11 22 00 27 00]);
  add_swarm_overpass(dt_swarm, dt_fr, swarm_id);

  
  swarm_id = "B";
  dt_swarm = datenum([2020 11 21 21 29 21]);
  add_swarm_overpass(dt_swarm, dt_fr, swarm_id);
    
  
  % set the mode to select the analyzed results: "manual" - a dialog box
  % will open and select the data folder manually; "autoselect": search 
  % /kaappi/EISCAT/RESULTS and select the data folder automatically.
  datasetinfo.EISCAT.filemode = "manual";
  
  datasetinfo.EISCAT.antenna = 'UHF'; % if filemode is "autoselect" must be specified.
  datasetinfo.EISCAT.pulsecode = 'beata'; % if filemode is "autoselect" must be specified.
  
  % Set parameters to be retrieved (refer to initial files)
  datasetinfo.paralist={'ne_lv0', 'Te_lv0', 'Ti_lv0', 'vi_lv0'};

  datasetinfo.save=0; % 0 - no save; 1 - save retrieved data
  
  % Set plotting pannels
  drawopt.order={1, 2, 3, 4}; % order of the panels to show in the figure for the parameters in datasetinfo.paralist
  drawopt.plottype=nan(size(drawopt.order));
  drawopt.visual='on';
  drawopt.figureclose=0; % 0 - keep fig in screen; 1- close figure after saving
  drawopt.figure.position=[1 1 20 25];
  drawopt.FigureFormat = {'png'};
  drawopt.save=1;
  % drawopt.addlines=getadditionallines(1);
  
  drawopt.xunit='UT';
elseif length(varargin)==1
  datasetinfo=varargin{1};
  drawopt.visual='off';
elseif length(varargin)==1
  datasetinfo=varargin{1};
  drawopt=varargin{2};
end


% Fetch and plot
for i=1:size(datelist,1)
 
  para_init; 
  
  
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


function [datelist]=getdatelist(id, varargin)

  if id==0
    prompt='Input starting date and time, e.g., [2005 09 01 00 00 00]:\n';
    dv=input(prompt);
    st_dn=datenum(dv);

    prompt='Input ending date and time, e.g., [2005 09 01 24 00 00]:\n';
    dv=input(prompt);
    ed_dn=datenum(dv);
  elseif id==1
    if isempty(varargin)
        st_dn=datenum([2021 02 09 00 00 00]);
        ed_dn=datenum([2021 02 09 04 00 00]);
    else
        st_dn = varargin{1};
        ed_dn = varargin{2};
    end
  elseif id==2
    dns=getdnlist(2); % dns in column
    st_time=[00 00 00 16 00 00]; % [HH MM SS];
    ed_time=[00 00 01 02 00 00];
    st_dn=dns+datenum(st_time);
    ed_dn=dns+datenum(ed_time);
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
    elseif id==2
    fn=[pwd '/events/' 'event_dates1.dat'];
    fid=fopen(fn);
    C=textscan(fid,'%s');
    fclose(fid);
    dns=datenum(C{1},'yyyymmdd');
  end
  
end

function [] = add_swarm_overpass(dt, dt_fr, id, varargin)
    global drawopt
    
    if ~isfield(drawopt, 'addlines')
        drawopt.addlines = {};
    end

    if strcmp(id.extract(1), 'B')
        color = [0.145, 0.615, 0.180];
    elseif strcmp(id.extract(1), 'A') || strcmp(id.extract(1), 'C')
        color = [0.615, 0.145, 0.486];
    end
    linewidth = 2;
    if ~isempty(varargin)
        if length(varargin) >= 1
            color = varargin{1};
        elseif length(varargin) >= 2
            linewidth = 2;
        end
    end
    
    
    h24 = (dt -floor(dt_fr))*24;
    drawopt.addlines = [drawopt.addlines;   ...
        {0,[h24 h24],NaN,   ...
        struct('LineStyle', '--', 'LineWidth', 2, 'Color',color, 'message', id)}];
end

function lines=getadditionallines(varargin)
  lines={1,NaN,[69.58 69.58],struct('LineStyle',':');   ...
         2,NaN,[69.58 69.58],struct('LineStyle',':','Color',[0 0 0]);   ...
         0,[2.3 2.3],NaN,struct('Color',[0.7 0 0]);   ...
         0,[23.46 23.46],NaN,struct('Color',[0.7 0 0]);   ...
         0,[22.79 22.79],NaN,[];   ...
         0,[23.57 23.57],NaN,[];       ...
         0,[22.14 22.14],NaN,[]};   
end

function shadings=getadditionalshadings(varargin)
  dn=datenum([2012 12 17 22 22 00]);
  t1=(dn-floor(dn))*24;
%   shadings={1,NaN,[69.58 69.58];   ...
%             0,[t1 t1],NaN};
end

