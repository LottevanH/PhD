function [xt,xtlb,xt_minor]=getxticksettings_timeline(tl_ran,tunit,varargin)
  % timeline for days, hours, minites, seconds not for month and years
  % tl in MATLAB date format (datenum)
  % iterval unit: 'dd' - day; 'HH' - hour; 'MM' - minute; 'SS' - second 
  format long
  n_var=length(varargin);
  t_st=tl_ran(1);
  t_ed=tl_ran(2);
  diff_t=t_ed-t_st;
  if n_var==0
    
    t_ran_dd=[1 2 3 5 10 15 20 30 50 100];
    tkl_dd=[6/24 8/24 12/24 1 2 3 5 10 20];
    tk_dd= [3/24 4/24 6/24 12/24 1 1 5 5 10];
    mtk_dd=[1/24 2/24 3/24 6/24 8/24 12/24 1 1 2];
    tmd_dd_unit={'dd';'dd';'dd';'dd';'dd';'dd';'dd';'dd';'dd'};
    tkmd_dd_unit={'HH';'HH';'HH';'dd';'dd';'dd';'dd';'dd';'dd'};
    tkminormd_dd_unit={'MM';'MM';'MM';'HH';'HH';'HH';'dd';'dd';'dd'};
    
    tmd_HH=[1 2 3 6 12 24]/24;
    tkmd_HH=[15/60 30/60 1 2 3]/24;
    tkminormd_HH=[3/60 5/60 10/60 20/60 30/60]/24;
    tmd_HH_unit={'HH';'HH';'HH';'HH';'HH'};
    tkmd_HH_unit={'MM';'MM';'HH';'HH';'HH'};
    tkminormd_HH_unit={'MM';'MM';'MM';'MM';'MM'};
    
    tmd_MM=[1 2 3 5 10 15 20 30 60]/1440;
    tkmd_MM=[15/60 30/60 1 1 2 3 5 10]/1440;
    tkminormd_MM=[3/60 5/60 10/60 20/60 30/60 30/60 1 2]/1440;
    tmd_MM_unit={'MM';'MM';'MM';'MM';'MM';'MM';'HH';'HH'};
    tkmd_MM_unit={'SS';'SS';'MM';'MM';'MM';'MM';'MM';'MM'};
    tkminormd_MM_unit={'SS';'SS';'SS';'SS';'SS';'SS';'SS';'SS'};
    
    tmd_SS=[1 2 3 5 10 15 20 30 60]/86400;
    tkmd_SS=[15/60 30/60 1 1 2 3 5 10]/86400;
    tkminormd_SS=[3/60 5/60 10/60 20/60 30/60 30/60 1 2]/86400;
    tmd_SS_unit={'SS';'SS';'MM';'MM';'MM';'MM';'MM';'MM'};
    tkmd_SS_unit={'FFF'; 'FFF'; 'SS'; 'SS';'SS';'SS';'SS';'SS'};
    tkminormd_SS_unit={'FFF'; 'FFF'; 'SS'; 'SS';'SS';'SS';'SS';'SS'};
    
    tmd=[tmd_SS(1:end-1) tmd_MM(1:end-1) tmd_HH(1:end-1) tmd_dd];
    tkmd=[tkmd_SS tkmd_MM tkmd_HH tkmd_dd];
    tkminormd=[tkminormd_SS tkminormd_MM tkminormd_HH tkminormd_dd];
    
    tmd_unit=[tmd_SS_unit;tmd_MM_unit;tmd_HH_unit;tmd_dd_unit];
    tkmd_unit=[tkmd_SS_unit;tkmd_MM_unit;tkmd_HH_unit;tkmd_dd_unit];
    tkminormd_unit=[tkminormd_SS_unit;tkminormd_MM_unit;tkminormd_HH_unit;tkminormd_dd_unit];

    ind=find(tmd>=diff_t);
    ind=ind(1)-1;
    if ind==length(tmd)
      disp('Timeline setting is out of range!')
      return
    end
    
    tunit=tmd_unit{ind};
    
    xt_itv.val=tkmd(ind);
    xt_itv.unit=tkmd_unit{ind};
%     xtlb_itv.val=tkmd(ind);
%     xtlb_itv.unit=tkmd_unit{ind};
    xtminor_itv.val=tkminormd(ind);
    xtminor_itv.unit=tkminormd_unit{ind};
    
    dstr=[];
    if strcmp(tunit,xt_itv.unit) && strcmp(tunit,'dd')
      dstr='dd-mm-yyyy';
    else
      if strcmp(tunit,xt_itv.unit) 
        dstr=tunit;
      else
        if strcmp(tunit,'dd')
          dstr=xt_itv.unit;
        elseif strcmp(xt_itv.unit,'FFF')
          dstr=[tunit '.' xt_itv.unit];
        else
          dstr=[tunit ':' xt_itv.unit];
        end
      end
      if ~strcmp(xt_itv.unit,xtminor_itv.unit)
        if strcmp(xtminor_itv.unit,'FFF')
          dstr=[dstr '.' xtminor_itv.unit];
        else
          dstr=[dstr ':' xtminor_itv.unit];
        end
      end
    end
  elseif n_var==1;
    xt_itv=varargin{1};
    %xtlb_itv=xt_itv;
    xtminor_itv=[];
  elseif n_var==2;
    xt_itv=varargin{1};
    %xtlb_itv=varargin{2};
    xtminor_itv=varargin{2};
%   elseif n_var==3;
%     xt_itv=varargin{1};
%     %xtlb_itv=varargin{2};
%     xtminor_itv=varargin{3};
  end
  
  st_dd=floor(t_st);
  sectime=round((tl_ran-st_dd)*86400);
  step=round(xt_itv.val*86400);
  step_lb=step;
  step_minor=round(xtminor_itv.val*86400);
  
  st=sectime(1);
  ed=sectime(2);
  xt=floor(st/step)*step:step:ceil(ed/step)*step;
  
  lenxt=length(xt);
  xtlb=cell(lenxt,1);
  xt_minor=[];
  for j=1:lenxt
    if floor(t_ed)>floor(t_st) && ~strcmp(xt_itv.unit,'dd')
      if mod(xt(j),86400)==0
         xtlb{j}=datestr(st_dd+xt(j)/86400,'mmm dd');
         continue;
      end
    end
    if mod((xt(j)-xt(1)),step_lb)==0
      xtlb{j}=datestr(xt(j)/86400,dstr);
    else
      xtlb{j}=[];
    end
    if j<lenxt
      st1=xt(j); ed1=xt(j+1);
      xt_minor=[xt_minor st1+step_minor:step_minor:ed1-step_minor];
    end
  end
  xt=xt/86400+st_dd;
  xt_minor=xt_minor/86400+st_dd;
  
end
