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
    t_ran_dd_unit={'dd';'dd';'dd';'dd';'dd';'dd';'dd';'dd';'dd'};
    tkl_dd_unit={'HH';'HH';'HH';'dd';'dd';'dd';'dd';'dd';'dd'};
    mtk_dd_unit={'MM';'MM';'MM';'HH';'HH';'HH';'dd';'dd';'dd'};
    
    t_ran_HH=[1 2 3 6 12 24]/24;
    tkl_HH=[20/60 30/60 1 2 3]/24;
    tk_HH=[10/60 15/60 30/60 1 1]/24;
    mtk_HH=[3/60 5/60 10/60 20/60 30/60]/24;
    t_ran_HH_unit={'HH';'HH';'HH';'HH';'HH'};
    tkl_HH_unit={'MM';'MM';'HH';'HH';'HH'};
    mtk_HH_unit={'MM';'MM';'MM';'MM';'MM'};
    
    t_ran_MM=[1 2 3 5 10 15 20 30 60]/1440;
    tkl_MM=[20/60 30/60 1 2 3 5 5 10]/1440;
    tk_MM=[10/60 15/60 30/60 1 1 1 5 5]/1440;
    mtk_MM=[3/60 5/60 10/60 20/60 30/60 30/60 1 1]/1440;
    t_ran_MM_unit={'MM';'MM';'MM';'MM';'MM';'MM';'HH';'HH'};
    tkl_MM_unit={'SS';'SS';'MM';'MM';'MM';'MM';'MM';'MM'};
    mtk_MM_unit={'SS';'SS';'SS';'SS';'SS';'SS';'MM';'MM'};
    
    t_ran_SS=[1 2 3 5 10 15 20 30 60]/86400;
    tkl_SS=[20/60 30/60 1 2 3 5 5 10]/86400;
    tk_SS=[10/60 15/60 30/60 1 1 1 5 5]/86400;
    mtk_SS=[3/60 5/60 10/60 20/60 30/60 30/60 1 1]/86400;
    t_ran_SS_unit={'SS';'SS';'MM';'MM';'MM';'MM';'MM';'MM'};
    tkl_SS_unit={'FFF'; 'FFF'; 'SS'; 'SS';'SS';'SS';'SS';'SS'};
    mtk_SS_unit={'FFF'; 'FFF'; 'SS'; 'SS';'SS';'SS';'SS';'SS'};
    
    t_ran=[t_ran_SS(1:end-1) t_ran_MM(1:end-1) t_ran_HH(1:end-1) t_ran_dd];
    tkl=[tkl_SS tkl_MM tkl_HH tkl_dd];
    tk=[tk_SS tk_MM tk_HH tk_dd];
    mtk=[mtk_SS mtk_MM mtk_HH mtk_dd];
    
    t_ran_unit=[t_ran_SS_unit;t_ran_MM_unit;t_ran_HH_unit;t_ran_dd_unit];
    tkl_unit=[tkl_SS_unit;tkl_MM_unit;tkl_HH_unit;tkl_dd_unit];
    mtk_unit=[mtk_SS_unit;mtk_MM_unit;mtk_HH_unit;mtk_dd_unit];

    ind=find(t_ran>diff_t);
    ind=ind(1)-1;
    if ind==length(t_ran)
      disp('Timeline setting is out of range!')
      return
    end
    
    tunit=t_ran_unit{ind};
    
    xtlb_itv.val=tkl(ind);
    xtlb_itv.unit=tkl_unit{ind};
    xt_itv.val=tk(ind);
%     xtlb_itv.val=tkmd(ind);
%     xtlb_itv.unit=tkmd_unit{ind};
    mxt_itv.val=mtk(ind);
    mxt_itv.unit=mtk_unit{ind};
    
    dstr=[];
    if strcmp(tunit,xtlb_itv.unit) && strcmp(tunit,'dd')
      dstr='dd-mm-yyyy';
    else
      if strcmp(tunit,xtlb_itv.unit) 
        dstr=tunit;
      else
        if strcmp(tunit,'dd')
          dstr=xtlb_itv.unit;
        elseif strcmp(xtlb_itv.unit,'FFF')
          dstr=[tunit '.' xtlb_itv.unit];
        else
          dstr=[tunit ':' xtlb_itv.unit];
        end
      end
      if ~strcmp(xtlb_itv.unit,mxt_itv.unit)
        if strcmp(mxt_itv.unit,'FFF')
          dstr=[dstr '.' mxt_itv.unit];
        else
          dstr=[dstr ':' mxt_itv.unit];
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
  step_lb=round(xtlb_itv.val*86400);
  step_minor=round(mxt_itv.val*86400);
  
  st=sectime(1);
  ed=sectime(2);
  xt=floor(st/step_lb)*step_lb:step:ceil(ed/step_lb)*step_lb;
  
  lenxt=length(xt);
  xtlb=cell(lenxt,1);
  xt_minor=[];
  for j=1:lenxt
    if floor(t_ed)>floor(t_st) && ~strcmp(xtlb_itv.unit,'dd')
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
