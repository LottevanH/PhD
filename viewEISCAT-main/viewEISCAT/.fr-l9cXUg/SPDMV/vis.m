function []=vis()
global drawopt dataset datasetinfo
fp_res=[pwd '/results/'];

[nrow,ncol]=size(drawopt.order);
npnl=nrow*ncol;
figure_layout;
axes_layout;

unifig=drawopt.figure.unit;
posfig=drawopt.figure.position;
paperpos=drawopt.figure.paperpos;
paperort=drawopt.figure.paperort;
figmode=drawopt.figure.mode;
fontsize_axes=drawopt.figure.axesfontsize;
fontsize_text=drawopt.figure.textfontsize;

ax=drawopt.axes.ax;
ay=drawopt.axes.ay;
ahei=drawopt.axes.ahei;
awid=drawopt.axes.awid;

hf=figure;
set(gcf,'Units',unifig,'Position',posfig,'DefaultAxesFontSize',fontsize_axes,...
  'DefaultAxesTickDir','in','DefaultTextFontSize',fontsize_text,'UserData',[],...
  'DefaultAxesXGrid','on','DefaultAxesYGrid','off','DefaultAxesYMinorTick','on',...
  'DefaultAxesXMinorTick','off','defaultaxesbox','on',...
  'DefaultAxesTickLength',[0.02 0.035],...
  'PaperType','A4','PaperOrientation',paperort,...
  'PaperUnits','centimeters', 'PaperPosition',paperpos,...
  'DefaultAxeslayer','top','DefaultsurfaceEdgeColor','none',...
  'DefaultTextHorizontalAlignment','center')

lsp=drawopt.axes.lsp;
ahei_tot=0;
for i=1:npnl
  ha(i)=axes('unit','normalized','position',[ax(i) ay(i) awid(i) ahei(i)]);
  ahei_tot=ahei_tot+ahei(i);
end

t_st=datasetinfo.dateran(1);
t_ed=datasetinfo.dateran(2);

for i=1:npnl
  
  [xtick,xticklabel,xminortick]=axes_xtick_timeline2([t_st t_ed],'HH');
  
  xtick=(xtick-floor(t_st))*24;
  xminortick=(xminortick-floor(t_st))*24;
  % settings 1
  set(gcf,'currentaxes',ha(i));
  if i==npnl && figmode==1 || figmode==0
    xticklabel0=xticklabel;
  elseif i<npnl && figmode==1
    xticklabel0=[];
  end
  
  hl=[];
  orderarr=drawopt.order{i};
  [m n]=size(orderarr);
  co=axes_legend_linecolor(m);
  
  plottype=drawopt.plottype(i);
    
  Ylim=[];
  %% plottype=1 : line plots
  for d=1
  if plottype==1
    for j=1:m
      para=dataset{orderarr(j,1)};
      
      if isempty(para.drawopt.xdataname)
        tl=para.tl;
      else
        tl=eval(['para.' para.drawopt.xdataname ';']);
      end
      
      
      data=para.val;
      if n==2
        err=para.err;
        if isempty(err)
          err=zeros(size(data));
        end
      else
        err=zeros(size(data));
      end
      
%       ind111=find((tl-t_st)<=0);
%       ind222=find((tl-t_ed)>=0);
%       if isempty(ind111)
%         ind_tst=1;
%       else
%         ind_tst=ind111(end);
%       end
%       if isempty(ind222)
%         ind_ted=length(tl);
%       else
%         ind_ted=ind222(1);
%       end
      ind_tst=1;ind_ted=length(tl);
      
      ydata=data(ind_tst:ind_ted); yerr=err(ind_tst:ind_ted);
      xdata=tl(ind_tst:ind_ted);
      ndata=length(xdata);
      
      xdata=(xdata-floor(t_st))*24;
      % build errorbar
      xx=reshape([xdata; xdata; nan(1,ndata)],3*ndata,1);
      yy=reshape([ydata-yerr; ydata+yerr; nan(1,ndata)],3*ndata,1);
      % plot data and errorbars
            
      hold on
      plot(xx,yy,'-','linewidth',1,'color',co(j,:)) % plot error bar first
      plot(xdata,ydata,lsp{j},'linewidth',1.5,'color',co(j,:))
      
      hold off
      
      if ~isempty(xdata)
        Ylim=[min([Ylim nanmin(ydata-yerr)]) nanmax([Ylim max(ydata+yerr)])];
        if isempty(find(isfinite(Ylim)==1))
          Ylim=[-1 1];
        end
      else
        Ylim=[-1 1];
      end
      userdata=get(gca,'UserData');
      if isempty(userdata)
        userdata=struct('xdata',{tl},'ydata',{ydata},'yerr',{yerr},     ...
          'datasetname',{['dataset{' num2str(orderarr(j,1)) '}']});
      else
        userdata.xdata=[userdata.xdata; {tl}];
        userdata.ydata=[userdata.ydata; {ydata}];
        userdata.yerr=[userdata.yerr; {yerr}];
        userdata.datasetname={userdata.datasetname;     ...
          {['dataset{' num2str(orderarr(j,1)) '}']}};
      end
      set(gca,'UserData',userdata)      
      
      % show legend
      if m>1
        p=floor(m/4)+1;
        k=floor(j/4)+1;
        pos_text=[1.005 1.2-(k-1)*0.125-j*.25+(p-2)*0.1];
        co_text=co(j,:);
        fz_text=fontsize_text-2*p-figmode*2;
        zlabel=para.label;
        hl=axes_legend(hl,zlabel,lsp{j},pos_text,co_text,fz_text,figmode);
      end
    end
    
    % xlim
    Xlim=[t_st t_ed];
    Xlim=(Xlim-floor(t_st))*24;
  end
  end
  %% plottype=2 : 2D plot
  for d=1
  if plottype==2
    % plot values or errors
    if orderarr<0
      para=dataset{-orderarr};
      data=para.err;
    else
      para=dataset{orderarr};
      data=para.val;
    end
    
    % some settings
    cmap=para.drawopt.colormap;
    Ylim=para.drawopt.ylim;
    Ylabel=para.drawopt.ylabel;
    zlim=para.drawopt.zlim;
    zscale=para.drawopt.zscale;
    zlabel=para.drawopt.zlabel;
    limval=[];
    if strcmp(zscale,'log')
      if ~isreal(zlim)
        limval= dataset{ind_para}.drawopt.limval;
        zlim(1)=limval-imag(zlim(1));
        zlim(2)=zlim(2)-limval;
        ind_0= find(data<=10^limval & data>=-10^limval);
        ind_nan=find(data==0 | isnan(data));
        ind_p=find(data>10^limval);
        ind_n=find(data<-10^limval);
        data(ind_0)=0;
        data(ind_p)=log10(data(ind_p))-limval;
        data(ind_n)=limval-log10(-data(ind_n));
        data(ind_nan)=NaN;
        rate=abs(zlim(1)/zlim(2));
        cmap=cat(1,colormap_mypalette(2,'negative',round(80*rate),0),cmap);
      else
        data=log10(data);
      end
    end
    zdata=data;
    
    if isempty(para.drawopt.xdataname)
      tl=para.tl;
    else
      tl=eval(['para.' para.drawopt.xdataname ';']);
    end
    
    if isstruct(tl)
      t1=para.tl.t1;
      t2=para.tl.t2;
      tl=para.tl.tmid;
    else
      t1=[];
      t2=[];
    end
    
    xdata=tl;
    if isempty(para.drawopt.ydataname)
      ydata=para.alt;
    else
      ydata=eval(['para.' para.drawopt.ydataname ';']);
    end
    Ylim=[ydata(1) ydata(end)];
    
%     if strcmp(para.drawopt.xdatagrid,'on')
   %     xres=para.drawopt.xdatagridres;
%     yres=para.drawopt.ydatagridres;
%     
%     if isempty(xres)
%       xres=median(diff(xdata));
%     end
%     if isempty(yres)
%       yres=median(diff(ydata));
%     end 
%     
%     [xdata, ydata, zdata]=vis_gridzdata(xdata,ydata,zdata,t_st,t_ed,xres,yres);
%     
%     imagesc()
%     end

    xres=para.drawopt.xdatares;
    yres=para.drawopt.ydatares;
    
    if isempty(xres)
      xdiff=diff(xdata);
      xres=median(xdiff);
    end
    if isempty(yres)
      ydiff=diff(ydata);
      yres=median(ydiff);
    end 
    
    [nx, ny]=size(zdata);
    if ~isempty(t1) && ~isempty(t2) && para.drawopt.xdatares==-1       
      zzz=[zdata;zdata;nan(nx,ny)];
      zdata=reshape(zzz(:),ny,3*nx);
      tadd=[t1(2:end)-t2(1:end-1) t2(end)+xres]/2;
      xxx=[t1;t2; t2+tadd];
      xdata=xxx(:);
      xres=median(diff(xdata));
    else      
      data=zdata;
      colvec=xdata;
      coldif=xdiff;
      colres=xres;
      indgaps=find(coldif>1.5*colres); % define the data gap
      for ng=1:length(indgaps)
        indg=indgaps(ng);
        addval=[data(:,indg) nan(ny,1)];
        [data,xdata]=array_addnewcolumn(data,colvec,indg,addval);
      end
    end
    
    zdata=[zdata, zdata(:,end)];
    zdata=[zdata; zdata(end,:)];
    xdata=[xdata xdata(end)+xres];
    ydata=[ydata; ydata(end)+yres];   
    
    xdata=xdata-xres/2;
    ydata=ydata-yres/2;
    
    xdata=(xdata-floor(t_st))*24;
    
    surface(xdata,ydata,zdata); shading flat;
    % xlim
    Xlim=([t_st t_ed]-floor(t_st))*24;
    % zlim
    zlim=axes_axis_limits([nanmin(data(:)) nanmax(data(:))],zlim);
    set(gca,'clim',zlim);
    % zlabel
    if isempty(zlabel)
      zlabel=axes_axis_label(para.label,para.unit);
    end
    
    %ylabel(Ylabel);
    % colorbar
    colormap(cmap); %
    axes_colorfreezing;
    axes_colorbar(zlabel,zscale,fontsize_text-2,limval);
    
  end
  end

  %% Plottype=-1 : vectors as a fuction of the timeline 
  for d=1
  if plottype==-1
    datascale=5;
    
    indvec=1:length(tl);
    x=dataset{paraind(1)}.tl(indvec);
    y=zeros(size(x));
    
    xlim=[t_st t_ed];
    ylim=[-1 1];
    
    plot_scale=diff(xlim)/15;
    fwid=posfig(3);fhei=posfig(4);
    
    udata=dataset{paraind(1)}.val(indvec);
    vdata=dataset{paraind(2)}.val(indvec);
    u=udata/datascale*plot_scale;
    v=vdata/datascale*plot_scale*diff(ylim)/diff(xlim)   ...
      *fwid/fhei*awid(i)/ahei(i);
       
    hq=quiver(x,y,u,v,0,'k-');
    %set(hq,'MaxHeadSize',1)
    set(hq,'ShowArrowHead','off')
    hold on
    plot(x,y,'k:')
    hold off
    
    % ylabel
%     text('String',{'IMF', 'GSM y-z plane'},'unit','normalized','position',[-0.1 0.5],   ...
%       'Rotation',90,'VerticalAlignment','middle','fontsize',10)
    set(gca,'xlim',xlim,'ylim',ylim,	...
      'xtick',xtick,'xticklabel',xticklabel0,   ...
      'ytick',ylim,'yticklabel',{},'yminortick','off')
    % draw legend
    ha_sp=axes('units','normalized','position',[ax(i)+awid(i) ay(i) awid(i) ahei(i)]);
    x=xlim(1)+0.08*diff(xlim);
    y=ylim(2)-0.3*diff(ylim);
    u=1*plot_scale;
    v=0;
    hq=quiver(x,y,u,v,0,'k-');
    set(hq,'ShowArrowHead','off')
%     text(x+plot_scale/2,y-0.15*diff(ylim),['B=' sprintf('%3.1f',Bt_scale) ' nT'],  ...
%       'fontsize',10,'fontweight','normal')
    set(gca,'xlim',xlim,'ylim',ylim,'visible','off')
  end
  end
  %% set axes and titles for line plots
  
  % ylim
  if ~isempty(para.drawopt.ylim)
    Ylim=axes_axis_limits(Ylim,para.drawopt.ylim);
  end
  % yscale
  if isempty(para.drawopt.yscale)
    Yscale='linear';
  else
    Yscale=para.drawopt.yscale;
  end
  % add zero line
  if Ylim(1)<0 && Ylim(2)>0
    hold on
    plot(Xlim,[0 0],'k--')
    hold off
  end
  
  % set gca
  set(gca,                              ...
        'xlim',         Xlim,           ...
        'ylim',         Ylim,           ...
        'yscale',       Yscale          ...
        );
  % ytick after ylim is set
  if isempty(para.drawopt.ytick)
    Ytick=get(gca,'ytick');
  else
    Ytick=para.drawopt.ytick;
  end
  % ytick label
  if isempty(para.drawopt.yticklabel)
    Yticklabel=get(gca,'yticklabel');
  else
    Yticklabel=para.drawopt.yticklabel;
  end
  % generate y label
  if plottype==1
  if m>1
    Ylabel=axes_axis_label(para.group,para.unit);
  else
    Ylabel=axes_axis_label(para.label,para.unit);
  end
  elseif plottype==2
    Ylabel=axes_axis_label(para.drawopt.ylabel,para.drawopt.yunit);
  end
  % set axis
  set(gca,                              ...
        'xtick',        xtick,          ...
        'xticklabel',   xticklabel0,    ... 
        'ytick',        Ytick,          ...
        'yticklabel',   Yticklabel      ...
        );
  
  % add y label
  haa=gca;
  ha1(i)=axes('unit','normalized','position',[ax(i) ay(i) awid(i) ahei(i)]);
  text('String',Ylabel,'unit','normalized','position',[-0.12 0.5],   ...
    'Rotation',90,'VerticalAlignment','baseline','fontsize',fontsize_text)
  set(gca,'visible','off')
  set(gcf,'currentaxes',haa)
  
  % add x minor ticks
  haa=gca;  
  tllen=get(gca,'ticklength');
  Ylim=get(gca,'Ylim');
  mtk_scale=tllen(1);
  xxx=reshape([xminortick;xminortick;nan(size(xminortick));   ...
    xminortick;xminortick;nan(size(xminortick));], ...
    length(xminortick)*6,1);
  ymtklen=ones(size(xminortick))*diff(Ylim)*mtk_scale;
  yyy=reshape([Ylim(2)*ones(size(ymtklen));Ylim(2)-ymtklen;nan(size(ymtklen));    ...
    Ylim(1)+ymtklen;Ylim(1)*ones(size(ymtklen));nan(size(ymtklen))],   ...
    length(xminortick)*6,1);
  ha1(i)=axes('unit','normalized','position',[ax(i) ay(i) awid(i) ahei(i)]);
  hold on
  plot(xxx,yyy,'k-','linewidth',0.5)
  set(gca,                              ...
        'xlim',         Xlim,           ...
        'ylim',         Ylim,           ...
        'xtick',        Xlim,           ...
        'ytick',        Ylim)
  hold off
  set(gca,'visible','off')
  set(gcf,'currentaxes',haa)
  
end
%% add title
ha_title=axes('unit','normalized','position',[ax(1) ay(1) awid(1) ahei(1)],'visible','off');
if floor(t_st)==floor(t_ed)
  Title=[datasetinfo.projname ': ' datestr(t_st,'yyyy-mm-dd HH:MM') 	...
    '-' datestr(t_ed,'HH:MM')];
else
  Title=[datasetinfo.projname ': ' datestr(t_st,'yyyy-mm-dd HH:MM') 	...
    ' - ' datestr(t_ed,'yyyy-mm-dd HH:MM')];
end
text('String',Title,'unit','normalized','position',[1.00 1.02],		...
  'HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',fontsize_axes);
set(gcf,'currentaxes',ha(1))
%% save and print
if ~isfield(drawopt,'save')
  drawopt.save=0;
end
if ~isfield(drawopt,'figureclose')
  drawopt.figureclose=0;
end
if drawopt.save==1
  fp_res=[pwd '/results/'];
  prefn=strrep(datasetinfo.projname,'/','-');
  if floor(t_st)==floor(t_ed)
    fn=[prefn '_' datestr(t_st,'yyyymmdd_HHMM') 	...
      '-' datestr(t_ed,'HHMM')];
  else
    fn=[prefn '_' datestr(t_st,'yyyymmdd_HHMM') 	...
      '-' datestr(t_ed,'yyyymmdd_HHMM')];
  end
  
  saveas(hf,[fp_res fn ],'fig');
  %print('-painters','-dtiff','-r300',[fp_res fn])
  print('-dpng','-r500',[fp_res fn]);
end
if drawopt.figureclose==1
  close(hf)
end
end
