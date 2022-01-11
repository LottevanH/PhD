function handle=my_colorbar(label,lg,fz,limval)

if nargin<2 | isempty(lg), lg='linear'; end
ax=[];

% Determine color limits by context.  If any axes child is an image
% use scale based on size of colormap, otherwise use current CAXIS.

ch=get(gca,'children');
hasimage=0; t=[];
for i=1:length(ch)
  if strcmp(get(ch(i),'type'),'image')
    hasimage=1;
    t=get(ch(i),'UserData'); % Info stored by imshow or imagesc
  elseif strcmp(get(ch(i),'Type'),'surface'), % Texturemapped surf?
    if strcmp(get(ch(i),'FaceColor'),'texturemap')
      hasimage=2;
      t=get(ch(i),'UserData'); % Info stored by imshow or imagesc
    end
  end
end
if hasimage
  if isempty(t), t=[0.5 size(colormap,1)+0.5]; end
else
  t=caxis;
end

h=gca;

% Search for existing colorbar
ch=get(gcf,'children'); ax=[];
for i=1:length(ch),
  d=get(ch(i),'userdata');
  if prod(size(d))==1 & d==h
    ax=ch(i); break
  end
end

if strcmp(get(gcf,'NextPlot'),'replace'),
  set(gcf,'NextPlot','add')
end


if isempty(ax)
  pos=get(h,'Position');
% stripe=0.03; edge=0.08; 
% [az,el]=view;
% if all([az,el]==[0 90]), space=0.02; else, space=.1; end
% set(h,'Position',[pos(1) pos(2) pos(3)*(1-stripe-edge-space) pos(4)])
% rect=[pos(1)+(1-stripe-edge)*pos(3) pos(2) stripe*pos(3) pos(4)];
  rect=[pos(1)+1.03*pos(3) pos(2) .045*pos(3) pos(4)-0.02];

  % Create axes for stripe
  ax=axes('Position', rect);
else
  axes(ax)
end

% Create color stripe
n=size(colormap,1);
if strcmp(lg,'log')
  tt=((0:n)-.5)'*diff(t)/(n-1)+t(1);
  surface([0 1],tt,[tt tt])
  set(ax,'CLim',t,'ylim',t,'tickdir','in','ytick',[],'xtick',[],'xlim',[0 1])
  freezeColors;
  ax2=axes('Position', rect);
  set(ax2,'xlim',[0 1],'ylim',t,'visible','off')
  if diff(t)<=2
    stept=1;
  else
    stept=2;
  end
  ll=[];yy=[];

  if isempty(limval)
    hold on
    for ii=floor(t(1)):floor(t(2))+1
      ll=[ll 0.25*(0.2:0.1:1)];
      yy=[yy log10(2:10)+ii];
      plot([0 1],[ii ii],'k:','linewidth',1)
    end
    xt1=[zeros(1,length(ll)); ll; nan(1,length(ll))];
    xt2=[ones(1,length(ll)); 1-ll; nan(1,length(ll))];
    yt=[yy; yy; nan(1,length(ll))];
    plot(xt1(:),yt(:),'k-')
    plot(xt2(:),yt(:),'k-')
    hold off
    for ii=ceil(t(1)):stept:floor(t(2))
      yticklabel=['10^{' num2str(ii) '}'];
      text(1.2,ii,yticklabel,'verticalalignment','baseline', ...
           'horizontalalignment','left','fontsize',fz-4)
    end
  else    
    ll=[]; yy=[];
    hold on
    for ii=0:floor(t(1))
      ll=[ll 0.25*(0.2:0.1:1)];
      yy=[yy ii-log10(2:10)];
      plot([0 1],[ii ii],'k:','linewidth',1)
    end
    for ii=0:ceil(t(2))
      ll=[ll 0.25*(0.2:0.1:1)];
      yy=[yy log10(2:10)+ii];
      plot([0 1],[ii ii],'k:','linewidth',1)
    end
    xt1=[zeros(1,length(ll)); ll; nan(1,length(ll))];
    xt2=[ones(1,length(ll)); 1-ll; nan(1,length(ll))];
    yt=[yy; yy; nan(1,length(ll))];
    plot(xt1(:),yt(:),'k-')
    plot(xt2(:),yt(:),'k-')
    hold off
    ytt=[fliplr(0:-stept:ceil(t(1))) stept:stept:floor(t(2))];
    for ii=ytt(1):stept:2
      if ii==0
        yticklabel=['\pm 10^{' num2str(limval) '}'];
      elseif ii<0
        yticklabel=['-10^{' num2str(limval-ii) '}'];
      elseif ii>0
        yticklabel=['+10^{' num2str(limval+ii) '}'];
      end
      text(1.2,ii,yticklabel,'verticalalignment','middle', ...
           'horizontalalignment','left','fontsize',fz-4)
    end
  end
else
  image([0 1],t,[1:n]')
% set(ax,'TickDir','in')
  set(ax,'ylim',t,'Ydir','normal')
  if abs(t(1)+180)<1
    set(ax,'ytick',[-180:45:180])
  end
  xplus=2;
% Justify color axis
% set(ax,'yticklabel',strjust(get(ax,'yticklabel')))
end
axes_colorfreezing;
set(ax,'userdata',h,'YaxisLoc','right','xtick',[],'fontsize',fz-2,'tickdir','in')


% ylabel(label,'Rotation',-90,'VerticalAlignment','baseline','fontsize',fz)
% hyl=get(ax,'Ylabel');
% set(hyl,'position',get(hyl,'position')+[xplus 0 0 ])
% if length(label)==2
%   set(hyl,'position',get(hyl,'position')+[0.47 0 0 ])
% end
set(gcf,'CurrentAxes',h)
if length(label)==2
  fz_add=-2;
else
  fz_add=0;
end
text('String',label,'Rotation',-90,'VerticalAlignment','baseline',...
      'fontsize',fz+fz_add,'position',[1.18 0.5],'Units','normalized',...
      'HorizontalAlignment','center');


set(gcf,'Nextplot','Replace')

if nargout>0, handle=ax; end
return
