function hal=axes_legend(hf,ha,legends,co,fz,varargin)
% set customed legends

layouts=[4 1];
if length(varargin)==1
  layouts=varargin{1};
end

set(0,'CurrentFigure',hf)
set(gcf,'currentaxes',ha)

pos_ha=get(gca,'position');

pos_hal=[pos_ha(1)+pos_ha(3)-0.12 pos_ha(2) 1-pos_ha(3)-pos_ha(1) pos_ha(4)];

hal=axes('position',pos_hal);

nl=length(legends);
m=layouts(1);
n=layouts(2);
hold on
for i=1:n
  for j=1:m
    if nl<j+(i-1)*m; continue; end
    xl=0.08+(i-1)/n;
    yl=0.9-(j-1)/m;
    
    xl1=xl;
    xl2=xl+0.4/n;
    plot([xl1 xl2],[yl yl],'-','color',co(j+(i-1)*m,:),'linewidth',1.5)
    text(xl2*1.08,yl,legends{j+(i-1)*m},'HorizontalAlignment','left',...
             'Color',co(j+(i-1)*m,:), ...
             'fontsize',fz,'VerticalAlignment','middle')
  end
end
set(hal,'xlim',[0 1],'ylim',[0 1],'visible','off')
set(gcf,'currentaxes',ha)

end
