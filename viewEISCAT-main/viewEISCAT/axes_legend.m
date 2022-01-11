function hal=axes_legend(hf,ha, hps, alphas, legends,co,fz,varargin)
% set customed legends

layouts=[4 2];
if length(varargin)==1
  layouts=varargin{1};
end

set(0,'CurrentFigure',hf)
set(gcf,'currentaxes',ha)

pos_ha=get(gca,'position');

pos_hal=[pos_ha(1)+pos_ha(3) pos_ha(2) 1-pos_ha(3)-pos_ha(1) pos_ha(4)];

hal=axes('position',pos_hal);

nl=length(legends);
m=layouts(1);
n=layouts(2);
hold on
for i=1:n
  for j=1:m
    if nl<j+(i-1)*m; continue; end
    xl=0.08+(i-1)/n*0.9;
    yl=0.9-(j-1)/m*0.9;
    
    xl1=xl;
    xl2=xl+0.2/n;
    np = j+(i-1)*m;
    % np = i + (j-1)*n;
    % plot([xl1 xl2],[yl yl], '-','color',co(j+(i-1)*m,:),'linewidth',1.5)
    hold on
    hp1 = plot([xl1 xl2], [yl yl], 'LineStyle', hps(np).LineStyle, ...
        'LineWidth', hps(np).LineWidth, 'Color', hps(np).Color);
    hp1.Color(4) = alphas(np);
    plot( (xl1+xl2)/2 , yl, 'Marker', hps(np).Marker, ...
        'MarkerSize',  hps(np).MarkerSize, 'Color', hps(np).Color)
    hold off
    text(xl2*1.08,yl,legends{j+(i-1)*m},'HorizontalAlignment','left',...
             'Color',hps(np).Color, ...
             'fontsize',fz,'VerticalAlignment','middle')
  end
end
set(hal,'xlim',[0 1],'ylim',[0 1],'visible','off')
set(gcf,'currentaxes',ha)

end
