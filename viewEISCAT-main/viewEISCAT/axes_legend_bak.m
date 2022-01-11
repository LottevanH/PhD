function ax=axes_legend(hl,Text,ls,pos,co,fz,mode)
h=gca;
if isempty(hl)
  pos_h=get(h,'position');
  rect=[pos_h(1)+1.02*pos_h(3) pos_h(2) .05*pos_h(3) pos_h(4)];
  ax=axes('Position', rect,'Ylim',[0 1],'Xlim',[0 1],'visible','off');
  
else
  ax=hl;
end
set(gcf,'currentaxes',ax)
hold on
x=[0 1];
y=[pos(2) pos(2)];
plot(x,y,ls,'color',co,'linewidth',1.5)
hold off
text('String',Text,'Units','Normalized','HorizontalAlignment','left',...
             'Position',[1.08 pos(2)],'Color',co, ...
             'fontsize',fz,'VerticalAlignment','middle')
set(gcf,'currentaxes',h)
end
