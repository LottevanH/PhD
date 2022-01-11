function [tl,indtst,indted]=ts_confinetimeline(tl0,tst,ted)

ind111=find((tl0-tst)<=0);
ind222=find((tl0-ted)>=0);
if isempty(ind111)
  indtst=1;
else
  indtst=ind111(end);
end
if isempty(ind222)
  indted=length(tl0);
else
  indted=ind222(1);
end
tl=tl0(1,indtst:indted);
end