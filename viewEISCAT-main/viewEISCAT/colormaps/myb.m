function [ f ]=myb(num,optc,nl,cut)
% create palette with nl levels
if nargin<2, cut=[]; end
if nargin==0, nl=[]; end
if isempty(cut), cut=0; end
if isempty(nl)
 nl=size(get(gcf,'colormap'),1);
end
funcn={'KBGYRM','DarkKBGYRM','GreenRed36','WarmCool','ColorWind','ColorWheel','WarmCool2'};

f=eval([funcn{num} '(optc)']);
nc=size(f,1);
n=nc-cut;
b=round([0:n-1]/(n-1)*(nl-1))+1;
%f=sin(interp1(b,f(1:n,:),1:nl)*pi/2);
f=tanh(interp1(b,f(1:n,:),1:nl))/tanh(1);
%f=interp1(b,f(1:n,:),1:nl);
return

end

