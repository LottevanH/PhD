function a=cal_FPIa(v,t,varargin)
  if isempty(varargin)
    N=3;
  else
    N=varargin{1};
  end
  aval=[];
  for i=1:length(v)-N+1
    x=t(i:i+N-1);
    y=v(i:i+N-1);
    p=polyfit(x,y,1);
    aval(i)=p(1);
    if mod(N,2)
      ta(i)=x(round(N/2));
    else
      ta(i)=(x(N/2)+x(N/2+1))/2;
    end
  end
  a.val=aval;
  a.tl=ta;
  a.err=[];
end