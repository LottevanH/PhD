function [X,Y,Z]=vis_gridzdata(x,y,z,Xstep,Ystep,interpmethod,varargin)
  stepscale=6;
  if Xstep==0
    X=x;
  else
    xdiff=diff(x);
    xres=median(xdiff);
    if isempty(Xstep)  
      Xstep=xres/stepscale;
    end
    X=x(1):Xstep:x(end);
  end
  
  if Ystep==0
    Y=y;
  else
    ydiff=diff(y);
    yres=median(ydiff);
    if isempty(Xstep)  
      Xstep=yres/stepscale;
    end
    Y=y(1):Xstep:y(end);
  end
  
  Z=griddata(x,y,z,X,Y,interpmethod);
  
  if Xstep~=0
    indgap=find(xdiff>1.5*xres);
    for i=1:length(indgap)
      indnan=find(X>x(indgap(i)) & X<x(indgap(i)+1));
      Z(indnan,:)=NaN;
    end
  end
  
  if Ystep~=0
    indgap=find(ydiff>1.5*yres);
    for i=1:length(indgap)
      indnan=find(Y>y(indgap(i)) & Y<y(indgap(i)+1));
      Z(:,indnan)=NaN;
    end
  end
  
  if ~isempty(varargin)
    x1=varargin{1};
    x2=varargin{2};
    
    indnan=[];
    for i=1:length(x1)-1
      ind1=find(X>x2(i) & X<x2(i-1));
      indnan=[indnan ind1];
    end
    Z(indnan,:)=NaN;
  end
end