function [co]=getlineplotcolor(m)
  co=[];
  if m==1
    co=[0 0 0];
  elseif m==2
    co=[0 76 153; 230 76 0]/255;
  elseif m==3
    co=[0 0 0; 1 0 0; 0 0 1];
  elseif m==4
    cl=othercolor('Set14',m);
    co(1,:)=cl(1,:);co(2,:)=cl(3,:);co(3,:)=cl(4,:);co(4,:)=cl(3,:);
  else
    co=othercolor('Set14',m);
  end
end
