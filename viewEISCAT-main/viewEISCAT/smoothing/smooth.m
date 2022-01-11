function s=smooth(varargin)
  a=varargin{1};
  w=varargin{2};
  type=1;
  edge=1;
  s=fastsmooth(a,w,type,edge);
  s=s';
end