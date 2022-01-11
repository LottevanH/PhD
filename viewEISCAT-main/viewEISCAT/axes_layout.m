%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global drawopt datasetinfo
npnl=prod(size(drawopt.order));
if ~isfield(drawopt, 'axes_height_scale')
  drawopt.axes_height_scale = ones(size(drawopt.order));
end

lsp={'-','-','-','-','-',':',':',':'};

spheiscale=2/3;
sppnl=[];
for i=1:npnl
  if ~isfinite(drawopt.plottype(i))
    para_ind=drawopt.order{i};
    if length(para_ind(:))==2 && para_ind(1)<0 && para_ind(2)<0
      sppnl=[sppnl i];
      drawopt.plottype(i)=-1;
    elseif size(para_ind,2)==2 && para_ind(1)>0 && para_ind(2)>0
      drawopt.plottype(i)=1;
    elseif ~isempty(dataset{para_ind(1)}.drawopt.plottype)
      drawopt.plottype(i)=dataset{para_ind(1,1)}.drawopt.plottype;
    else
      drawopt.plottype(i)=dataset{para_ind(1,1)}.ndim;
    end
  end  
end
nsppnl=length(sppnl);

if npnl==1
  ax=0.14; ay=0.1; ahei=0.78; awid=0.68;
else

  ax=zeros(npnl,1);ay=zeros(npnl,1);ahei=zeros(npnl,1);awid=zeros(npnl,1);
  
  tothei=0.8;
  spheiscale=2/3;
  x0=0.165;y0=0.9; wid0=0.65;
  %ahei0=tothei/(npnl-nsppnl+spheiscale*nsppnl);
  ahei0 = tothei/sum(drawopt.axes_height_scale);
  for i=1:npnl
%     if drawopt.plottype(i)==-1
%       ahei(i)=spheiscale*ahei0;
%     else
%       ahei(i)=ahei0;
%     end
    ahei(i) = drawopt.axes_height_scale(i) * ahei0;
    awid(i)=wid0;
    ax(i)=x0;
    ay(i)=y0-sum(ahei(1:i));
  end

end

% drawopt.aheiscale=[1.5 1 1 1 1 1 1 1 1];
% aheibase=tothei/sum(drawopt.aheiscale);
% ahei=drawopt.aheiscale*aheibase;
% for i=1:npnl
%    awid(i)=wid0;
%    ax(i)=x0;
%    ay(i)=y0-sum(ahei(1:i));
% end
 
drawopt.axes.ax=ax;
drawopt.axes.ay=ay;
drawopt.axes.ahei=ahei;
drawopt.axes.awid=awid;
drawopt.axes.lsp=lsp;