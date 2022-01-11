%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global drawopt datasetinfo
npnl=prod(size(drawopt.order));

lsp={'-','-','-','-',':',':',':',':'};

spheiscale=2/3;
sppnl=[];
for i=1:npnl
  if ~isfinite(drawopt.plottype(i))
    para_ind=drawopt.order{i};
    if length(para_ind(:))==2 && para_ind(1)>0 && para_ind(2)>0
      sppnl=[sppnl i];
      draopt.plottype(i)=-1;
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
  x0=0.15;y0=0.9; wid0=0.65;
  ahei0=tothei/(npnl-nsppnl+spheiscale*nsppnl);
  for i=1:npnl
    if drawopt.plottype(i)==-1
      ahei(i)=spheiscale*ahei0;
    else
      ahei(i)=ahei0;
    end
    awid(i)=wid0;
    ax(i)=x0;
    ay(i)=y0-sum(ahei(1:i));
  end

end
 
drawopt.axes.ax=ax;
drawopt.axes.ay=ay;
drawopt.axes.ahei=ahei;
drawopt.axes.awid=awid;
drawopt.axes.lsp=lsp;