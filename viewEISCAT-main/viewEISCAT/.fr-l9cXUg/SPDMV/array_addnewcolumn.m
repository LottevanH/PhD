function [newdata, newcolvec]=array_addnewcolumn(data,colvec,loc,val)
  [nrow, ncol]=size(data);
  data1=data(:,1:loc);
  data2=data(:,loc+1:end);
  
  col1=colvec(1,1:loc);
  col2=colvec(1,loc+1:end);
  difcol=diff(col);
  if length(val)==1
    data3=ones(nrow,1)*val;
  elseif size(val,1)==nrow
    
    data3=val;

  else
    disp('Error: the new column cannot be added')
    return
  end
  ncol_val=size(val,2);
  step=difcol(loc)/ncol_vec;
  tempvec=1:ncol_val;
  col3=colvec(loc)+step*tempvec;
  
  newdata=[data1 data3 data2];
  newcolvec=[col1 col3 col2];
  
  return
end