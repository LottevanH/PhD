%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global drawopt

npnl=prod(size(drawopt.order));
default_pos=0;
if isfield(drawopt,'figure')
  if isfield(drawopt.figure,'position')
    if ~isfield(drawopt.figure,'unit')
      drawopt.figure.unit='centimeters';
    end
    drawopt.figure.paperort='portrait';
    drawopt.figure.paperpos=[2 2 18/drawopt.figure.position(4)*drawopt.figure.position(3) 18];
  else
  default_pos=1;
  end
elseif ~isfield(drawopt,'figure')
  default_pos=1;
end
if default_pos
  drawopt.figure.unit='centimeters';
  drawopt.figure.position=[10 5 18  15];
  drawopt.figure.paperort='portrait';
  drawopt.figure.paperpos=[2 2 18 15 ];
end

if ~isfield(drawopt.figure, 'axesfontsize')
    drawopt.figure.axesfontsize=10;
end
if ~isfield(drawopt.figure, 'textfontsize')
    drawopt.figure.textfontsize=10;
end

