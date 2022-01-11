%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global drawopt

npnl=prod(size(drawopt.order));

if npnl==1 || drawopt.figure.mode==0
  drawopt.figure.unit='centimeters';
  drawopt.figure.position=[10 5 20 15];
  drawopt.figure.paperort='landscape';
  drawopt.figure.paperpos=[2.5 2 20 15];
else
  drawopt.figure.unit='pixels';
  drawopt.figure.position=[200 100 600 600];
  drawopt.figure.paperort='portrait';
  drawopt.figure.paperpos=[3 5 15 15];
end

drawopt.figure.axesfontsize=10;
drawopt.figure.textfontsize=12;

