
% A(1)=Signal level 
% A(2)=T1 Relaxation Time
% A(3)=constant

function Y = T1Map_Fit(A,X)
  Y = A(1)*(1-2*exp(-X./A(2)))+A(3);

