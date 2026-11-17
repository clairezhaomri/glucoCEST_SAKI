
%A(1)=signal level (M0),more stable than M0 from T1Fit
%A(2)=T2 value
%A(3)=constant


function Y = T2Map_Fit(A,X)
  Y = A(1)*exp(-X./A(2))+A(3);

