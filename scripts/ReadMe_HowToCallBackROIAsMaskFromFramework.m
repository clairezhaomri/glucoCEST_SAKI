%exported roi from framework, to call it back
% roi; nx1 struct

%ex: 
a=roi(1).volume.vol;

b=roi(2).volume.vol;

figure,imagesc(b); %when plot b, you'll figure out it's the mask of roi we desired.