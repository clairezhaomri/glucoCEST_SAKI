%preparation for flirt 2D
%Wan-Ting Zhao, 2019-07-03, UKJ MedPhy. Medizinische Mikrobiologie.
%because MGE share different geometry and resolution with WASSR/CEST
%data,we need to make some preparation to get the image ready for
%registration.
%flirt assume same FOV and matrix size of two sets of images,rescale didnt work for
%our case (probably because the majority of brain size dont variate much).

% flirt in 3D dataset didnt work either because our phantom didn't have structural landmark for better identification
% increase degree of freedom didn't help at all
% so we prepare the images to get ready for this.

% resize
% data passing: "template" for 2D CEST matrix,'delta21'for magnitude data of MGE
% deltaB0map

I=delta21;
J1=imresize(I,[size(template,1),size(template,2)]);
figure,imagesc(J1);

%go to fslview to check edge with cross hair