
for i = 1:handles.CEST.CEST.method.Number_Offsets
    handles.CEST.CEST_ROI(i,:,:) = ((squeeze(handles.CEST.CEST.dataset(i,:,:)) .* handles.CEST.CEST.visu_pars.VisuCoreDataSlope(1,i)).* handles.ROI.BI(:,:))./(squeeze(handles.CEST.CEST_S0.dataset(1,:,:)).*handles.CEST.CEST_S0.visu_pars.VisuCoreDataSlope);
    
    calcmatrix = squeeze(handles.CEST.CEST_ROI(i,:,:));
    handles.CEST.mittel(i) = mean(calcmatrix(calcmatrix > 0));
end

handles.CEST.xWerte = handles.CEST.CEST.method.MT_Offsets_NoM0;
handles.CEST.yWerte = handles.CEST.mittel;
[handles.CEST.xSort, handles.CEST.xIndex] = sort(handles.CEST.xWerte);
handles.CEST.Ergebnis = [(handles.CEST.xWerte(handles.CEST.xIndex(:)))*0.0025; handles.CEST.yWerte(handles.CEST.xIndex(:))];

figure;plot(handles.CEST.Ergebnis(1,:), handles.CEST.Ergebnis(2,:));
title('CEST-Spektrum');
xlabel('offset in ppm');
ylabel('S_s_a_t/S_0');

for j = 1:handles.CEST.CEST.method.Number_Offsets/2+0.5
    handles.CEST.asymm_Spectra_ROI(1,j) = handles.CEST.Ergebnis(2,j) - handles.CEST.Ergebnis(2,end-j+1);
    handles.CEST.asymm_Spectra_ROI(2,j) = (handles.CEST.xSort(j)*-1)*0.0025;
    %if handles.CEST.mittel(j) == handles.CEST.mittel(end-j+1), break, end
end

figure;plot(handles.CEST.asymm_Spectra_ROI(2,:),handles.CEST.asymm_Spectra_ROI(1,:));
title('Asymmetrie-Spektrum');
xlabel('offset in ppm');
ylabel('MTR_a_s_y_m');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% splined curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added by PB, october, 16th 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = handles.CEST.Ergebnis; % Z-Spektrum
b = handles.CEST.asymm_Spectra_ROI; % Asymm-Spektrum
aa = spline(a(1,:),a(2,:));
bb = spline(b(2,:),b(1,:));
z = linspace(min(a(1,:)),max(a(1,:)),100);
asymm = linspace(0,max(a(1,:)),100);
%figure;plot(a(1,:), a(2,:));
%figure;plot(z,ppval(aa,z));
%figure;plot(b(2,:), b(1,:));
%figure;plot(asymm,ppval(bb,asymm));
figure;plot(z,ppval(aa,z),'b',asymm,ppval(bb,asymm));
xlabel('offset [ppm]');
ylabel('Relative Intensity S_s_a_t/S_0 [%] / MTR_a_s_y_m [%]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Data save
% a = handles.CEST.Ergebnis; % Z-Spektrum
% b = handles.CEST.asymm_Spectra_ROI; % Asymm-Spektrum
% save ('test.mat', 'a', 'b');