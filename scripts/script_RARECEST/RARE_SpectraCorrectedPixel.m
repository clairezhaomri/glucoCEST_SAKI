test = 1+1;
%handles.CEST.ContrastMatrix_corrected(:,:)

roi_cs = handles.ROI.BI(:,:);
roi_cs = double(roi_cs);
[poszeile,posspalte] = find(roi_cs == 1);
B = size(poszeile);


%a CEST-Spectrum will be plotted for each pixel
%to avoid too many windows, a query is implemented
breaker = input(['Do you want to plot a spectrum for each pixel? (' num2str(B(1,1)) ' spectra pairs in total [y/n]'],'s');
% ---> this query might not be reasonable because the subprogramm does not
% have any output if you select 'n'; but to plot 300 spectra is not
% reasonable anyway...
for k = 1:B(1);
            for i = 1:handles.CEST.CEST.method.Number_Offsets
                x(1,i) = handles.CEST.CEST.method.MT_Offsets_NoM0(1,i);
                y(1,i) = (handles.CEST.CEST.dataset(i,poszeile(k),posspalte(k)).*handles.CEST.CEST.visu_pars.VisuCoreDataSlope(1,i))./(handles.CEST.CEST_S0.dataset(1,poszeile(k),posspalte(k)).*handles.CEST.CEST_S0.visu_pars.VisuCoreDataSlope);
            end
            
            
            [xSort, xIndex] = sort(x);
            CEST = [(x(xIndex(:)))*0.0025; y(xIndex(:))];
            
            CEST_corrected(1,:) = CEST(1,:) - handles.WASSR.ppmOffsetMatrix(poszeile(k),posspalte(k));
            CEST_corrected(2,:) = CEST(2,:);
            lastValue = CEST(1,end) - abs(handles.WASSR.ppmOffsetMatrix(poszeile(k),posspalte(k)));
            
            pp = spline(CEST_corrected(1,:)', CEST_corrected(2,:)');
            tt = linspace(-lastValue,lastValue,101);
            zSpectra_corrected = [tt; ppval(pp,tt)];
            for m = 1:100
                pos = lastValue * m/100;
                asymm_corrected(1,m) = pos;
                asymm_corrected(2,m) = ppval(pp,-pos) - ppval(pp,pos);
            end
            
            if breaker == 'y'
                %figure;plot(zSpectra_corrected(1,:),zSpectra_corrected(2,:));
                %figure;plot(asymm_corrected(1,:),asymm_corrected(2,:));
                figure;plot(zSpectra_corrected(1,:),zSpectra_corrected(2,:),'b',asymm_corrected(1,:),asymm_corrected(2,:),'b');
            end
            
end











