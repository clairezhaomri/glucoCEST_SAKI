test = 1+1;

roi_cs = handles.ROI.BI(:,:);
roi_cs = double(roi_cs);
[poszeile,posspalte] = find(roi_cs == 1);
B = size(poszeile);



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
    
    correctedCEST(k,:,:) = zSpectra_corrected(:,:);
    correctedA(k,:,:) = asymm_corrected(:,:);
    
    
    if correctedCEST(k,1,end) >= 3 % Grenzwert bis zu welchem Wert mindestens ein Z-Spektrum konstruiert werden soll in ppm
        calccorrectedCEST(k,1,:) = ppval(pp,linspace(-3,3,100)); % min und max wert sollten dem Grenzwert entsprechen
        calccorrectedCEST(k,2,:) = linspace(-3,3,100); % min und max wert sollten dem Grenzwert entsprechen
    end
    
    if correctedA(k,1,end) >= 3 % Grenzwert bis zu welchem Wert mindestens ein A-Spektrum konstruiert werden soll in ppm
        pos = linspace(-3,3,100);
        aa = spline(asymm_corrected(1,:)', asymm_corrected(2,:)');
        calccorrectedA(k,1,:) = ppval(aa,pos(pos>0)); % min und max wert sollten dem Grenzwert entsprechen
        calccorrectedA(k,2,:) = pos(pos>0); % min und max wert sollten dem Grenzwert entsprechen
    end
    
end



plotcorrectedCEST = mean(calccorrectedCEST,1);
plotcorrectedA = mean(calccorrectedA,1);

figure;plot(squeeze(plotcorrectedCEST(1,2,:)),squeeze(plotcorrectedCEST(1,1,:)));
figure;plot(squeeze(plotcorrectedA(1,2,:)),squeeze(plotcorrectedA(1,1,:)));