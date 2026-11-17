roi_cs = handles.ROI.BI(:,:);
roi_cs = double(roi_cs);
[poszeile,posspalte] = find(roi_cs == 1);
B = size(poszeile);

%Erweiterung PB, 23.10.2015 |Teil 1/3
spacer = 0;
%Erweiterung Ende
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%Bereich ist auskommentiert - Die berechneten Daten werden an keiner
%anderen Stelle verwendet, nur an einer ganz am Ende, die bei Übernahme des
%Projekts bereits auskommentiert war.
%---> kann demnach voraussichtlich gelöscht werden...
    %Das ergibt an dieser Stelle keinen Sinn! Nur wenn Werte jenseits von
    %3ppm aufgenommen werden, werden die Variablen erzeugt, die geplottet
    %werden sollen... Das ist aber grundsätzlich nicht möglich!
    %--> Die Plotbereiche müssen automatisch adjustiert werden.
%     if correctedCEST(k,1,end) >= 3 % Grenzwert bis zu welchem Wert mindestens ein Z-Spektrum konstruiert werden soll in ppm
%         calccorrectedCEST(k,1,:) = ppval(pp,linspace(-3,3,100)); % min und max wert sollten dem Grenzwert entsprechen
%         calccorrectedCEST(k,2,:) = linspace(-3,3,100); % min und max wert sollten dem Grenzwert entsprechen
%     end
%     
%     if correctedA(k,1,end) >= 3 % Grenzwert bis zu welchem Wert mindestens ein A-Spektrum konstruiert werden soll in ppm
%         pos = linspace(-3,3,100);
%         aa = spline(asymm_corrected(1,:)', asymm_corrected(2,:)');
%         calccorrectedA(k,1,:) = ppval(aa,pos(pos>0)); % min und max wert sollten dem Grenzwert entsprechen
%         calccorrectedA(k,2,:) = pos(pos>0); % min und max wert sollten dem Grenzwert entsprechen
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Erweiterung PB, 23.10.2015 | Teil 2/3
    if (correctedCEST(k,1,1) <= -2.6)
        %spectra_Zsingles_y(k-spacer,:) = correctedCEST(k,2,:);
        spectra_Zsingles_x(k-spacer,:) = correctedCEST(k,1,:);
        %spectra_Asingles_y(k-spacer,:) = correctedA(k,2,:);
        %spectra_Asingles_x(k-spacer,:) = correctedA(k,1,:);
        CEST_spectra_storage(k-spacer) = pp; %Ergibt Sinn, da die entsprechenden Punkte später wieder mit ppval berechnet werden!
    else
        spacer = spacer + 1;
    end
    %Sinn des Spacers:
    % Nur Daten im passenden Wertebereich werden in spectra_Zsingles_x
    % geschrieben.
    
    
    %Erweiterung Ende
end

%Erweiterung PB, 23.10.2015 | Teil 3/3
%Problem: Alle Spektren haben verschiedene Wertebereiche im Intervall [-3,
%3] ppm! Entsprechend sind die Abstände der Punkte zwischen den Spektren
%nicht gleich.
%Idee: Ermitteln des passenden Wertebereichs für alle Spektren und dann
%erneute Berechnung der Spektren innerhalb dieses Wertebereichs aus dem
%Spline, damit alle Spektren die gleiche Anzahl an äquidistanten Punkten
%haben.
minimum_x = max(spectra_Zsingles_x(:,1));
maximum_x = min(spectra_Zsingles_x(:,size(spectra_Zsingles_x,2)));
%Erstellung der einzelnen Spektren
for k=1:size(spectra_Zsingles_x,1)
    tt = linspace(minimum_x, maximum_x, 101);
    Spectra_Z_x(k,:) = tt; %Werte sollen alle gleich sein!
    %Spectra_Z_x(:) = tt;
    Spectra_Z_y(k,:) = ppval(CEST_spectra_storage(k),tt);
end
%Mittelwertbildung der Spektren
Spectra_Zmean = mean(Spectra_Z_y,1);
%Berechnung des Asymmetriespektrums des gemittelten Spektrums
for m = 1:floor(size(Spectra_Zmean,2)/2)
    pos = maximum_x * m/floor(size(Spectra_Zmean,2)/2);
    Spectra_A_x(m) = pos; %Werte sollen alle gleich sein!
    %Spectra_A_x(m) = pos;
    %Spectra_A_y(m) =
    %Spectra_Zmean(m)-Spectra_Zmean(size(Spectra_Zmean,2)-m); %Das müsste
    %ein spiegelverkehrtes Asymmetrie-Spektrum erzeugen
    Spectra_A_y(floor(size(Spectra_Zmean,2)/2)+1-m) = Spectra_Zmean(m)-Spectra_Zmean(size(Spectra_Zmean,2)+1-m);
end
figure;plot(Spectra_Z_x(1,:), Spectra_Zmean,'b',Spectra_A_x, Spectra_A_y,'--b');
hold on
%Berechnung der Marker für das Z-Spektrum
marker_max = round(handles.CEST.CEST.method.Final_MT_Offset/400); %Umrechnung von Hz in ppm
marker_steps = handles.CEST.CEST.method.Number_Offsets; %in vivo 21, in vitro 41
positions_marker = linspace(marker_max*(-1),marker_max,marker_steps); %variabel halten? - wahrscheinlich ja, da weniger als 21 Punkte gezeigt werden... (siehe min und max wegen Korrektur und Mittelwerte...
for k=1:marker_steps
    if (positions_marker(k) < minimum_x)
        threshold_min = k+1;
    end
    if (exist('threshold_max') ==0)
        if (positions_marker(k) > maximum_x)
            threshold_max = k-1;
        end
    end
end
positions_marker = positions_marker(threshold_min:threshold_max);
for k=1:size(spectra_Zsingles_x,1)
    Spectra_Zmarker_singles(k,:) = ppval(CEST_spectra_storage(k), positions_marker);
end
Spectra_Zmarker = mean(Spectra_Zmarker_singles,1);
clear Spectra_Zmarker_singles
plot(positions_marker, Spectra_Zmarker, 'db');
xlabel('offset [ppm]');
ylabel('Relative Intensity S_s_a_t/S_0 [%] / MTR_a_s_y_m [%]');
hold off
Spectra_Z_x = squeeze(Spectra_Z_x(1,:));
save ('CESTspectra.mat', 'Spectra_Z_x', 'Spectra_Zmean', 'Spectra_A_x', 'Spectra_A_y', 'positions_marker', 'Spectra_Zmarker');
%Erweiterung Ende


%plotcorrectedCEST = mean(calccorrectedCEST,1);
%plotcorrectedA = mean(calccorrectedA,1);

%figure;plot(squeeze(plotcorrectedCEST(1,2,:)),squeeze(plotcorrectedCEST(1,1,:)));
%figure;plot(squeeze(plotcorrectedA(1,2,:)),squeeze(plotcorrectedA(1,1,:)));
