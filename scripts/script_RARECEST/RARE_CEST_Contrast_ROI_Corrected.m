etest = 1+1;

N_ROI = input('Index des ROI');


if isfield(handles.CEST,'ContrastMatrix_corrected') == 0
    
    for k = 1:handles.CEST.CEST.reco.RECO_size(1)
        display(num2str(k));
        for l = 1:handles.CEST.CEST.reco.RECO_size(2)
            if squeeze(handles.CEST.CEST_S0.dataset(1,k,l)) > 5*10^3 % Filter Rauschen %%Standardeinstellung Annika: 10^4
                for i = 1:handles.CEST.CEST.method.Number_Offsets
                    x(1,i) = handles.CEST.CEST.method.MT_Offsets_NoM0(1,i);
                    y(1,i) = (handles.CEST.CEST.dataset(i,k,l).*handles.CEST.CEST.visu_pars.VisuCoreDataSlope(1,i))./(handles.CEST.CEST_S0.dataset(1,k,l).*handles.CEST.CEST_S0.visu_pars.VisuCoreDataSlope);
                end
                
                
                [xSort, xIndex] = sort(x);
                CEST = [(x(xIndex(:)))*0.0025; y(xIndex(:))];
                
                CEST_corrected(1,:) = CEST(1,:) - handles.WASSR.ppmOffsetMatrix(k,l);
                CEST_corrected(2,:) = CEST(2,:);
                lastValue = CEST(1,end) - abs(handles.WASSR.ppmOffsetMatrix(k,l));
                
                pp = spline(CEST_corrected(1,:)', CEST_corrected(2,:)');
                tt = linspace(-lastValue,lastValue,101);
                zSpectra_corrected = [tt; ppval(pp,tt)];
                for m = 1:1000
                    pos = lastValue * m/1000;
                    asymm_corrected(1,m) = pos;
                    asymm_corrected(2,m) = ppval(pp,-pos) - ppval(pp,pos);
                end
                
                sat_pos = find(zSpectra_corrected(1,:) == 0);
                sat = zSpectra_corrected(2,sat_pos);
                
                start_Contrast = 1; % start der Auswertung bei 1 ppm
                startPositions_Contrast = find(asymm_corrected(1,:)<start_Contrast);
                startValue_Contrast = asymm_corrected(1,startPositions_Contrast(end)+1);
                
                stop_Contrast = 2.4; % stop der Auswertung bei 2 ppm (laut Paper 1.2-2.4 ppm...)
                stopPositions_Contrast = find(asymm_corrected(1,:)<stop_Contrast);
                stopValue_Contrast = asymm_corrected(1,stopPositions_Contrast(end));
                
                if stop_Contrast - stopValue_Contrast > 0.1
                    handles.CEST.ContrastMatrix_corrected(k,l) = 0;
                else
                    
                    asymm_corrected = asymm_corrected(:,asymm_corrected(1,:)>=startValue_Contrast&asymm_corrected(1,:)<=stopValue_Contrast);
                    
                    AsymmValue = abs(trapz(asymm_corrected(1,:), asymm_corrected(2,:))) / (stopValue_Contrast - startValue_Contrast);
                    
                    if sat < 0.3 & AsymmValue<1 % Abfrage nach Sättigungsschwelle 0.3 = 70%
                        handles.CEST.ContrastMatrix_corrected(k,l) = AsymmValue;
                    else
                        handles.CEST.ContrastMatrix_corrected(k,l) = 0;
                    end
                end
            else
                handles.CEST.ContrastMatrix_corrected(k,l) = 0;
            end
        end
    end
    
    
    
    handles.ROI.ROI_matrix_cs(N_ROI,:,:) = handles.CEST.ContrastMatrix_corrected(:,:) .* handles.ROI.BI(:,:);
    
    %roi_cs = handles.ROI.BI(:,:);
    %roi_cs = double(roi_cs);
    %roi_cs(find(roi_cs == 0)) = 100;
    %handles.ROI.ROI_calc_cs(N_ROI,:,:) = (handles.CEST.ContrastMatrix_corrected(:,:) + roi_cs(:,:))-1;
    B = size(handles.ROI.ROI_matrix_cs);
    for w = 1:B(1)
        calcmatrix_cs = squeeze(handles.ROI.ROI_matrix_cs(w,:,:));
        mittel = mean(calcmatrix_cs(calcmatrix_cs > 0));
        display(['Mittelwert:' num2str(mittel)]);
    end
    
    
    
    matrix = zeros(B(2),B(3));
    for a = 1:B(1)
        matrix = matrix + squeeze(handles.ROI.ROI_matrix_cs(a,:,:));
    end
    h2 = figure;
    imagesc(matrix);
    colormap('hot');
    cb = colorbar;
    title(cb, 'MTR_a_s_y_m_m');
    
    
else
    handles.ROI.ROI_matrix_cs(N_ROI,:,:) = handles.CEST.ContrastMatrix_corrected(:,:) .* handles.ROI.BI(:,:);
    
    %roi_cs = handles.ROI.BI(:,:);
    %roi_cs = double(roi_cs);
    %roi_cs(find(roi_cs == 0)) = 100;
    %handles.ROI.ROI_calc_cs(N_ROI,:,:) = (handles.CEST.ContrastMatrix_corrected(:,:) + roi_cs(:,:))-1;
    B = size(handles.ROI.ROI_matrix_cs);
    for w = 1:B(1)
        calcmatrix_cs = squeeze(handles.ROI.ROI_matrix_cs(w,:,:));
        mittel = mean(calcmatrix_cs(calcmatrix_cs > 0));
        display(['Mittelwert:' num2str(mittel)]);
    end
    
    matrix = zeros(B(2),B(3));
    for a = 1:B(1)
        matrix = matrix + squeeze(handles.ROI.ROI_matrix_cs(a,:,:));
    end
    h2 = figure;
    imagesc(matrix);
    colormap('hot');
    cb = colorbar;
    title(cb, 'MTR_a_s_y_m_m');
end




