for k = 1:handles.CEST.CEST.reco.RECO_size(1)
    for l = 1:handles.CEST.CEST.reco.RECO_size(2)
        if squeeze(handles.CEST.CEST_S0.dataset(1,k,l)) > 1*10^4 % Filter Rauschen
            for i = 1:handles.CEST.CEST.method.Number_Offsets
                x(1,i) = handles.CEST.CEST.method.MT_Offsets_NoM0(1,i);
                y(1,i) = (handles.CEST.CEST.dataset(i,k,l).*handles.CEST.CEST.visu_pars.VisuCoreDataSlope(1,i))./(handles.CEST.CEST_S0.dataset(1,k,l).*handles.CEST.CEST_S0.visu_pars.VisuCoreDataSlope);
            end
            
            
            [xSort, xIndex] = sort(x);
            CEST = [(x(xIndex(:)))*0.0025; y(xIndex(:))];
            
            %         for j = 1:handles.CEST.CEST.method.Number_Offsets/2+0.5
            %             handles.CEST.asymm_Spectra_ROI(1,j) = handles.CEST.Ergebnis(2,j) - handles.CEST.Ergebnis(2,end-j+1);
            %             handles.CEST.asymm_Spectra_ROI(2,j) = (handles.CEST.xSort(j)*-1)*0.0025;
            %             %if handles.CEST.mittel(j) == handles.CEST.mittel(end-j+1), break, end
            %         end
            
            for j = 1:handles.CEST.CEST.method.Number_Offsets/2+0.5
                asymm_Contrast(j,1) = CEST(2,j) - CEST(2,end-j+1);
                asymm_Contrast(j,2) = (xSort(j)*-1)*0.0025;
            end
            
            sat_pos = find(CEST(1,:) == 0);
            sat = CEST(2,sat_pos);
            
            start_Contrast = 1;
            startPositions_Contrast = find(asymm_Contrast(:,2)<start_Contrast);
            startValue_Contrast = asymm_Contrast(startPositions_Contrast(1)-1,2);
            
            stop_Contrast = 2;
            stopPositions_Contrast = find(asymm_Contrast(:,2)>stop_Contrast);
            stopValue_Contrast = asymm_Contrast(stopPositions_Contrast(end)+1,2);
            
            asymm_Contrast = asymm_Contrast(asymm_Contrast(:,2)>=startValue_Contrast&asymm_Contrast(:,2)<=stopValue_Contrast,:);
            
            asymm_s = size(asymm_Contrast);
            
            if asymm_s(1) == 1
                AsymmValue = asymm_Contrast(1)
            else
                AsymmValue = abs(trapz(asymm_Contrast(:,2), asymm_Contrast(:,1))) / (stopValue_Contrast - startValue_Contrast);
            end
            
            if sat < 0.3 & AsymmValue<1
                handles.CEST.ContrastMatrix(k,l) = AsymmValue;
            else
                handles.CEST.ContrastMatrix(k,l) = 0;
            end
        else
            handles.CEST.ContrastMatrix(k,l) = 0;
        end
        
    end
end


N_ROI = input('Index des ROI');
handles.ROI.ROI_matrix(N_ROI,:,:) = handles.CEST.ContrastMatrix(:,:) .* handles.ROI.BI(:,:);

%roi = handles.ROI.BI(:,:);
%roi = double(roi);
%roi(find(roi == 0)) = 100;
%handles.ROI.ROI_calc(N_ROI,:,:) = (handles.CEST.ContrastMatrix(:,:) + roi(:,:))-1;
B = size(handles.ROI.ROI_matrix);
    for w = 1:B(1)
        calcmatrix = squeeze(handles.ROI.ROI_matrix(w,:,:));
        mittel = mean(calcmatrix(calcmatrix > 0));
        display(['Mittelwert:' num2str(mittel)]);
    end


matrix = zeros(B(2),B(3));
for a = 1:B(1)
    matrix = matrix + squeeze(handles.ROI.ROI_matrix(a,:,:));
end
h2 = figure;
imagesc(matrix);
colormap('hot');
cb = colorbar;
title(cb, 'MTR_a_s_y_m_m');

%figure;imagesc(handles.CEST.ContrastMatrix);colormap('hot');colorbar;

%figure;imagesc(squeeze(handles.ROI.ROI_matrix(N_ROI,:,:)));colormap('hot');colorbar;






