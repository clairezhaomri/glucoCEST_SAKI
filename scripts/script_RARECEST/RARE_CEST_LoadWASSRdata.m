%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Changes 2016 01 21 - PB
% Adaption of changes in main program code for an automated file loading
%%%%%
% input('Bitte WASSR Daten mit allen Offsets wählen [Enter drücken]')
% handles.WASSR.WASSR = read_paravision_data;
% input('Bitte WASSR Daten ohne Offset [Enter drücken]')
% handles.WASSR.WASSR_S0 = read_paravision_data;
fprintf('%s \n','Please specify the scan IDs of the WASSR measurements');
handles.numbers.WASSR_ao = input('Scan ID of the WASSR measurement with all offsets: ');
handles.numbers.WASSR_so = input('Scan ID of the WASSR measurement with S_0 single offset acquisition: ');
handles.path_WASSR_so = [handles.pathname '\' num2str(handles.numbers.WASSR_so) '\pdata\1\2dseq'];
handles.path_WASSR_ao = [handles.pathname '\' num2str(handles.numbers.WASSR_ao) '\pdata\1\2dseq'];
handles.WASSR.WASSR = read_paravision_data(handles.path_WASSR_ao);
handles.WASSR.WASSR_S0 = read_paravision_data(handles.path_WASSR_so);
% End of changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


handles.WASSR.MatrixSize = size(handles.WASSR.WASSR.dataset);
handles.WASSR.OffsetMatrix = zeros(handles.WASSR.MatrixSize(2),handles.WASSR.MatrixSize(3));
tic
for k = 1:handles.WASSR.MatrixSize(2)
    for l = 1:handles.WASSR.MatrixSize(3)

        for i = 1:handles.WASSR.WASSR.method.Number_Offsets
            x(1,i) = handles.WASSR.WASSR.method.MT_Offsets_NoM0(1,i);
            y(1,i) = (handles.WASSR.WASSR.dataset(i,k,l).*handles.WASSR.WASSR.visu_pars.VisuCoreDataSlope(1,i))./(handles.WASSR.WASSR_S0.dataset(1,k,l).*handles.WASSR.WASSR_S0.visu_pars.VisuCoreDataSlope); 
        end

        qq = spline(x(1,:),y(1,:));
        xx = -400:0.1:400;

        yMinValue = min(ppval(qq,xx));
        xPos = find(ppval(qq,xx) == yMinValue);
        xValue = xx(xPos);
        handles.WASSR.OffsetMatrix(k,l) = xValue;
    end 
end

handles.WASSR.ppmOffsetMatrix = 0.0025*handles.WASSR.OffsetMatrix;
handles.WASSR.HzOffsetMatrix = handles.WASSR.OffsetMatrix;
toc
test = 1+1;