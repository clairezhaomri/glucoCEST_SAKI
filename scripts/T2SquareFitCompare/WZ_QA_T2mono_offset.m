
addpath('/Users/clairezhao/brukermatlabscripts/T1T2MapFit')

T2AcqDataPerEcho=squeeze(Data(:,:,:,9));%TR5k
% T2AcqDataPerEcho=squeeze(Data(:,:,:,4));%TR500


EffectiveTE=Data.headerHistory{1}.acqp.ACQ_echo_time;
EffectiveTE=EffectiveTE';
%%
T2AcqDataPerEcho=T2AcqDataPerEcho.*ROICtx;%
% T2AcqDataPerEcho=T2AcqDataPerEcho.*ROIMedulla;%


%%
T2Map=zeros(size(Data,1),size(Data,2));
RMSEMap=zeros(size(Data,1),size(Data,2));
SSEMap=zeros(size(Data,1),size(Data,2));
T2fitResults=zeros(3,1,size(Data,1),size(Data,2));

T2AcqData=zeros(length(EffectiveTE),1,size(Data,1),size(Data,2));


% T2AcqDataPerEchoMask=T2AcqDataPerEcho.*mask_ctx;



DataPath=pwd;
cd('/Users/clairezhao/brukermatlabscripts/T1T2MapFit/');
figure1 = figure,
axes1 = axes('Parent',figure1);
hold(axes1,'on');
    
for xi=1:size(T2AcqDataPerEcho,1)
            fprintf('%d/128\n',xi);

        for yi=1:size(T2AcqDataPerEcho,2)

                FitVector=squeeze(T2AcqDataPerEcho(xi,yi,:)); 
 %               figure,plot(TI,FitVector,'-o');
        if FitVector(1,1)==0
            continue;
%         elseif FitVector(1,1)<12
%             continue;
% % %             
%         elseif FitVector(4,1)>8
%             continue;
%         elseif FitVector(2,1)>12.5
%             continue;
%         elseif FitVector(3,1)>10
%             continue;
%         elseif FitVector(5,1)>5
%             continue;
%              elseif FitVector(6,1)>5
%             continue;
        else
                try                    
                   
                    % Set up fittype and options.
                    ft = fittype( 'a.*exp(-x./b)+c', 'independent', 'x', 'dependent', 'y' );
                    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
                    opts.Display = 'Off';
                    opts.StartPoint = [15 115 1];%[22 24 2.3] for TR500


                    % Fit model to data.
                    [fitresult, gof] = fit( EffectiveTE, FitVector, ft, opts );
                    RMSEMap(xi,yi)=gof.rmse;
                    T2Map(xi, yi) = fitresult.b;
                    T2AcqData(:,1,xi,yi)=FitVector;
%                     T2fitResults(:,1,xi,yi)=fitresult;
                    SSEMap(xi, yi)=gof.sse;
                    
%                      %Plotting
%                     plot(EffectiveTE,FitVector,'bo');hold on;
%                     plot(linspace(EffectiveTE(1),EffectiveTE(end),100),ft((fitresult.a),(fitresult.b),(fitresult.c),linspace(EffectiveTE(1),EffectiveTE(end),100)),'-','Color',[0.600000023841858 0.800000011920929 1]'); 

                    
                    plot(EffectiveTE,log10(FitVector),'bo');hold on;
                    plot(linspace(EffectiveTE(1),EffectiveTE(end),100),log10(ft((fitresult.a),(fitresult.b),(fitresult.c),linspace(EffectiveTE(1),EffectiveTE(end),100))),'-','Color',[0.600000023841858 0.800000011920929 1]'); 

                    
                catch Me

                    T2Map(xi, yi) = 0;
                    T2fitResults(:,1,xi,yi)=0;
                    T2AcqData(:,1,xi,yi,1)=0;
                    RMSEMap(xi, yi)=0;
                    SSEMap(xi, yi)=0;
                end
        end


       end
    
end
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',30);
xlim(axes1,[14 154]);
% ylim(axes1,[0 25]);%[0 120] for S,[1.5 5]for squared data,[0.5 3]for monoexponential
% ylim(axes1,[0 4]);%, [0 2]for S,[0 4] for log

ylabel('log(SI)');
% ylabel('signal intensity(a.u.)');
xlabel('Effective TE');



% raw data
FitResult.mean=mean2(nonzeros(T2Map));
FitResult.std=std2(nonzeros(T2Map));
FitResult.RMSE=mean2(nonzeros(RMSEMap));
FitResult.SSE=mean2(log(nonzeros(SSEMap)));


cd(DataPath);

save('FitResult.mat','-struct','FitResult');

