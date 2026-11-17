function WZ_QA_T2(DataPath)

Data=MRIDataLoader.LoadFrom2DSeq([DataPath filesep 'pdata' filesep '1' filesep '2dseq']);

T2AcqDataPerEcho=squeeze(Data(:,:,:,10));%selected echo

EffectiveTE=Data.headerHistory{1}.acqp.ACQ_echo_time;
EffectiveTE=EffectiveTE';


T2Map=zeros(size(Data,1),size(Data,2));
MSEMap=zeros(size(Data,1),size(Data,2));

T2fitResults=zeros(3,1,size(Data,1),size(Data,2));

T2AcqData=zeros(length(EffectiveTE),1,size(Data,1),size(Data,2));


tic

% cd('/Users/clairezhao/brukermatlabscripts/T1T2MapFit/');

    
for xi=1:size(T2AcqDataPerEcho,1)
            fprintf('%d/128\n',xi);

        for yi=1:size(T2AcqDataPerEcho,2)

                FitVector=squeeze(T2AcqDataPerEcho(xi,yi,:)); 
 %               figure,plot(TI,FitVector,'-o');
                try
                    A0=[65,29,9];% [65 29,9]for TR500 and TR3000 rat %[86,25,5] for mice
                    fun=@(A,X)(A(1)*exp(-X./A(2))+A(3));
                    [Params,MSE] = nlinfit(EffectiveTE, FitVector,fun,A0);
                    T2fitResults(:,1,xi,yi)=[Params(1) Params(2) Params(3)];
                    T2Map(xi, yi) = Params(2);
                    T2AcqData(:,1,xi,yi)=FitVector;
                    MSEMap(xi,yi)=mean(MSE.^2);
%                     GlobalM0Map(xi, yi, sl) = Params(1);
                     % Debugging Plots
%                    figure,plot(linspace(EffectiveTE(1),EffectiveTE(end),100),T2Fit(Params,linspace(EffectiveTE(1),EffectiveTE(end),100)),'-'); 

                catch Me
%                     EffectiveTE';
%                     FitVector';
%                    plot(TI,FitVector)
                    T2Map(xi, yi) = 0;
                    T2fitResults(:,1,xi,yi)=0;
                    T2AcqData(:,1,xi,yi,1)=0;
                    MSEMap(xi, yi)=0;
                end


       end
    
end


toc

cd(DataPath);
save('T2AcqDataAndFit.mat','T2Map', 'T2fitResults','T2AcqData','EffectiveTE','MSEMap');
end