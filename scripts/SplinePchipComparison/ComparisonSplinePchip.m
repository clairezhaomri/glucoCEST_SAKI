% test on spline and pchip on 
x = -3:3; 
y = [-1 -1 -1 0 1 1 1]; 
xq1 = -3:.01:3;
p = pchip(x,y,xq1);
s = spline(x,y,xq1);
P1=figure,plot(x,y,'o',xq1,p,'-',xq1,s,'-.','Linewidth',2.5)
legend({'Sample Points','pchip','spline'},'FontSize',15,'Location','SouthEast')
grid on


%% oscillatory dataset
x = 0:15;
y = besselj(1,x);
xq2 = 0:0.01:15;
p = pchip(x,y,xq2);
s = spline(x,y,xq2);
P2=figure,plot(x,y,'o',xq2,p,'-',xq2,s,'-.','Linewidth',2.5)
legend({'Sample Points','pchip','spline'},'FontSize',15,'Location','NorthEast')
grid on;

%% my data (KHTCH-09415)

x = linspace(-2000,2000,31);
y = [0.6302577106
0.6457131
0.6233167413
0.599094508
0.5780852923
0.5595550492
0.5371171905
0.5201071475
0.5119214241
0.4945804184
0.4813685457
0.4526386026
0.3966436929
0.312434406
0.1920497437
0.07684328643
0.07118711154
0.1552506749
0.2702653301
0.3634462729
0.435890984
0.4855608603
0.5141057489
0.5382800104
0.5697589695
0.5806416394
0.5926256412
0.6211656601
0.6542399559
0.6811259699
0.7063160613]
% SamplePoint=zeros(20,1);
% SamplePoint=[100;200;300;400;500;600;700;800;900;1000;2000;3000;4000;6000;8000;10000;15000;20000;25000;30000];
% MinimumListHeader=[100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 6000 8000 10000 15000 20000 25000 30000];
% MinimumList=zeros(2,20);

SamplePoint=zeros(16,1);
SamplePoint=[100;200;300;400;500;600;700;800;900;1000;2000;3000;4000;6000;8000;10000];
MinimumListHeader=[100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 6000 8000 10000];
MinimumList=zeros(2,16);



figure,
for i=1:numel(SamplePoint)
%     counter=ceil(i/4);
%     Num=mod(i,4);
%     if Num == 0
%         Num = 4;
%     end
    xq3=linspace(-2000,2000,SamplePoint(i,1));
    p = pchip(x,y,xq3);
    s = spline(x,y,xq3);
    subplot(4,4,i);
    plot(x,y,'o',xq3,p,'-',xq3,s,'-.','Linewidth',2.5)
grid on;
min_p=min(p);
MinimumList(1,i)=min_p;
min_s=min(s);
MinimumList(2,i)=min_s;
axis([-300 300 0 0.16 ]);
title(sprintf('%3.f pt',SamplePoint(i,1)));
legend({'Sample Pt',sprintf('Pchip_m_i_n=%3.5f',min_p),sprintf('Spline_m_i_n=%3.5f',min_s)},'Location','South');
end

cd('/Users/clairezhao/Documents/MATLAB/CEST/SplinePchipComparison/')
save('Result.mat','MinimumList','MinimumListHeader')

% calculate simulation Results
MinimumListInPercent=MinimumList;
MinimumListInPercent(1,:)=MinimumListInPercent(1,:)./MinimumListInPercent(1,end);
MinimumListInPercent(2,:)=MinimumListInPercent(2,:)./MinimumListInPercent(2,end);

figure,plot(log10(MinimumListHeader),MinimumListInPercent(1,:),'b-*');hold on;
plot(log10(MinimumListHeader),MinimumListInPercent(2,:),'r-*');
legend({'pchip','spline'},'FontSize',16);grid on;
xlabel('log_1_0(sample number)','FontSize',16);
ylabel('normalized signal(a.u.)','FontSize',16);


%% %% my data (KHTCH-09415) for asym spectrum (0-2 ppm)

x = linspace(-2000,2000,31);
y = [0.6302577106
0.6457131
0.6233167413
0.599094508
0.5780852923
0.5595550492
0.5371171905
0.5201071475
0.5119214241
0.4945804184
0.4813685457
0.4526386026
0.3966436929
0.312434406
0.1920497437
0.07684328643
0.07118711154
0.1552506749
0.2702653301
0.3634462729
0.435890984
0.4855608603
0.5141057489
0.5382800104
0.5697589695
0.5806416394
0.5926256412
0.6211656601
0.6542399559
0.6811259699
0.7063160613]
% SamplePoint=zeros(20,1);
% SamplePoint=[100;200;300;400;500;600;700;800;900;1000;2000;3000;4000;6000;8000;10000;15000;20000;25000;30000];
% MinimumListHeader=[100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 6000 8000 10000 15000 20000 25000 30000];
% MinimumList=zeros(2,20);

SamplePoint=zeros(16,1);
SamplePoint=[100;200;300;400;500;600;700;800;900;1000;2000;3000;4000;6000;8000;10000];
MinimumListHeader=[100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 6000 8000 10000];
IntegralList=zeros(2,16);
% Minimal list for (1:16);integral list for (2:16)


figure,
for i=1:numel(SamplePoint)

    xaxis=linspace(-2000,2000,SamplePoint(i,1));
    xq3=linspace(-2000,2000,SamplePoint(i,1));

    s = spline(x,y,xq3);
    
   min_s=min(s);
   IntegralList(1,i)=min_s;
   

%prepare data for aymmetric spectrum
  %make symmetric
  [min_value,min_idx]=min(s);
  a=xaxis(1,min_idx);
  xaxis_corrected1=xaxis-a;
  
  %generate symmetric datapoint
  synData_x_base=linspace(xaxis_corrected1(1),xaxis_corrected1(end),SamplePoint(i,1));
  synData_y_base=s;

  synData_y1_base= synData_y_base;
  FillNr=(length(synData_y_base)./2)-min_idx;
%         FillNr=abs(min_idx-(length(synData_y_base)./2));
  

  if isrow(synData_y_base)
  synData_y_base=synData_y_base';
  end
  
  if isrow(synData_y1_base)
  synData_y1_base=synData_y1_base';
  end

  %make sure half pt in the left hand side
  if FillNr> 0
  synData_y1_base(1:FillNr-1)=synData_y_base(1:FillNr-1);
  synData_y1_base(FillNr:end)=synData_y_base(1:end-FillNr+1);% now symmetric
  synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip

  
  elseif FillNr< 0
  qFillNr=abs(FillNr);       
  synData_y1_base(1:qFillNr)=0;
  synData_y1_base(1:end-qFillNr)=synData_y_base(qFillNr+1:end);
  synData_y1_base(end-qFillNr+1:end)=synData_y_base(end-qFillNr+1:end);% now symmetric
  synData_y1_base(end/2+1:end)=flipud(synData_y1_base(1:end/2));%now flip
  
  elseif FillNr == 0
  synData_y1_base(1:end/2)=synData_y_base(1:end/2);
  synData_y1_base(end/2+1:end)=flipud(synData_y_base(1:end/2));
  end
  

  
  
  %find display idx
  start_Contrast = 200; %1.2 ppm from previous
  startPositions_Contrast = find(synData_x_base<start_Contrast);
  startIdx = max(startPositions_Contrast);
%                 
  stop_Contrast = 800; % 2.4 ppm
  stopPositions_Contrast = find(synData_x_base<stop_Contrast);
  stopIdx = max(stopPositions_Contrast);

        pp1=spline(synData_x_base',synData_y_base);
        pp2=spline(xaxis', synData_y1_base);

        y1=ppval(pp1,synData_x_base);
        y2=ppval(pp2,synData_x_base);
        
        %calculate asymSpectrum
        asym_Base=abs((y2(startIdx:stopIdx)-(y1(startIdx:stopIdx))));%mirror trace-signal trace

%           asym_BaseForDisplay=zeros(length(synData_x_base),1);
%           asym_BaseForDisplay(startIdx1:stopIdx1)=asym_Base;

% %        correct display
%         asym_Base(1:startIdx1)=0;
%         asym_Base(stopIdx1:end)=0;
  
          asym_integral=trapz(asym_Base)./SamplePoint(i,1);% normalize for area per bin
          IntegralList(2,i)=asym_integral;
          
    subplot(4,4,i);
%     plot(x,y,'o',xq3,s,'-.',xaxis(startIdx:stopIdx),asym_Base(startIdx:stopIdx),'m','Linewidth',2.5);
    plot(x-a,y,'o',xq3-a,s,'-.',xaxis(startIdx:stopIdx),asym_Base,'m-*','Linewidth',2.5);

    grid on;
          

axis([200 800 0 0.15]);
title(sprintf('%3.f pt',SamplePoint(i,1)));
legend({'Sample Pt',sprintf('Spline_m_i_n=%3.4f',min_s),sprintf('AsymIntegral=%3.4f/pt',asym_integral)},'Location','North');

end

cd('/Users/clairezhao/Documents/MATLAB/CEST/SplinePchipComparison/')
save('SplinePtResult.mat','IntegralList','MinimumListHeader')

% calculate simulation Results
IntegralListInPercent=IntegralList;
IntegralListInPercent(1,:)=IntegralListInPercent(1,:)./IntegralListInPercent(1,end);
IntegralListInPercent(2,:)=IntegralListInPercent(2,:)./IntegralListInPercent(2,end);

figure,plot(log10(MinimumListHeader),IntegralListInPercent(1,:),'b-*');
% legend({sprintf('AsymIntegral=%3.5f',asym_integral)},'FontSize',16);
grid on;
xlabel('log_1_0(sample number)','FontSize',16);
ylabel('normalized signal(a.u.)','FontSize',16);