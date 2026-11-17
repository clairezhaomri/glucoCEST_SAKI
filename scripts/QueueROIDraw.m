function [ROI] = QueueROIDraw(path)
%Queue ROI Drawing as a new pop out figure for ROI
    if isfile('BaseMap+-5ppm:spline.mat')==1
        load('BaseMap+-5ppm:spline.mat');
    else
        sprintf('%s does not has data!', path)
    end
    
    AnatomyDataPath=fullfile(path,filesep,'pdata',filesep,'1');
    AnatomyData=MRIDataLoader.LoadFrom2DSeq(AnatomyDataPath); 
    AnatomySlice=squeeze(AnatomyData(:,:,end));
    
    if isfile('CESTCorrectData.mat')==1
        load('CESTCorrectData.mat');
        AnatomySlice=squeeze(CESTCorrectData(:,:,end));
    end
    
    if isfolder('coregistration')
    cd('coregistration')
    load('coregCEST.mat')
    AnatomySlice(:,:,:)=squeeze(coregCEST(:,:,end));
    cd ..
    end
    
    I=figure,imagesc(AnatomySlice);colormap gray;
%     load('BaseMap+-5ppm:spline.mat', 'Base_SplineIntegralResults');
%     I=figure,imagesc(Base_SplineIntegralResults,[0,60]);colormap jet;

    ROI1=roipoly;
    ROI2=roipoly;
    ROI=ROI1>0|ROI2>0;
    waitfor(I);
end

