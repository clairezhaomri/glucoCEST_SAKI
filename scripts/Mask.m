%% creat binary mask by manual drawing. Wan-Ting Zhao.

ssImg=Data(:,:,1);% input slice for drawing (must be 2D single slice)
I=figure,imagesc(ssImg);
mask=roipoly;
waitfor(I);

mData=Data.* mask;% masked data


save('mask','mask');
save('mData','mData');