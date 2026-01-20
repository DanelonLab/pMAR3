clear all
close all


Nsample=1; % Number of samples
Nimg=1; % Number of images in each sample (length should match the value of Nsample)


for sampleNumber = 1:1
    
    clearvars -except sampleNumber Nsample Nimg
    
    %% variables you need to set manually for each experiment
    Cy5_channel='C3'; %suffix for Cy5 channel
    dsGreen_channel='C1'; %suffix for dsGreen channel
    mCherry_channel='C2'; %suffix for mCherry channel

    Nimages=Nimg(sampleNumber); %assigns the number of images in this sample as defined above
    sampleName=['sample' num2str(sampleNumber)]; %name of the sample
    
    %folder in which the cropped liposome images will be saved
    %adjust the directory to 'yourfilepath/foldercontainingthecodefile/Recognized Liposomes/'
    experimentDate='2025_05_14'; %date of experiment is used to name the datafile
    folder = ['C:\Users\federicoramire\Desktop\Test-smeldit\Rebutal\Segment\Compare_SMELDiT\Modify_SMELDiT' sampleName];
     
    %%rest of the code
    
    if exist(folder, 'dir')
    else
      mkdir(folder)
    end
    
    StoreRecognizedLiposomes= true; %save cropped liposome images (default to true unless already present)
    CropSize=30; %radius of the cropped images (setting to 30 means 60x60 pixel image size)
    %Values that are used in saving the data from image analysis
    Metrics={'Area','Perimeter','Radius','Circularity','DsGreen Intensity','DsGreen variance','mCherry intensity','mCherry variance', 'Cy5 intensity', 'Cy5 variance'};
    varNames={'Area','Perimeter','R','Circ','dsGreen','dsGreen_var','mCherry','mCherry_var','Cy5','Cy5_var'};
    Units={'# pixels', '# pixels','\mum','NaN','a.u.', 'a.u.', 'a.u.', 'a.u.', 'a.u.', 'a.u.'};
    PlotScaleLog=[false false false false true true true true true true];
    f = waitbar(0, ['Processed ' num2str(0) ' out of ' num2str(Nimages) ' images.']);
    for i = 1:Nimages

    I_dsGreen = imread([sampleName '_'  num2str(i,'%03.f') '_' dsGreen_channel '.tif']);
    I_mCherry = imread([sampleName '_'  num2str(i,'%03.f') '_' mCherry_channel '.tif']);
    I_Cy5 = imread([sampleName '_'  num2str(i,'%03.f') '_' Cy5_channel '.tif']);

    [Area{i}, Perimeter{i}, dsGreen{i}, dsGreen_var{i}, centroids{i}, mCherry{i}, mCherry_var{i}, Cy5{i}, Cy5_var{i}] = image_analysis(I_Cy5, I_mCherry, I_dsGreen);
    ImInx{i}=ones(1,size(centroids{i},2)).*i;

    %   Make the images for each recognized liposome
    if StoreRecognizedLiposomes
        padded_dsGreen = padarray(I_dsGreen, [CropSize CropSize], 0);
        padded_Cy5 = padarray(I_Cy5, [CropSize CropSize], 0);
        padded_mCherry = padarray(I_mCherry, [CropSize CropSize], 0);
        centroids{i}=round(centroids{i});
        for jj = 1:size(centroids{i},2)

            IC_dsGreen=imcrop(padded_dsGreen,[centroids{i}(:,jj)' CropSize*2-1 CropSize*2-1]).*2^4;
            IC_Cy5=imcrop(padded_Cy5,[centroids{i}(:,jj)' CropSize*2-1 CropSize*2-1]).*2^4;
            IC_mCherry=imcrop(padded_mCherry,[centroids{i}(:,jj)' CropSize*2-1 CropSize*2-1]).*2^4;
            IC = cat(3, IC_mCherry, IC_dsGreen, zeros(2*CropSize,2*CropSize)); % RGB image only to see the dsGreen channel, mCherry, Cy5 are set to zero
            IC2 = cat(3, IC_Cy5,IC_Cy5,IC_Cy5); % Gray scale image for Cy5 (set to three same channels)
            IC=imfuse(IC,IC2,'blend');
            baseFileName=[sampleName '_00'  num2str(i) '_' sprintf( '%05d', jj ) '.tif'];
            fullFileName = fullfile(folder, baseFileName);
            imwrite(IC, fullFileName);

        end
    end
    waitbar(i/Nimages,f, ['Processed ' num2str(i) ' out of ' num2str(Nimages) ' images.'])
    end

    Area = cell2mat(Area);
    Perimeter = cell2mat(Perimeter);
    dsGreen= cell2mat(dsGreen);
    mCherry= cell2mat(mCherry);
    Cy5= cell2mat(Cy5);
    Cy5_var= cell2mat(Cy5_var);
    mCherry_var= cell2mat(mCherry_var);
    dsGreen_var= cell2mat(dsGreen_var);
    ImInx= cell2mat(ImInx);
    R = (2+sqrt(3*Area./(2*pi)))*0.25; %erosion disk 2
    Circ = Perimeter.^2./(4*pi*Area);

  save(['data_' experimentDate '_' sampleName],'Area', 'Perimeter','Circ','R','dsGreen','dsGreen_var','mCherry','mCherry_var','Cy5','Cy5_var','centroids','ImInx','Nimages','folder','sampleName')
end

%% image analysis

function [Area, Perimeter, dsGreen, dsGreen_var, Centroids, mCherry, mCherry_var, Cy5, Cy5_var] = image_analysis(I_Cy5,I_mCherry, I_dsGreen)

h = [-1 -1 -1; -1 12 -1; -1 -1 -1]./4; 
I_sharp = imfilter(I_Cy5,h);
I_filled = imfill(I_sharp);
I_inside = I_filled-I_sharp;
fgm = zeros(size(I_Cy5));
fgm(I_inside<=100)= 0;
fgm(I_inside>100)= 1;

fgm = imfill(fgm,'holes');
se = strel('disk',2);
fgm = imerode(fgm,se);
fgm2=zeros(size(fgm));

CC = bwconncomp(fgm);
pixels = CC.PixelIdxList;
props = regionprops(CC,'PixelIdxList','Perimeter','Area','Centroid');
passed = logical(ones(1,length(pixels)));

% select for kind of circular things
for j = 1:length(pixels)
    P(j) = props(j).Perimeter;
    A(j) = props(j).Area;

    C(j) = P(j)^2/(4*pi*A(j));
    if C(j) >2 & C(j) < 0.5
        fgm(CC.PixelIdxList{j}) = 0;
        passed(j)=false;
    elseif (2+sqrt(3*A(j)/(2*pi)))*0.125 < 0.85
        fgm(CC.PixelIdxList{j}) = 0; 
        passed(j)=false;
    elseif mean(I_mCherry(CC.PixelIdxList{j}))>500
         fgm(CC.PixelIdxList{j}) = 0; 
         passed(j)=false;
    elseif mean(I_Cy5(CC.PixelIdxList{j}))>600 % || mean(I_Cy5(CC.PixelIdxList{j}))<200
         fgm(CC.PixelIdxList{j}) = 0; 
         passed(j)=false;
    end
    
    
end

clear CC
clear pixels

% determine observables
CC = bwconncomp(fgm);
x_size=size(fgm,1);
y_size=size(fgm,2);
pixels = CC.PixelIdxList;
pixels2 = pixels;


props=props(passed);
Area=zeros(1,length(pixels));
dsGreen=zeros(1,length(pixels));
mCherry=zeros(1,length(pixels));
dsGreen_var=zeros(1,length(pixels));
Centroids=zeros(2,length(pixels));
Perimeter=zeros(1,length(pixels));

for j = 1:length(pixels)
    Area(j) = length(pixels{j});
    dsGreen(j) = mean(I_dsGreen(pixels{j}));
    dsGreen_var(j)= var(single(I_dsGreen(pixels{j})));
    Centroids(:,j)=props(j).Centroid;
    Perimeter(j) = props(j).Perimeter;

    %expand each element kk=6 pixels away from border
    pix=cell2mat(pixels(j));
    for kk=1:6
        isnotborder=(pix>x_size) & (pix<x_size*(y_size-1)+1 & ((rem(pix,x_size)~=0)&(rem(pix,x_size)~=1)));
        pix_notborder=pix(isnotborder);
        connect=[pix_notborder-x_size; pix_notborder-1; pix_notborder+1; pix_notborder+x_size];
        pix=[pix; connect];
        pix=unique(pix);
    end
    pixels2{j}=pix;
    pix=cell2mat(pixels2(j));
    pixborder{j}=pix(~ismember(pix,cell2mat(pixels(j))));
    mCherry(j) = mean(I_mCherry(pixborder{j}));
    Cy5(j) = mean(I_Cy5(pixborder{j}));
    mCherry_var(j) = var(single(I_mCherry(pixborder{j})));
    Cy5_var(j) = var(single(I_Cy5(pixborder{j})));
end
end