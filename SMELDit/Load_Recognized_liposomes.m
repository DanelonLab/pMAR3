clear all 
close all
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor','DefaultAxesLineWidth'},{'k','k','k',1})

%% Parameters to set
experimentDate='2025_05_14'; %date of experiment is used to name the datafile

sampleNumber=1; %The number of the sample to load
sampleName=['Sample' num2str(sampleNumber)]; %Sample name

%folder that corresponds to the directory that contains the saved cropped liposome images
    folder = ['C:\Users\federicoramire\Desktop\Test-smeldit\Rebutal\Segment\Compare_SMELDiT\Modify_SMELDiT' sampleName];
    
%% Body of the code 
load(['data_' experimentDate '_' sampleName],'Area','R','Perimeter','Circ','dsGreen','dsGreen_var','mCherry', ...
    'mCherry_var','Cy5','Cy5_var','centroids','ImInx','Nimages')
Metrics={'Area','Radius','Perimeter','Circularity','DsGreen Intensity','DsGreen variance','mCherry intensity','mCherry variance', 'Cy5 intensity','Cy5 variance'};
varNames={'Area','R','Perimeter','Circ','dsGreen','dsGreen_var','mCherry','mCherry_var','Cy5','Cy5_var'};
Units={'# pixels', '\mum','# pixels','NaN', 'a.u.', 'a.u.', 'a.u.', 'a.u.', 'a.u.', 'a.u.', 'a.u.'};
PlotScaleLog=[false false false false true true true true true true];

Nmontage=64; %Number of images in a single montage
Nbins=100; %Bins in the 2d histogram of the heatmap

x_metric =3; y_metric=4; %The metrics that are displayed initially


figure(1)
plotHeatmap(x_metric,y_metric,Nbins)
pause(0.00001);
% frame_h = get(handle(gcf),'JavaFrame');
set(gca,'FontSize',20);
% set(frame_h,'Maximized',1);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 9/16, 1]);

%% UI buttons defined

c = uicontrol;
c.String = 'New ROI';
c.Callback = @selectROI;
c.Position = [5 65 90 20];

c2 = uicontrol;
c2.String = 'Save ROI';
c2.Callback = @saveROI;
c2.Position = [5 40 90 20];

c3 = uicontrol;
c3.String = 'Load ROI';
c3.Callback = @loadROI;
c3.Position = [5 15 90 20];

c4 = uicontrol;
c4.String = 'Export Data Selection';
c4.Callback = @exportData;
c4.Position = [105 40 120 20];

c5 = uicontrol;
c5.String = 'Import Data Selection';
c5.Callback = @importData;
c5.Position = [105 15 120 20];

c_x = uicontrol('Style','popupmenu');
c_x.String = Metrics;
c_x.Position = [350 15 90 20];
c_x.Callback = @selection_x;

c_y = uicontrol('Style','popupmenu');
c_y.String = Metrics;
c_y.Position = [450 15 90 20];
c_y.Callback = @selection_y;

set(c_y,'Value',y_metric)
set(c_y,'Callback',{@selection_y,c_x,c_y,x_metric,Nbins})

set(c_x,'Value',x_metric)
set(c_x,'Callback',{@selection_x,c_x,c_y,y_metric,Nbins})

%% UI functions

function selection_x(src,event,c_x,c_y,y_metric,Nbins)
    x_metric = c_x.Value;
    plotsurf=evalin('base', 'plotsurf');
    set(c_x,'Callback',{@selection_x,c_x,c_y,y_metric,Nbins})
    set(c_y,'Callback',{@selection_y,c_x,c_y,x_metric,Nbins})
    assignin('base','x_metric',x_metric)
    plotHeatmap(x_metric,y_metric,Nbins,plotsurf)
end

function selection_y(src,event,c_x,c_y,x_metric,Nbins)
    y_metric = c_y.Value;
    plotsurf=evalin('base', 'plotsurf');
    set(c_x,'Callback',{@selection_x,c_x,c_y,y_metric,Nbins})
    set(c_y,'Callback',{@selection_y,c_x,c_y,x_metric,Nbins})
    assignin('base','y_metric',y_metric)
    plotHeatmap(x_metric,y_metric,Nbins,plotsurf)
end

function loadROI(src,event)
    pgon=evalin('base', 'pgon');
    delete(pgon);
    [ROIfile,ROIpath] = uigetfile;
    load([ROIpath ROIfile],'roiMask','position','edgeX','edgeY','x_metric','y_metric')
    
    plotsurf=evalin('base', 'plotsurf');
    Nbins=evalin('base', 'Nbins');
    plotHeatmap(x_metric,y_metric,Nbins,plotsurf,edgeX,edgeY)
    
    pgon = polyshape(position);
    figure(1),hold on, pgon=plot(pgon);
    assignin('base','pgon',pgon)
    
    c_x=evalin('base', 'c_x');
    c_y=evalin('base', 'c_y');
    
    set(c_x,'Value',x_metric)
    set(c_y,'Value',y_metric)
    set(c_x,'Callback',{@selection_x,c_x,c_y,y_metric,Nbins})
    set(c_y,'Callback',{@selection_y,c_x,c_y,x_metric,Nbins})
    
    assignin('base','c_x',c_x)
    assignin('base','c_y',c_y)
    
    makeMontage(roiMask)
end

function saveROI(src,event)
    [ROIfile,ROIpath] = uiputfile('*.mat','Workspace File','ROIfile.mat');
    roiMask=evalin('base', 'roiMask');
    position=evalin('base', 'position');
    edgeX=evalin('base', 'edgeX');
    edgeY=evalin('base', 'edgeY');
    x_metric=evalin('base', 'x_metric');
    y_metric=evalin('base', 'y_metric');
    
    save([ROIpath ROIfile],'roiMask','position','edgeX','edgeY','x_metric','y_metric')
end

function exportData(src,event)
    [datafile,datapath] = uiputfile('*.mat','Workspace File','dataSelection.mat');
    Area=evalin('base', 'Area(testtot)');
    R=evalin('base', 'R(testtot)');
    Perimeter=evalin('base', 'Perimeter(testtot)');
    Circ=evalin('base', 'Circ(testtot)');
    dsGreen=evalin('base', 'dsGreen(testtot)');
    dsGreen_var=evalin('base', 'dsGreen_var(testtot)');
    mCherry=evalin('base', 'mCherry(testtot)');
    mCherry_var=evalin('base', 'mCherry_var(testtot)');
    Cy5=evalin('base', 'Cy5(testtot)');
    Cy5_var=evalin('base', 'Cy5_var(testtot)');
    centroids=evalin('base', 'centroids');
    ImInx=evalin('base', 'ImInx');
    Nimages=evalin('base', 'Nimages');
    testtot=evalin('base', 'testtot');
    inx=evalin('base', 'inx');
    folder=evalin('base', 'folder');
    sampleName=evalin('base', 'sampleName');
    start=1;
    
    save([datapath datafile],'Area', 'Perimeter','Circ','R','dsGreen','dsGreen_var','mCherry','mCherry_var','Cy5','Cy5_var','centroids','ImInx','Nimages','inx','sampleName','folder')
end

function importData(src,event)
    [datafile,datapath] = uigetfile;
    load([datapath datafile],'Area', 'Perimeter','Circ','R','dsGreen','dsGreen_var','mCherry','mCherry_var','Cy5','Cy5_var','centroids','ImInx','Nimages','inx','sampleName','folder')
    
    assignin('base', 'Area',Area);
    assignin('base', 'R',R);
    assignin('base', 'Perimeter',Perimeter);
    assignin('base', 'Circ',Perimeter);
    assignin('base', 'dsGreen',dsGreen);
    assignin('base', 'dsGreen_var',dsGreen_var);
    assignin('base', 'mCherry',mCherry);
    assignin('base', 'mCherry_var',mCherry_var);
    assignin('base', 'Cy5',Cy5);
    assignin('base', 'Cy5_var',Cy5_var);
    assignin('base', 'centroids',centroids);
    assignin('base', 'ImInx',ImInx);
    assignin('base', 'Nimages',Nimages);
    assignin('base', 'inxBase',inx);
    assignin('base', 'folder',folder);
    assignin('base', 'sampleName',sampleName);
    
    x_metric=evalin('base', 'x_metric');
    y_metric=evalin('base', 'y_metric');
    Nbins=evalin('base', 'Nbins');
    plotsurf=evalin('base', 'plotsurf');
    plotHeatmap(x_metric,y_metric,Nbins,plotsurf)
end

%% -- auxilliary functions

function selectROI(src,event)
        ise = evalin( 'base', 'exist(''pgon'',''var'') == 1' );
        if ise
            pgon = evalin('base', 'pgon');
            delete(pgon);
        end
        roi = impoly;
        position = getPosition(roi);
        delete(roi)
        
        Nbins=evalin('base', 'Nbins');
        edgeX=evalin('base', 'edgeX');
        edgeY=evalin('base', 'edgeY');
        position(position(:,1)<min(edgeX),1)=min(edgeX);
        position(position(:,2)<min(edgeY),2)=min(edgeY);
        position(position(:,1)>max(edgeX),1)=max(edgeX);
        position(position(:,2)>max(edgeY),2)=max(edgeY);

        pgon_out = polyshape(position);
        hold on, pgon_out=plot(pgon_out);
        assignin('base','pgon',pgon_out)
        assignin('base','position',position)
        [~,~,~,polyX,polyY] = histcounts2(position(:,1),position(:,2),edgeX,edgeY);
        polyX(end+1)=polyX(1);
        polyY(end+1)=polyY(1);
        roiMask = poly2mask(polyX,polyY,Nbins,Nbins);
        roiMask = roiMask';
        SE = strel('square',3);
        roiMask = imdilate(roiMask,SE);
        assignin('base','roiMask',roiMask)
        makeMontage(roiMask);
end

function plotHeatmap(x_metric,y_metric,Nbins,plotsurf,edgeX,edgeY)
    
    varNames=evalin('base', 'varNames');
    Metrics=evalin('base', 'Metrics');
    Units=evalin('base', 'Units');
    PlotScaleLog=evalin('base', 'PlotScaleLog');
    
    X = [evalin('base', varNames{x_metric}); evalin('base', varNames{y_metric})];
    assignin('base','X', X)
    
    
    if exist('edgeX','var') == 0
        if PlotScaleLog(x_metric)
            %edgeX=logspace(min(log(X(1,X(1,:)>0)/10)/log(10)),max(log(X(1,:)*10)/log(10)),Nbins+1);
            edgeX=logspace(min(log((X(1,:)+10^-6)/10)/log(10)),max(log(X(1,:)*10)/log(10)),Nbins+1);
        else
            edgeX=linspace(min(X(1,:)/1.05),max(X(1,:)*1.05),Nbins+1);
        end

        if PlotScaleLog(y_metric)
            %edgeY=logspace(min(log(X(2,X(2,:)>0)/10)/log(10)),max(log(X(2,:)*10)/log(10)),Nbins+1);
            edgeY=logspace(min(log((X(2,:)+10^-6)/10)/log(10)),max(log(X(2,:)*10)/log(10)),Nbins+1);
        else
            edgeY=linspace(min(X(2,:)/1.05),max(X(2,:)*1.05),Nbins+1);
        end
    
    end
    assignin('base','edgeX',edgeX)
    assignin('base','edgeY',edgeY)

    [N,~,~,BinX,BinY]=histcounts2(X(1,:)+10^-6,X(2,:)+10^-6,edgeX,edgeY);
    assignin('base','BinX',BinX)
    assignin('base','BinY',BinY)
    
    if exist('plotsurf','var')
        delete(plotsurf);
        cb=evalin('base', 'cb');
        delete(cb);
    end
    
    ise = evalin( 'base', 'exist(''pgon'',''var'') == 1' );
    if ise
        pgon = evalin('base', 'pgon');
        delete(pgon);
    end
    
    sub3 = subplot(4,4,[2 3 4 6 7 8 10 11 12]);
    plotsurf = surf(edgeX(1:end-1),edgeY(1:end-1),log(N'+1));
    assignin('base','plotsurf',plotsurf)
    
    if PlotScaleLog(x_metric)
        set(gca,'xscale','log')
    else
        set(gca,'xscale','linear')
    end
    
    if PlotScaleLog(y_metric)
        set(gca,'yscale','log')
    else
        set(gca,'yscale','linear')
    end
    
    view(0,-90)
    set(gca, 'Ydir', 'reverse')
    set(gca,'FontSize',20);
    shading flat
    map=cat(2,linspace(1,0.05,100)',linspace(1,0.6,100)',linspace(1,0.05,100)');
    map2=cat(2,linspace(0.05,0.3,100)',linspace(0.6,0,100)',linspace(0.05,0.3,100)');
    map=cat(1,map,map2);
    xlim([min(edgeX) max(edgeX(1:end-1))])
    ylim([min(edgeY) max(edgeY(1:end-1))])
%     xlabel([Metrics{x_metric} ' [' Units{x_metric} ']'])
%     ylabel([Metrics{y_metric} ' [' Units{y_metric} ']'])
%     title('SMELD it v1.1')
    pbaspect([1 1 1])
    
    subplot(4,4,[1 5 9]);
    histogram(X(2,:),edgeY,'FaceColor',[0.05, 0.6, 0.05]),camroll(90)
    set(gca,'FontSize',20);
    xlim([min(edgeY) max(edgeY(1:end-1))])
    xticks('');
%     pbaspect([3 1 1])
    if PlotScaleLog(y_metric)
        set(gca,'xscale','log')
    else
        set(gca,'xscale','linear')
    end
    subplot(4,4,[14 15 16]);
    histogram(X(1,:),edgeX,'FaceColor',[0.05, 0.6, 0.05]),camroll(180),set(gca, 'Ydir', 'reverse')
    set(gca,'FontSize',20);
    xlim([min(edgeX) max(edgeX(1:end-1))])
    xticks('');
%     pbaspect([3 1 1])
    if PlotScaleLog(x_metric)
        set(gca,'xscale','log')
    else
        set(gca,'xscale','linear')
    end
    subplot(4,4,[2 3 4 6 7 8 10 11 12])
    pos = get(sub3,'Position');
    cb = colorbar('Position',[pos(1)+pos(4) pos(2) 0.01 pos(3)]);
    colormap(gca,map)
    assignin('base','cb',cb)
    set(cb,'Ticks',cb.Ticks,...
        'TickLabels',ceil(exp(cb.Ticks))-1)

end

function plotButtonPushed(src,event,pgon,X,centroids,ImInx,c,folder,Nimages,sampleName,Nmontage,Nbins)
        if exist('pgon','var')
            delete(pgon);
        end
        roi = impoly;
        position = getPosition(roi);
        delete(roi)
        
        edgeX=evalin('base', 'edgeX');
        edgeY=evalin('base', 'edgeY');
        BinX=evalin('base', 'BinX');
        BinY=evalin('base', 'BinY');
        position(position(:,1)<min(edgeX),1)=min(edgeX);
        position(position(:,2)<min(edgeY),2)=min(edgeY);
        position(position(:,1)>max(edgeX),1)=max(edgeX);
        position(position(:,2)>max(edgeY),2)=max(edgeY);

        pgon_out = polyshape(position);
        hold on, pgon_out=plot(pgon_out);
        [~,~,~,polyX,polyY] = histcounts2(position(:,1),position(:,2),edgeX,edgeY);
        polyX(end+1)=polyX(1);
        polyY(end+1)=polyY(1);
        roiMask = poly2mask(polyX,polyY,Nbins,Nbins);
        roiMask = roiMask';
        for ii = 1:length(X(1,:))
            if (BinX(ii)==0)&&(BinY(ii)>0)
                testtot(ii) = roiMask(BinX(ii)+1,BinY(ii));
            elseif (BinY(ii)==0) && (BinX(ii)>0)
                testtot(ii) = roiMask(BinX(ii),BinY(ii)+1);
            elseif (BinY(ii)==0) && (BinX(ii)==0)
                testtot(ii) = roiMask(BinX(ii)+1,BinY(ii)+1);
            else
                testtot(ii) = roiMask(BinX(ii),BinY(ii));
            end
        end
        set(c,'Callback',{@plotButtonPushed,pgon_out,X,centroids,ImInx,c,folder,Nimages,sampleName,Nmontage,Nbins})
        assignin('base','c',c)

        inx=find(testtot);

        im_recognitions=[0];
        for ii=1:Nimages-1
            im_recognitions=[im_recognitions size(centroids{ii},2)];
        end    
        im_offset=cumsum(im_recognitions);
        
        Nmontage_used=Nmontage;
        if Nmontage > length(inx)
            Nmontage_used = length(inx);
        end

        filenamesmontage=cell(1,Nmontage_used);
        for ii=1:Nmontage_used
            jj=ii;
            Ix=inx(jj);
            baseFileName=[sampleName '_00' num2str(ImInx(inx(jj))) '_' sprintf( '%05d', Ix-im_offset(ImInx(inx(jj)))) '.tif'];
            fullFileName = fullfile(folder, baseFileName);
            filenamesmontage{ii}= fullFileName;
        end

        figure(5), montage(filenamesmontage)%, title(['Showing ' num2str(Nmontage_used) ' example liposomes'])
end

function makeMontage(roiMask)
        X=evalin('base', 'X');
        centroids=evalin('base', 'centroids');
        ImInx=evalin('base', 'ImInx');
        folder=evalin('base', 'folder');
        Nimages=evalin('base', 'Nimages');
        sampleName=evalin('base', 'sampleName');
        BinX=evalin('base', 'BinX');
        BinY=evalin('base', 'BinY');
        R = evalin('base', 'R'); 

        
        
        for ii = 1:length(X(1,:))
            if (BinX(ii)==0)&&(BinY(ii)>0)
                testtot(ii) = roiMask(BinX(ii)+1,BinY(ii));
            elseif (BinY(ii)==0) && (BinX(ii)>0)
                testtot(ii) = roiMask(BinX(ii),BinY(ii)+1);
            elseif (BinY(ii)==0) && (BinX(ii)==0)
                testtot(ii) = roiMask(BinX(ii)+1,BinY(ii)+1);
            else
                testtot(ii) = roiMask(BinX(ii),BinY(ii));
            end
        end
        f1=figure(1);
        t = uicontrol(f1,'Style','text',...
                'String',['Percentage of liposomes inside Gate: ' num2str(sum(testtot)/length(testtot)*100)],...
                'Position',[5 870 80 60]);

        assignin('base','testtot',testtot)
        
        ise = evalin( 'base', 'exist(''inxBase'',''var'') == 0' );
        if ise
            inxBase=linspace(1,length(testtot),length(testtot));
        else
            inxBase=evalin('base', 'inxBase');
        end
            
        inx=inxBase(testtot);
        
                assignin('base','inx',inx)
        
        
        im_recognitions=[0];
        for ii=1:Nimages-1
            im_recognitions=[im_recognitions size(centroids{ii},2)];
        end    
        im_offset=cumsum(im_recognitions);
        

        Nmontage_used = length(inx);

        filenamesmontage=cell(1,Nmontage_used);
        
        for ii=1:Nmontage_used
            jj=ii;
            Ix=inx(jj);
            baseFileName=[sampleName '_00' num2str(ImInx(inx(jj))) '_' sprintf( '%05d', Ix-im_offset(ImInx(inx(jj)))) '.tif'];
            fullFileName = fullfile(folder, baseFileName);
            filenamesmontage{ii}= fullFileName;
        end
        
        assignin('base','filenamesmontage',filenamesmontage)

        figure(2)
        montageSize = ceil(sqrt(min(Nmontage_used, 64))); % Maximum 8x8 grid
        for ii = 1:min(Nmontage_used, 64)
            subplot(montageSize, montageSize, ii)
    
            % Read and enhance image
            imageData = imread(filenamesmontage{ii});
            if size(imageData, 3) == 1
                % Enhance grayscale image
                imageData = imadjust(imageData, stretchlim(imageData, [0.01 0.99]), []);
            else
                % Enhance RGB image
                for channel = 1:3
                    imageData(:,:,channel) = imadjust(imageData(:,:,channel), stretchlim(imageData(:,:,channel), [0.01 0.99]), []);
                end
            end
    
            % Display enhanced image
            imshow(imageData, [])
        end


         % --------------------------------------------
        % New Block for Annotated Montage
        % --------------------------------------------
        % Cap to a maximum of 64 images
        Nmontage_used = length(inx);
        Nmontage_used = min(Nmontage_used, 64);
        
        % Display annotated montage
        figure(3)
        montageSize = ceil(sqrt(Nmontage_used)); % Grid size for montage
        for ii = 1:Nmontage_used
            subplot(montageSize, montageSize, ii)
            % Read image
            imageData = imread(filenamesmontage{ii});
            imshow(imageData, [])
            hold on
        
            % Enhance image: Adjust brightness, contrast, and saturation
            if size(imageData, 3) == 1
                % Grayscale image: Adjust brightness and contrast
                imageData = imadjust(imageData, stretchlim(imageData, [0.01 0.99]), []); % Contrast stretching
            else
                % RGB image: Adjust each channel
                for channel = 1:3
                    imageData(:,:,channel) = imadjust(imageData(:,:,channel), stretchlim(imageData(:,:,channel), [0.01 0.99]), []);
                end
            end
        
            % Show enhanced image
            imshow(imageData, [])
            hold on

            % Image center (fixed for 60x60 images)
            imageCenter = [30, 30]; % [x, y]
        
            % Get the radius for the current image
            Ix = inx(ii); % ID of the current object
            
            % Find the index of the current ID (Ix) in inxBase
            radiusIdx = find(inxBase == Ix, 1);
            
            if ~isempty(radiusIdx) && length(R) >= radiusIdx
                radius = R(radiusIdx)*4.2; % Get the radius from the mapped position
                
                % Plot center point
                plot(imageCenter(1), imageCenter(2), 'y+', 'MarkerSize', 2, 'LineWidth', 1);
            
                % Plot circle contour
                theta = linspace(0, 2 * pi, 100); % 100 points to form a smooth circle
                x = imageCenter(1) + radius * cos(theta);
                y = imageCenter(2) + radius * sin(theta);
                plot(x, y, 'y-', 'LineWidth', 2);
            else
                % If radius data is missing, display a warning in the subplot
                text(10, 10, 'No radius data', 'Color', 'red', 'FontSize', 40);
            end
            hold off
        end
        sgtitle('Annotated Montage with Circle Contours and Centroids');
end