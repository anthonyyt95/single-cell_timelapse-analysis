%% 1. Separates 340 and 380 files into separate image planes for Cellpose

keepVars = {};
keepVars = [who', keepVars];
load('Input.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imageDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230904_Ca-2a2+2a3-\Images';
posTerm = input;
suffix1 = '_380';        % 380 signal
suffix2 = '_340';        % 340 signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirList = dir(imageDir);
mkdir('stackImages')
for i = [1:length(posTerm(:,1))]
    position = posTerm{i,1};
    disp(['Processing: ',position])


    % Searches for 340, 380, GFP, and Amt files
    filesFound = 0;
    for j = [1:length(dirList)]
        filename = dirList(j).name;
        loc = strfind(filename, position);
        
        if not(isempty(loc))
            is340 = strfind(filename, suffix1);
            is380 = strfind(filename, suffix2);

            if not(isempty(is340))
                filename = [imageDir,'\',filename];
                imInfo = imfinfo(filename);
                for k = [1:length(imInfo)]
                    f340Im(:,:,k) = imread(filename,k);
                end
                filesFound = filesFound + 1;
            elseif not(isempty(is380))
                filename = [imageDir,'\',filename];
                imInfo = imfinfo(filename);
                for k = [1:length(imInfo)]
                    f380Im(:,:,k) = imread(filename,k);
                end
                filesFound = filesFound + 1;
            end
        end
    end
    
    if filesFound < 2
        disp("Some image files were not found")
        continue
    end

    % Combines 340 and 380 images into individual files for Cellpose
    imPlane = 1;
    parent = ['stackImages\',position];
    for j = [1:length(f380Im(1,1,:))]
        planeIm = f380Im(:,:,j);
        planeIm(:,:,2) = f340Im(:,:,j);
        planeIm = sum(planeIm,3);
        planeIm = medfilt2(planeIm,[10,10]);
        planeIm = rescale(planeIm,0,1);
        %imshow(planeIm,[]);


        % Writes image
        fileOutName = [parent,'_imPlane',char(string(imPlane)),'_cp.tif'];
        imwrite(planeIm, fileOutName);
        imPlane = imPlane + 1;
    end
    a = 1;
end


keepVarList = '';
for i = [1:length(keepVars)]
    keepVarList = [keepVarList, ' ',keepVars{i}];
    keepVarList = strtrim(keepVarList);
end
eval(['clearvars -except ',keepVarList,';']);
disp('fin')

%% 2. Combines separated Cellpose masks into image stacks

load('Input.mat')

keepVars = {};
keepVars = [who', keepVars];
%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cellMaskDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230904_Ca-2a2+2a3-\Cellpose\cellMasks';
imageDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230904_Ca-2a2+2a3-\Images';
data = input;

%%%%%%%%%%%%%%%%%%%% Code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = [1:length(data(:,1))]
    i
    fileKey = data{i,1};
    fileKey2 = ['Fura2 - ', fileKey];
    
    % Searches for combined stack image (in 'imageDir')
    dirList = dir(imageDir);
    for j = [1:length(dirList)]
        filename = dirList(j).name;
        loc = strfind(filename,fileKey2);

        if not(isempty(loc))
            fileIm = [imageDir,'\',filename];

            stackNum = imfinfo(fileIm);
            stackNum = length(stackNum)*2;
            break
        end
    end

    % Acquires each plane of the associated cell mask and combines them in
    % order
    dirList = dir(cellMaskDir);
    fileNameOut = [fileKey,'_cp_masks.tif'];
    for j = [1:stackNum]
        fileNameOI = [fileKey,'_imPlane',char(string(j)),'_cp_cp_masks.tif'];
        
        for k = [1:length(dirList)]
            filename = dirList(k).name;
            if isequal(filename,fileNameOI)
                filename = [cellMaskDir,'\',filename];
                tmp = imread(filename);

                if j == 1
                    imwrite(tmp,fileNameOut);
                else
                    imwrite(tmp,fileNameOut,'WriteMode','append');
                end
            end
        end
    end
end
disp('fin')

keepVarList = '';
for i = [1:length(keepVars)]
    keepVarList = [keepVarList, ' ',keepVars{i}];
    keepVarList = strtrim(keepVarList);
end
eval(['clearvars -except ',keepVarList,';']);


%% 3. Performs timelapse analysis

addpath('C:\Program Files\MATLAB\R2022b\toolbox\personalScripts');

load('Input.mat');

keepVars = {'input','singleOutput', 'header'};
keepVars = [who', keepVars];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cpDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230904_Ca-2a2+2a3-\Cellpose\cellMasks';
imDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230904_Ca-2a2+2a3-\Images';
sampPrefixes = input;

manualBackgroundCorrection = 0;
outName = '20230911_singleOutput';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 1 -  Normalizes objects belonging to the same cell across cell mask stacks
dirList = dir;

if not(isfile('cpObjectNums.mat'))
    cpDirList = dir(cpDir);
    cpObjNums = sampPrefixes(:,1);
    for i = [1:length(sampPrefixes(:,1))]
        sampPrefix = sampPrefixes{i,1};

        % Searches for associated file
        foundFile = 0;
        for j = [1:length(cpDirList)]
            filename = cpDirList(j).name;
            loc1 = strfind(filename,'.tif');
            loc2 = strfind(filename,[sampPrefix,'_cp_masks']);
            if not(isempty(loc1)) & not(isempty(loc2))
                cpFilename = filename;
                foundFile = 1;
                break
            end
        end

        if foundFile == 1
            
            printText = ['Processing file:',...
                cpFilename];
            disp(printText)
    
            % Load file
            cpFilename = [cpDir,'\',cpFilename];
            cpInfo = imfinfo(cpFilename);
            for j = [1:length(cpInfo)]
                if j == 1
                    cpImage = imread(cpFilename,j);
                else
                    cpImage(:,:,j) = imread(cpFilename,j);
                end
            end
            
            % For each object in image plane #1: determines if subsequent
            % planes have an overlapping object with generally similar
            % perimeter, circularity, and area
            uniqueObjs = unique(cpImage(:,:,1));
            uniqueObjs(1) = [];
            uniqueObjsOut = {};
            cellStep = 1;
            exclude = [];
            tmp = [];
            for j = [1:length(uniqueObjs(:,1))]
                j
                % Steps through each image plane to determine if there exists
                % an object that: (1) shares similar centroid (2) overlaps (3)
                % shares similar perimeter, and (4) shares similar area with
                % objMask 1
                excludeCurrentObj = 0;
                for k = [1:length(cpImage(1,1,:))]
                    if k == 1
                        objNum1 = uniqueObjs(j);
                        objMask1 = cpImage(:,:,1);
                        objMask1(not(objMask1==objNum1)) = 0;
                        objMask1(objMask1==objNum1) = 1;
                        objMask1 = logical(objMask1);
                        
                        objMaskInfo1 = regionprops(objMask1,'Perimeter','Centroid','Area');
                        objMaskInfo1 = objMaskInfo1(1);
                        objNums = objNum1;
                        continue
                    elseif k >3
                        objNum1 = objNum2;
                        objMask1 = logical(objMask2);
                        objMaskInfo1 = objMaskInfo2;
                    end
    
                    objMask2 = cpImage(:,:,k);
                    objNum2 = objMask2(objMask1);  
                    objNum2(find(objNum2==0)) = [];
    
                    % Determines whether any object is overlapping
                    if isempty(objNum2)
                        excludeCurrentObj = 1;
                        break
                    end
    
                    objNum2 = mode(objNum2);
                    objMask2(not(objMask2==objNum2)) = 0;
                    objMask2(objMask2==objNum2) = 1;
                    objMask2 = logical(objMask2);
                    objMaskInfo2 = regionprops(objMask2, 'Perimeter','Centroid','Area');
                    objMaskInfo2 = objMaskInfo2(1);
    
                    % Determines if objMask2 overlaps with at least 75% of
                    % objectMask1
                    objOverlapPercent = objMask1 + objMask2;
                    numerator = length(find(objOverlapPercent==2));
                    denom = length(find(objOverlapPercent==1)) + numerator;
                    objOverlapPercent = (numerator/denom)*100;
                    if objOverlapPercent < 20
                        excludeCurrentObj = 1;
                    end
    
                    % Determines whether centroids area close between the 2
                    % objects
                    centDist = pdist2(objMaskInfo1.Centroid,...
                        objMaskInfo2.Centroid,...
                        'euclidean');
                    if centDist > 20
                        excludeCurrentObj = 1;
                    end
    
                    % Determines whether the perimiters area similar
                    if abs(objMaskInfo1.Perimeter-objMaskInfo2.Perimeter) > 30
                        excludeCurrentObj = 1;
                    end
    
                    % Determines whether the areas are similar
                    if abs(objMaskInfo1.Area-objMaskInfo2.Area) > 400
                        excludeCurrentObj = 1;
                    end
    
                    if excludeCurrentObj == 1
                        break
                    else
                        objNums = [objNums, objNum2];
                    end
                end
    
                if excludeCurrentObj == 0
                    uniqueObjsOut{cellStep,1} = objNums;
                    cellStep = cellStep + 1;
                end
            end
    
            cpObjNums{i,2} = uniqueObjsOut;
        end
    end
    save('cpObjectNums.mat','cpObjNums');
else
    load('cpObjectNums.mat');
end


%%%%%% 2 - Performs time-lapse analysis
imDirList = dir(imDir);
cpDirList = dir(cpDir);
singleOutput = {};
backCoord = {};
for i = [1:length(sampPrefixes(:,1))]
    %%%%% 1 - Acquires original, object mask, and outline images for each sample
    sampPrefix = sampPrefixes{i,1};
    disp(['Processing sample ', char(string(i)), '/', char(string(length(sampPrefixes(:,1)))), ...
        ': ', sampPrefix]);
 
    % Searches for .tif files (340, 380, and RFP files)
    for j = [1:length(imDirList)]
        imFile = imDirList(j).name;
        isTif = strfind(imFile,'.tif');
        is340 = strfind(imFile,'_340');
        is380 = strfind(imFile,'_380');
        isAmt = strfind(imFile,'_Zn - FRET (CFP-YFP)');
        isGFP = strfind(imFile,'_Zn - YFP')
        isSamp = strfind(imFile,sampPrefix);
        if not(isempty(isTif)) & not(isempty(isSamp))
            imPath = [imDir, '/', imFile];
            
            info = imfinfo(imPath);
            imPlane = 1;
            
            if not(isempty(is340))
                for k = [1:length(info)]
                    f340Im(:,:,k) = imread(imPath,k);
                end
            elseif not(isempty(is380))
                for k = [1:length(info)]
                    f380Im(:,:,k) = imread(imPath,k);
                end
            elseif not(isempty(isAmt))
                amtIm = imread(imPath);
            elseif not(isempty(isGFP))
                gfpIm = imread(imPath);
            end
          
        end
    end
    
   % Searches for object file & mask 
    maskNotFound = 1;
    outlineNotFound = 1;
    for j = [1:length(cpDirList)]
        imFile = cpDirList(j).name;
        isTIF = strfind(imFile, '.tif');
        isObjMask = strfind(imFile,'_cp_masks');
        isSamp = strfind(imFile,sampPrefix);

        if not(isempty(isTIF)) & not(isempty(isSamp))
            if not(isempty(isObjMask))
                imPath = [cpDir, '/', imFile];

                objMaskImInfo = imfinfo(imPath);
                for k = [1:length(objMaskImInfo)]
                    if k == 1
                        objMaskIm = imread(imPath, k);
                    else
                        objMaskIm(:,:,k) = imread(imPath,k);
                    end
                end
                maskNotFound = 0;
            end
        end
    end
    if maskNotFound == 1
        disp('No object mask file found for sample');
        continue
    end
    
    %%%% 2 - Acquires background correction
    if manualBackgroundCorrection == 1
        mergedIm = mean(f380Im,3);
        
        imshow(mergedIm/255,[]);
        [xBackCoord,yBackCoord,button] = ginputWhite(100000);
        close
        
        coordBack = {}
        for j = [1:length(xBackCoord)]
            xval = xBackCoord(j);
            yval = yBackCoord(j);
    
            xMin = ceil(xval-3);
            xMax = ceil(xval+3);
            yMin = ceil(yval-3);
            yMax = ceil(yval+3);
            
            f340Back = [];
            f380Back = [];
            for k = [1:length(f340Im(1,1,:))]
                tmp1 = mean(mean(f340Im([yMin:yMax],[xMin:xMax],k)));
                tmp2 = mean(mean(f380Im([yMin:yMax],[xMin:xMax],k)));

                f340Back = [f340Back; tmp1];
                f380Back = [f380Back; tmp2];
            end
            
            amtBack = mean(mean(amtIm([yMin:yMax],[xMin:xMax])));
            gfpBack = mean(mean(gfpIm([yMin:yMax],[xMin:xMax])));
            

            coordBack{j,1} = xval;
            coordBack{j,2} = yval;
            coordBack{j,3} = f340Back;
            coordBack{j,4} = f380Back;
            coordBack{j,5} = amtBack;
            coordBack{j,6} = gfpBack;
        end

        backCoord{i,1} = sampPrefix;
        backCoord{i,2} = coordBack;
    else
        load('backCoord.mat');
    end

    if manualBackgroundCorrection == 1
        continue
    end
    
    %%%%% 3 - Performs ratiometic calculations (and applies background
    %%%%% coorection by the coordinate)
    
    % 3.1 - Acquires unique objects according to 'cpObjNums'
    posObjNums = cpObjNums{i,2};
    if isempty(posObjNums)
        continue
    end

    uniqueObjs = {};
    for j = [1:length(posObjNums(:,1))]
        uniqueObjs = [uniqueObjs;...
            posObjNums{j,1}(1)];
    end

    for j = [1:length(uniqueObjs)]
        objNum = uniqueObjs{j};
        objMask = objMaskIm(:,:,1);
        objMask(not(objMaskIm(:,:,1)==objNum)) = 0;
        objMask(objMaskIm(:,:,1)==objNum) = 1;
        objInfo = regionprops(objMask, 'Area','BoundingBox','Centroid');
        objPerim = bwboundaries(objMask);
        
        uniqueObjs{j,2} = objInfo.Area;
        uniqueObjs{j,3} = objInfo.BoundingBox;
        uniqueObjs{j,4} = objInfo.Centroid;
        uniqueObjs{j,5} = objPerim{1};
    end

    % 3.2 - For each object: (1) applies local background correction (2)
    % acquires timelapse ratiometric calculations (3) acquires NGFR
    % intensity
    outlineIm = {};
    for j = [1:length(uniqueObjs(:,1))]
        if rem(j,5) == 0
            disp(['--Processing object ',char(string(j)),'/', char(string(length(uniqueObjs(:,1))))])
        end

        % 3.2.1 - Acquire appropriate backgrouund values for the given
        % object
        sampBack = backCoord{i,2};
        backCoordSamp = [cell2mat(sampBack(:,1)),cell2mat(sampBack(:,2))];
        minDistInd = pdist2(uniqueObjs{j,4},backCoordSamp,'Euclidean')';
        minDistInd = find(minDistInd==min(minDistInd));
        
        f340Backs = sampBack{minDistInd,3};
        f380Backs = sampBack{minDistInd,4};
        
        % 3.2.2 - For the object, for each frame: (1) applies background
        % correction (2) caluclates ratio
        objNums = posObjNums{j,1};
        f340Int = [];
        f380Int = [];
        objRatios = [];
        objSum = [];
        for k = [1:length(f340Im(1,1,:))]
            % Acquire relevant object mask number for each plane
            objNum = objNums(k);
            objMask = objMaskIm(:,:,k);
            objMask(not(objMaskIm(:,:,k)==objNum)) = 0;
            objMask(objMaskIm(:,:,k)==objNum) = 1;

            % Acquires f340 & f380 Intensity and Ratio measurements
            f340Back = f340Backs(k);
            f380Back = f380Backs(k);

            objf340 = f340Im(:,:,k);
            objf380 = f380Im(:,:,k);

            objf340 = objf340(logical(objMask)) - f340Back;
            objf380 = objf380(logical(objMask)) - f380Back;
            objf340 = mean(objf340);
            objf380 = mean(objf380);

            ratio = objf340/objf380;

            f340Int = [f340Int; objf340];
            f380Int = [f380Int; objf380];
            objSum = [objSum; objf340 + objf380];
            objRatios = [objRatios; ratio];
        end

        % 3.2.3: Acquires Amt & GFP intensity
        objNum = objNums(1);
        objMask = objMaskIm(:,:,1);
        objMask(not(objMaskIm(:,:,1)==objNum)) = 0;
        objMask(objMaskIm(:,:,1)==objNum) = 1;

        amtBack = sampBack{minDistInd,5};
        objAmt = amtIm(logical(objMask)) - amtBack;
        amtInt = mean(objAmt);

        gfpBack = sampBack{minDistInd,6};
        objGFP = gfpIm(logical(objMask)) - gfpBack;
        gfpInt = mean(objGFP);

        % uniqueObjs ammendment #2: columns 6, 7, 8, 9, 10
        uniqueObjs{j,6} = f340Int;
        uniqueObjs{j,7} = f380Int;
        uniqueObjs{j,8} = objSum;
        uniqueObjs{j,9} = objRatios;   
        uniqueObjs{j,10} = amtInt;
        uniqueObjs{j,11} = gfpInt;
    end

    objSamps = repelem({sampPrefix},length(uniqueObjs(:,1)));
    uniqueObjs = [objSamps', uniqueObjs];
    singleOutput = [singleOutput; uniqueObjs];
end

if manualBackgroundCorrection == 1
    save('backCoord.mat','backCoord');
else
    
    header = {'SampleName',
        'ObjNum',
        'ObjArea',
        'BoundingBox',
        'ObjCenter',
        'ObjPerimeter',
        'f340Intensity',
        'f380Intensity',
        'f340-380Sum',
        'f380-340Ratio',
        'AmtIntensity',
        'GFPIntensity'};
    singleOutput = [header'; singleOutput];
    
    save([outName,'.mat'],'singleOutput');
end
disp('fin')

keepVarList = '';
for i = [1:length(keepVars)]
    keepVarList = [keepVarList, ' ',keepVars{i}];
    keepVarList = strtrim(keepVarList);
end
eval(['clearvars -except ',keepVarList,';']);

%% 2 - Create plot or histogram of various metrics (and collect threshold values)


load('20230911_singleOutput.mat')
load('Input.mat')

keepVars = {'posData', 'negData1', 'negData2', 'negData4',...
    'header','avgOutput','singleOutput','timeVal'};
keepVars = [who', keepVars];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = singleOutput;
layout = input;
createGatingFigures = 1;
outName = 'output';
valueThreshold = [0 10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
header = data(1,:);
data(1,:) = [];

%%%%% 1- Plots certain parameters of interest and allows application of gating

% 1.1: Object area gating
colOI = 'ObjArea'
loc = find(strcmp(header,colOI));
val = data(:,loc);
val = cell2mat(val);
f = figure;
histogram(val,100);
[xval,yval,button] = ginput(2);
include = find(val < max(xval) & val > min(xval));
disp(['Included ', ...
    char(string(length(include))), ...
    '/', char(string(length(data(:,1)))), ...
    ' cells']);
negData1 = data;
negData1(include,:) = [];
data = data(include,:);
if createGatingFigures==1
    hold
    xlimVal = get(gca,'XLim');
    yval = get(gca,'YLim');
    pgon = polyshape([min(xval),max(xval),max(xval),min(xval)],...
        [min(yval),min(yval),max(yval),max(yval)]);
    pgon = plot(pgon);
    pgon.FaceColor = [0.5 1 0.5];
    pgon.FaceAlpha = 0.2;
    pgon.EdgeColor = [0 0 0];
    pgon.LineWidth = 1;
    xlim(xlimVal);
    ylim(yval);
    set(gca,'TickLength',[0,0]);
    set(gca,'FontSize',25);
    set(gca,'FontName','Arial');
    set(gca,'LineWidth',2);

    exportgraphics(f, ...
        [outName,'_',colOI,'_gatingGraph.pdf'], ...
        'ContentType','vector');
end
close

% 1.2: Total Fura2 intensity gating
colOI = 'f340-380Sum'
loc = find(strcmp(header,colOI));
val = [];
for i = [1:length(data(:,1))]
    val = [val; data{i,loc}(1)];
end
val = log2(val);
f = figure;
histogram(val,100);
[xval,yval,button] = ginput(2);
include = find(val < max(xval) & val > min(xval));
disp(['Included ', ...
    char(string(length(include))), ...
    '/', char(string(length(data(:,1)))), ...
    ' cells']);
negData2 = data;
negData2(include,:) = [];
data = data(include,:);
if createGatingFigures==1
    hold
    xlimVal = get(gca,'XLim');
    yval = get(gca,'YLim');
    pgon = polyshape([min(xval),max(xval),max(xval),min(xval)],...
        [min(yval),min(yval),max(yval),max(yval)]);
    pgon = plot(pgon);
    pgon.FaceColor = [0.5 1 0.5];
    pgon.FaceAlpha = 0.2;
    pgon.EdgeColor = [0 0 0];
    pgon.LineWidth = 1;
    xlim(xlimVal);
    ylim(yval);
    set(gca,'TickLength',[0,0]);
    set(gca,'FontSize',25);
    set(gca,'FontName','Arial');
    set(gca,'LineWidth',2);

    exportgraphics(f, ...
        [outName,'_',colOI,'_gatingGraph.pdf'], ...
        'ContentType','vector');
end
close

% 1.3: Gating according to Fura2 trend
colOI = 'f340-380Sum'
loc = find(strcmp(header,colOI));
f = figure;
hold
objInd = [1:length(data(:,1))];
objCoordInd = [];
objCoord = [];
xval = [1:length(data{1,loc})];
for i = [1:length(data(:,1))]
    yval = data{i,loc};
    
    p = plot(xval,yval,'Color',[0.5 0.5 0.5]);

    tmp = repelem(objInd(i),length(yval))';
    objCoordInd = [objCoordInd; tmp];
    objCoord = [objCoord;...
        xval',yval];
end
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',25);
set(gca,'FontName','Arial');
set(gca,'LineWidth',2);

roi = drawpolygon('LineWidth',2,...
    'Color',[0.5 1 0.5],...
    'FaceAlpha',0.1);
roi = roi.Position;
objKeep = find(inpolygon(objCoord(:,1),objCoord(:,2),...
    roi(:,1),roi(:,2)));
objKeep = objCoordInd(objKeep,1);
objKeep = unique(objKeep,'stable');
negData3 = data;
negData3(objKeep,:) = [];
data = data(objKeep,:);
if createGatingFigures==1
    exportgraphics(f, ...
        [outName,'_',colOI,'_gatingGraph.pdf'], ...
        'ContentType','vector');
end
close

% 1.4 Gate according to NGFR intensity
colOIX = 'GFPIntensity'
colOIY = 'AmtIntensity'
locX = find(strcmp(header,colOIX));
locY = find(strcmp(header,colOIY));
f = figure
hold
objInd = [1:length(data(:,1))];
xval = cell2mat(data(:,locX));
yval = cell2mat(data(:,locY));
for i = [1:length(yval(:,1))]               % Performs compensation
    yval(i) = yval(i)-0*xval(i);
    xval(i) = xval(i)-0*yval(i);
end
obj = logicleTransform(100,2,4,0);          % Performs biexponential transform
yval = obj.transform(yval);
xval = obj.transform(xval);
objCoord = [xval,yval];

s = scatter(xval,yval,'.');
set(gca,'TickLength',[0,0]);
set(gca,'FontSize',25);
set(gca,'FontName','Arial');
set(gca,'LineWidth',2);
%
roi = drawpolygon('LineWidth',2,...
    'Color',[0.5 1 0.5],...
    'FaceAlpha',0.1);
roi = roi.Position;
objKeep = find(inpolygon(objCoord(:,1),objCoord(:,2),...
    roi(:,1),roi(:,2)));
negData4 = data;
negData4(objKeep,:) = [];
data = data(objKeep,:);
if createGatingFigures==1
    exportgraphics(f, ...
        [outName,'_',colOIX,'+',colOIY,'_gatingGraph.pdf'], ...
        'ContentType','vector');
end
close

%%%%% 2 - Acquires metrics on gated cells (& averages based on sample)
colOI = {'f380-340Ratio'};
colInd = find(strcmp(header, colOI));

% Organizes data and acquire means
avgOutput = [];
singleOutput = layout(:,1);
exclude = [];
for i = [1:length(layout(:,1))]
    position = layout{i,1};

    loc = find(strcmp(data(:,1),position));
    if isempty(loc)
        exclude = [exclude; i];
        continue
    end
    
    timelapseValues = [];
    subdata = data(loc,:);
    for j = [1:length(subdata(:,1))]
        tmp = subdata{j,colInd};

        % Excludes cell if Inf, -Inf, zero/negative
        anyInf = find(tmp==Inf);
        anyNegInf = find(tmp==-Inf);
        anyZero = find(tmp<=0);
        if not(isempty(anyInf)) | not(isempty(anyNegInf)) | not(isempty(anyZero))
            continue
        end
        
        aboveThresh = find(tmp > valueThreshold(2));
        belowThresh = find(tmp < valueThreshold(1));
        if not(isempty(aboveThresh)) | not(isempty(belowThresh))
            continue
        end
        
        timelapseValues = [timelapseValues, tmp];
    end
    avgValues = timelapseValues;

    if i == 1
        timeVal = [0:30:length(avgValues)*30];
    end
    
    avgOutput = [avgOutput, nanmean(avgValues,2)];
    singleOutput{i,2} = timelapseValues;
end
sampHeader = layout(:,1);
sampHeader(exclude) = [];
avgOutput = [sampHeader'; num2cell(avgOutput)];

posData = [header; data];
negData1 = [header; negData1];
negData2 = [header; negData2];
negData3 = [header; negData3];
negData4 = [header; negData4];

keepVarList = '';
for i = [1:length(keepVars)]
    keepVarList = [keepVarList, ' ',keepVars{i}];
    keepVarList = strtrim(keepVarList);
end
eval(['clearvars -except ',keepVarList,';']);

%% Calculates various metrics based on selected points

data = input;       % X-by-Y cell: col.1 - time values
                    %              row.1 - column names

timeVal = cell2mat(data([2:end],1));
cond = data(1,[2:end]);
values = cell2mat(data([2:end],[2:end]));
performNorm = 1;

% Normalization parameters if checked
if performNorm == 1
    for i = [1:length(values(1,:))]
        normVal = mean(values([1:4],i));
        values(:,i) = values(:,i)./normVal;
    end
end

% Creates figure to select points
figure
hold
for i = [1:length(cond(1,:))]
    val = values(:,i);
    p = plot(timeVal,val);
end
ylim([0,5])

% Select breakpoints
[breaks,y,button] = ginput(10000);
breaks = [breaks; timeVal(end)];

% Acqures time windows
windows = {};
for i = [1:length(breaks(:,1))-1]
    if i == 1
        loc = find(timeVal <= breaks(i));
        windows = [windows; loc];
        
        loc = find(timeVal > breaks(i) & timeVal <= breaks(i+1));
        windows = [windows; loc];
    else
        loc = find(timeVal > breaks(i) & timeVal <= breaks(i+1));
        windows = [windows; loc];
    end
end

% Acquires metrics for each window
metricsOut = [];
header = {};
for i = [1:length(windows(:,1))]
    ind = windows{i,1};
    valWin = values(ind,:);
    timeWin = timeVal(ind,:);
    
    metrics = [];
    for j = [1:length(valWin(1,:))]
        val = valWin(:,j);
        
        % Calculates average
        avg = mean(val);
        
        % Calculates peak
        peak = max(val);
        
        % Calculates slope
        coefs = polyfit(timeWin,val,1);
        slope = coefs(1);
        
        % Calculates AUC
        tmp = val-val(1);
        auc = trapz(timeWin,tmp);
        
        tmp = [avg, peak, slope, auc];
        metrics = [metrics; tmp];
    end
    
    metricsOut = [metricsOut, metrics];
    colName = {'Average','Peak','Slope','AUC'};
    %colName = repelem({['window',num2str(i)]},length(metrics(1,:)));
    header = [header, colName];
end
metricsOut = [cond',num2cell(metricsOut)];

% Organizes data matrix
uniqueSamp = unique(cond,'stable');
output = {};
for i = [1:length(uniqueSamp)]
    samp = uniqueSamp{i};
    loc = find(strcmp(cond, samp));
    tmp = metricsOut(loc,:);
    output = [output; tmp];
end
header = ['Sample',header];
output = [header; output];
  
clear auc avg button coefs colName cond data header i ind j loc metrics 
clear metricsOut p peak samp slope timeVal timeWin tmp uniqueSamp
clear val values valWin windows y


%% 5 - Validate cellmask alignment algorithm

load('cpObjectNums.mat');
load('Input.mat')

keepVars = {};
keepVars = [who', keepVars];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230824_Ca-SERCA2\Images';
cpDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230824_Ca-SERCA2\Cellpose\cellMasks';
sampPrefixes = input;

suffix1 = '_340';
suffix2 = '_380';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imDirList = dir(imDir);
cpDirList = dir(cpDir);
for i = [1:length(sampPrefixes(:,1))]
    %%%%% 1 - Acquires original, object mask, and outline images for each sample
    sampPrefix = sampPrefixes{i,1};
    disp(['Processing sample ', char(string(i)), '/', char(string(length(sampPrefixes(:,1)))), ...
        ': ', sampPrefix]);
 
    % Searches for .tif file
    for j = [1:length(imDirList)]
        imFile = imDirList(j).name;
        isTif = strfind(imFile,'.tif');
        isSamp = strfind(imFile,sampPrefix);
        isF340 = strfind(imFile,suffix1);
        isF380 = strfind(imFile,suffix2);
        if not(isempty(isTif)) & not(isempty(isSamp))
            if not(isempty(isF340))
                f340File = [imDir,'/',imFile];
            elseif not(isempty(isF380))
                f380File = [imDir,'/',imFile];
            end
        end
    end

    % Creates image
    imPlane = 1;
    info = imfinfo(f340File);
    for j = [1:length(info)]
        f340Im = imread(f340File,j);
        f340Im = medfilt2(f340Im,[10,10]);
        f340Im = rescale(f340Im,0,1);
        f380Im = imread(f380File,j);
        f380Im = medfilt2(f380Im,[10,10]);
        f380Im = rescale(f380Im,0,1);

        tmpIm = f340Im;
        tmpIm(:,:,2) = f380Im;
        tmpIm = sum(tmpIm,3);

        if j == 1
            mergedIm = tmpIm;
        else
            mergedIm(:,:,j) = tmpIm;
        end
    end
    
   % Searches for object file & mask 
    maskNotFound = 1;
    outlineNotFound = 1;
    for j = [1:length(cpDirList)]
        imFile = cpDirList(j).name;
        isTIF = strfind(imFile, '.tif');
        isObjMask = strfind(imFile,'_cp_masks');
        isSamp = strfind(imFile,sampPrefix);

        if not(isempty(isTIF)) & not(isempty(isSamp))
            if not(isempty(isObjMask))
                imPath = [cpDir, '/', imFile];

                objMaskImInfo = imfinfo(imPath);
                for k = [1:length(objMaskImInfo)]
                    if k == 1
                        objMaskIm = imread(imPath, k);
                    else
                        objMaskIm(:,:,k) = imread(imPath,k);
                    end
                end
                maskNotFound = 0;
            end
        end
    end
    if maskNotFound == 1
        disp('No object mask file found for sample');
        continue
    end

    % Acquires uniqueObjs abd creates image for each image plane
    posObjNums = cpObjNums{i,2};
    if isempty(posObjNums)
        continue
    end
    
    for j = [1:length(mergedIm(1,1,:))]
        intImage = mergedIm(:,:,j);
        maskImage = objMaskIm(:,:,j);

        for k = [1:length(posObjNums(:,1))]
            objNum = posObjNums{k,1}(j);

            maskImage(maskImage==objNum) = 10000;
        end
        maskImage(not(maskImage==10000)) = 0;
        maskImage(maskImage==-1) = 1;
        maskImage = logical(maskImage);

        intImage(not(maskImage)) = 0;

        % Saves mask image and intensity image
        if j == 1
            imwrite(intImage,...
                [sampPrefix,'_intensityImage.tif']);
        else
            imwrite(intImage,...
                [sampPrefix,'_intensityImage.tif'],'WriteMode','append');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keepVarList = '';
for i = [1:length(keepVars)]
    keepVarList = [keepVarList, ' ',keepVars{i}];
    keepVarList = strtrim(keepVarList);
end
eval(['clearvars -except ',keepVarList,';']);

%% 6 - Generates Excel containing single cell values for each condition

load('20230924_DP_output.mat')
load('Input.mat')

data = singleOutput;
sampPrefixes = input;
sampPrefixColInd = 2;

colOut = {};
valOut = [];
nanSep = nan(length(data{1,2}),1);
for i = [1:length(data(:,1))]
    samp = sampPrefixes{i,sampPrefixColInd};
    values = data{i,2};
    
    valOut = [valOut, values, nanSep];

    colTmp = repelem({samp},length(values(1,:)));
    colTmp = [colTmp, {'NaN'}];
    colOut = [colOut, colTmp];
end
output = [colOut; num2cell(valOut)];

xlswrite('singleCellOutput.xlsx', output)

%% 6. Organizes output data (single cell) for graphPad

keepVars = [who'];
%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = singleOutput;
ind = 2-1;

samps = data(:,1);
data = data(:,[2:end]);
uniqueSamps = unique(samps,'stable');
for i = [1:length(uniqueSamps)]
    samp = uniqueSamps{i};
    loc = find(strcmp(samps,samp));

    if not(isempty(samp))
        sampVal = [];
        for j = [1:length(loc)]
            val = data{loc(j),ind};

            sampVal = [sampVal; val];
        end
    end

    uniqueSamps{i,2} = sampVal;
end

% Creates table for excel
maxArrayLen = [];
for i = [1:length(uniqueSamps(:,1))]
    maxArrayLen = [maxArrayLen;
        length(uniqueSamps{i,2})];
end
maxArrayLen = max(maxArrayLen);

valOut = [];
for i = [1:length(uniqueSamps(:,1))]
    val = uniqueSamps{i,2};
    
    appendNaN = maxArrayLen - length(val);
    appendNaN = nan(appendNaN,1);
    val = [val; appendNaN];

    valOut = [valOut, val];
end

% Saves & opens excel file
xlswrite('output.xlsx',valOut,'Sheet1');
winopen('output.xlsx');

%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

keepVarList = '';
for i = [1:length(keepVars)]
    keepVarList = [keepVarList, ' ', keepVars{i}];
end
eval(['clearvars -except', keepVarList, ';']);




%% Tmp

imageDir = 'C:\Users\antho\OneDrive\Desktop\Misc\20230707_Ca-2a2-\Images\tmp';

dirList = dir(imageDir);
mkdir([imageDir,'\newfolder']);
parentNew = [imageDir,'\newfolder'];
for i = [1:length(dirList)]
    filename = dirList(i).name
    
    isTIF = strfind(filename,'.tif');
    if not(isempty(isTIF))
        fileOrigin = [imageDir,'\',filename];
    
        filename = strsplit(filename,'_');
        filename = [filename{1}, {'Iono'}, filename([2:end])];
        filename = strjoin(filename,'_');
        
        fileNew = [parentNew,'\',filename];
    
        copyfile(fileOrigin,fileNew);
    end
end















