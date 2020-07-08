%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Title:    Target Embedding Algorithm (T.E.A)
% Author:   Adam Clayden
% Date:     02/05/2020
% Note:		Currently supports 4 target sizes only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TargetEmbeddingAlgorithm()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START OF CUSTOMISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global columnName    
    imagePath = ''; % path where your images are stored
    exportPath = ''; % path you wish to export final images to
    img_info = readtable(''); % path containing a list of image names from (e.g. .xls file)
    columnName = '';  % name of the column to read image names from
    
	% image width and height
    desiredImageWidth = 800;
    desiredImageHeight = 600;

    % monitor values
    x_monitor = 40.5; % visible monitor width (cm)
    y_monitor = 30.5; % visible monitor height (cm)
    
    % distance between monitor and subject (cm)
    MonitorSubjectDistance = 90;

    % radius of the central circular exclusion mask in degrees
    radius_deg = 3;
    % size of a region buffer to place around the target to define the
    % region size
    region_buffer = 3;
    
    %thickness of the letter
    %These are example values that you can change
    horizontalWidth = [1 2 2 4];
    verticalWidth = [1 3 3 5];
    
    %actual width and height
    %These are example values that you can change
    widthVector = [7 13 19 33];
    heightVector = [9 16 23 40];
    
    %the index of each target size, with the length being the count
    letterMatrix = [1 2 3 4];
    
    valueType = 'Median'; %Min, Lower Quartile, Median, Upper Quartile, Max
    
    %Change this to find suitable locations based on a certain target size
    %Largest size is recommended and is defined by the number within
    %letterMatrix. E.g. 1, 2, 3, etc
    placeWithTypeThreshold = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CUSTOMISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global contrast_map_difference finalImageData nrScenes summedMaps
    global newMins newLQs newMedians newUQs newMaxs CircleMask
    radius_pix = ComputeCircleMask(MonitorSubjectDistance, desiredImageWidth, desiredImageHeight, x_monitor, radius_deg);
    PrepareImages(img_info, desiredImageWidth, desiredImageHeight, imagePath);
    noOfLetters = GenerateDifferenceMaps(desiredImageWidth, desiredImageHeight, region_buffer, letterMatrix, heightVector, widthVector, verticalWidth, horizontalWidth);
    selectedValue = GenerateSummedMaps(valueType, noOfLetters);
    WriteImages(letterMatrix, exportPath, noOfLetters, heightVector, widthVector, verticalWidth, horizontalWidth, selectedValue, desiredImageWidth, desiredImageHeight, radius_pix, region_buffer);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates a central circular mask to use as the exclusion zone later on
%
% Input: Distance between monitor and subject, image width/height,
%        x dimension of monitor, radius of central mask in degrees
%
% Output: radius_pix - radius of the circular exclusion zone in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function radius_pix = ComputeCircleMask(MonitorSubjectDistance, desiredImageWidth, desiredImageHeight, x_monitor, radius_deg)
    global CircleMask
    % compute degree-per-pixel value                    
    ang = atan2(1,MonitorSubjectDistance)*180/pi; % how many degrees correspond to 1 cm on the subject screen?
    r = desiredImageWidth/x_monitor;     % 1280 pixel (x axis) equal 40.5 cm => 'r' carries the ratio between pixel and cm
    DEG_PER_PIX = ang / r; 
    PIX_PER_DEG = 1/DEG_PER_PIX;
    radius_pix = radius_deg * PIX_PER_DEG;
    N = 256;
    t = (0:N)*2*pi/N;
    x = radius_pix*cos(t)+desiredImageWidth/2;
    y = radius_pix*sin(t)+desiredImageHeight/2;
    CircleMask = poly2mask(x,y,desiredImageHeight,desiredImageWidth);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in images via name from the submitted image list and then resizes
% them and converts to greyscale as needed. After images are stored in an
% array
%
% Input: image names, image width/height, the image path
%
% Output: Nothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PrepareImages(img_info, res_x, res_y, imagePath)
    global finalImageData nrScenes columnName
    imageList = img_info{:,columnName};
    nrScenes = size(img_info,1);
    
    finalImageData = cell(nrScenes,1);
    
    for i = 1:nrScenes
        imgFile = imageList{i};
        imdata = imread([imagePath imgFile], 'jpg');
        if size(imdata,2) ~= res_x || size(imdata,1) ~= res_y
            imdata = imresize(imdata,[res_y res_x]);
            fprintf(1,'\nwrong image resolution! image resized.\n');
        end

        if size(imdata,3) == 3 % if colored image
            imdata = rgb2gray(imdata);
        end
        
        finalImageData{i}.image = imdata;
        finalImageData{i}.name = imgFile;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loops through all images over all target sizes to calculate contrast 
% difference maps using RMS contrast
%
% Input: image width/height, region buffer, letter indices, letter thickness and size
%
% Output: Number of letter sizes (length of letterMatrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function noOfLetters = GenerateDifferenceMaps(desiredImageWidth, desiredImageHeight, region_buffer, letterMatrix, heightVector, widthVector, verticalWidth, horizontalWidth)
    global contrast_map_difference nrScenes finalImageData 
    noOfLetters = length(letterMatrix);
    contrast_map_difference = cell(nrScenes, noOfLetters);
    for iIndex = 1 : nrScenes
        for fontIndex = 1 : noOfLetters
            hWidth = horizontalWidth(fontIndex);
            vWidth = verticalWidth(fontIndex);
            myHeight = heightVector(fontIndex);
            myWidth = widthVector(fontIndex);
            mySize = max([myHeight myWidth]);
            if mod(mySize,2) ~= 0 % if odd number
                mySize = mySize + 1; % add 1 to make it even number
            end
            centraliseLetter = floor(hWidth/2);
            region_size = mySize + 2*region_buffer;
            colStart = 1 + mySize/2+region_buffer;
            rowStart = 1 + mySize/2+region_buffer;
            colabsEnd = desiredImageWidth - region_size/2;
            rowabsEnd = desiredImageHeight - region_size/2;
            % initialize output matrices for contrast values
            contrast_map_withLet = ones(desiredImageHeight,desiredImageWidth)*NaN;
            contrast_map_withoutLet = ones(desiredImageHeight,desiredImageWidth)*NaN;
            currentImage = finalImageData{iIndex}.image;
            currentName  = finalImageData{iIndex}.name; % get current image name
            meanCimg     = quantile(currentImage(:),0.50);
            fprintf(['processing image ',num2str(iIndex),currentName])
            % move letter pixelwise and calculate local contrast
            for rowIndex = rowStart:rowabsEnd % rowStart:desiredImageHeight % row
                if mod(rowIndex,desiredImageHeight/10) == 0
                    fprintf([num2str(rowIndex/(desiredImageHeight/10)*10) ' done.\n'])
                end
                for colIndex = colStart:colabsEnd
                    newImage = currentImage;
                    newImage(rowIndex-floor(mySize/2):rowIndex-floor(mySize/2) + mySize,colIndex-centraliseLetter:colIndex-centraliseLetter + (vWidth-1)) = 0;
                    newImage((rowIndex-floor(mySize/2)):(rowIndex-floor(mySize/2)) + (hWidth-1),colIndex-floor(myWidth/2):colIndex-floor(myWidth/2) + myWidth-1) = 0;
                    % cut out evaluation region = box around the letter
                    newSection = newImage(rowIndex-region_size/2:rowIndex+region_size/2,colIndex-region_size/2:colIndex+region_size/2);
                    % compute local contrast for box *with* letter
                    vec = reshape(newSection,size(newSection,1)*size(newSection,2),1);
                    dev = std(double(vec));
                    lcontrast    = dev/quantile(newImage(:),0.50);
                    contrast_map_withLet(rowIndex,colIndex) = lcontrast;
                    % compute local contrast for box *without* letter
                    oldSection = currentImage(rowIndex-region_size/2:rowIndex+region_size/2,colIndex-region_size/2:colIndex+region_size/2);
                    vec2 = reshape(oldSection,size(oldSection,1)*size(oldSection,2),1);
                    dev2 = std(double(vec2));
                    lcontrast2   = dev2/meanCimg;
                    contrast_map_withoutLet(rowIndex,colIndex) = lcontrast2;
                end
            end
            contrast_map_difference{iIndex, fontIndex} = contrast_map_withLet - contrast_map_withoutLet;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates the summed difference maps for each image and also calculates
% the value type depending on the user's preference (e.g. median etc)
%
% Input: value type (e.g. median, min etc), the number of letter sizes (
% output from GenerateDifferenceMaps)
%
% Output: selectedValue (numeric form of valueType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function selectedValue = GenerateSummedMaps(valueType, noOfLetters)
    global contrast_map_difference nrScenes summedMaps
    global newMins newLQs newMedians newUQs newMaxs
    
    selectedValue = -1;
    switch valueType
    case 'Min'
        selectedValue = 1;
    case 'Lower Quartile'
        selectedValue = 2;
    case 'Median'
        selectedValue = 3;
    case 'Upper Quartile'
        selectedValue = 4;
    case 'Max'
        selectedValue = 5;
    otherwise
        return;
    end
    
    newMins = cell(1,nrScenes);
    newLQs = cell(1,nrScenes);
    newMedians = cell(1,nrScenes);
    newUQs = cell(1,nrScenes);
    newMaxs = cell(1,nrScenes);
    summedMaps = cell(1,nrScenes);
    
    for sumIndex = 1:nrScenes
        currentMap = contrast_map_difference{sumIndex,1};
        for cFont = 2:noOfLetters
            currentMap = abs(currentMap + contrast_map_difference{sumIndex,cFont});
            if(cFont == noOfLetters)
                summedMaps{sumIndex} = currentMap;
            end
        end
    end
    
    switch selectedValue
    case 1
        for sumIndex = 1:nrScenes
            newMins{sumIndex} = quantile(summedMaps{sumIndex}(:),0);
        end
    case 2
        for sumIndex = 1:nrScenes
            newLQs{sumIndex} = quantile(summedMaps{sumIndex}(:),0.25);
        end
    case 3
        for sumIndex = 1:nrScenes
            newMedians{sumIndex} = quantile(summedMaps{sumIndex}(:),0.50);
        end
    case 4
        for sumIndex = 1:nrScenes
            newUQs{sumIndex} = quantile(summedMaps{sumIndex}(:),0.75);
        end
    case 5
        for sumIndex = 1:nrScenes
            newMaxs{sumIndex} = quantile(summedMaps{sumIndex}(:),1.0);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loops through all images at all target sizes, calls 'PlaceWithType' function
% to determine final location for the letter, and writes images to export path
%
% Input: letter indices, image export path, number of letter sizes, letter
% thickness and size, selected value type (numeric), x/y resolution of image,
% radius of circle in pixels, region buffer
%
% Output: nothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteImages(letterMatrix, exportPath, noOfLetters, heightVector, widthVector, verticalWidth, horizontalWidth, selectedValue, res_x, res_y, radius_pix, region_buffer)
    global finalImageData nrScenes summedMaps CircleMask
    global newMins newLQs newMedians newUQs newMaxs
    
    switch selectedValue
    case 1
        currentValue = newMins;
    case 2
        currentValue = newLQs;
    case 3
        currentValue = newMedians;
    case 4
        currentValue = newUQs;
    case 5
        currentValue = newMaxs;
    end
    
    for imageIndex = 1:nrScenes
        for fontIndex = 1:noOfLetters
            myHeight = heightVector(fontIndex);
            myWidth = widthVector(fontIndex);
            mySize = max([myHeight myWidth]);
            if mod(mySize,2) ~= 0 % if odd number
                mySize = mySize + 1; % add 1 to make it even number
            end
            switch fontIndex
                case 1
                    final1Image = finalImageData{imageIndex}.image;
                    hWidth1 = horizontalWidth(fontIndex);
                    vWidth1 = verticalWidth(fontIndex);
                    centraliseLetter1 = floor(vWidth1/2);
                    myHeight1 = myHeight;
                    myWidth1 = myWidth;
                    mySize1 = mySize;
                    region_size1 = mySize1+2*region_buffer;
                    if(letterMatrix(end) == 1)
                        finHwidth = horizontalWidth(fontIndex);
                        finVwidth = verticalWidth(fontIndex);
                        finCentraliseLetter = floor(vWidth1/2);
                        finHeight = myHeight;
                        finWidth = myWidth;
                        finSize = mySize;
                        finRegionSize = mySize1+2*region_buffer;
                    end
                case 2
                    final2Image = finalImageData{imageIndex}.image;
                    hWidth2 = horizontalWidth(fontIndex);
                    vWidth2 = verticalWidth(fontIndex);
                    centraliseLetter2 = floor(vWidth2/2);
                    myHeight2 = myHeight;
                    myWidth2 = myWidth;
                    mySize2 = mySize;
                    region_size2 = mySize2+2*region_buffer;
                    if(letterMatrix(end) == 2)
                        finHwidth = horizontalWidth(fontIndex);
                        finVwidth = verticalWidth(fontIndex);
                        finCentraliseLetter = floor(vWidth2/2);
                        finHeight = myHeight;
                        finWidth = myWidth;
                        finSize = mySize;
                        finRegionSize = mySize2+2*region_buffer;
                    end
                case 3
                    final3Image = finalImageData{imageIndex}.image;
                    hWidth3 = horizontalWidth(fontIndex);
                    vWidth3 = verticalWidth(fontIndex);
                    centraliseLetter3 = floor(vWidth3/2);
                    myHeight3 = myHeight;
                    myWidth3 = myWidth;
                    mySize3 = mySize;
                    region_size3 = mySize3+2*region_buffer;
                    if(letterMatrix(end) == 3)
                        finHwidth = horizontalWidth(fontIndex);
                        finVwidth = verticalWidth(fontIndex);
                        finCentraliseLetter = floor(vWidth3/2);
                        finHeight = myHeight;
                        finWidth = myWidth;
                        finSize = mySize;
                        finRegionSize = mySize3+2*region_buffer;
                    end
                case 4
                    final4Image = finalImageData{imageIndex}.image;
                    hWidth4 = horizontalWidth(fontIndex);
                    vWidth4 = verticalWidth(fontIndex);
                    centraliseLetter4 = floor(vWidth4/2);
                    myHeight4 = myHeight;
                    myWidth4 = myWidth;
                    mySize4 = mySize;
                    if(letterMatrix(end) == 4)
                        finHwidth = horizontalWidth(fontIndex);
                        finVwidth = verticalWidth(fontIndex);
                        finCentraliseLetter = floor(vWidth4/2);
                        finHeight = myHeight;
                        finWidth = myWidth;
                        finSize = mySize;
                        finRegionSize = mySize4+2*region_buffer;
                    end
                otherwise
                return;
            end
            if(fontIndex == length(letterMatrix))
                [fxCoor,fyCoor,~,~,~,~,~,~,~,~] = PlaceWithType(currentValue{imageIndex},finWidth,finHeight,summedMaps{imageIndex},CircleMask,res_x,res_y,finSize,radius_pix);
                imageName = finalImageData{imageIndex}.name;
                final4Image(fyCoor-floor(mySize4/2):fyCoor-floor(mySize4/2) + mySize4,fxCoor-centraliseLetter4:fxCoor-centraliseLetter4 + (vWidth4-1)) = 0;
                final4Image((fyCoor-floor(mySize4/2)):(fyCoor-floor(mySize4/2)) + (hWidth4-1),fxCoor-floor(myWidth4/2):fxCoor-floor(myWidth4/2) + myWidth4-1) = 0;
                final3Image(fyCoor-floor(mySize3/2):fyCoor-floor(mySize3/2) + mySize3,fxCoor-centraliseLetter3:fxCoor-centraliseLetter3 + (vWidth3-1)) = 0;
                final3Image((fyCoor-floor(mySize3/2)):(fyCoor-floor(mySize3/2)) + (hWidth3-1),fxCoor-floor(myWidth3/2):fxCoor-floor(myWidth3/2) + myWidth3-1) = 0;
                final2Image(fyCoor-floor(mySize2/2):fyCoor-floor(mySize2/2) + mySize2,fxCoor-centraliseLetter2:fxCoor-centraliseLetter2 + (vWidth2-1)) = 0;
                final2Image((fyCoor-floor(mySize2/2)):(fyCoor-floor(mySize2/2)) + (hWidth2-1),fxCoor-floor(myWidth2/2):fxCoor-floor(myWidth2/2) + myWidth2-1) = 0;
                final1Image(fyCoor-floor(mySize1/2):fyCoor-floor(mySize1/2) + mySize1,fxCoor-centraliseLetter1:fxCoor-centraliseLetter1 + (vWidth1-1)) = 0;
                final1Image((fyCoor-floor(mySize1/2)):(fyCoor-floor(mySize1/2)) + (hWidth1-1),fxCoor-floor(myWidth1/2):fxCoor-floor(myWidth1/2) + myWidth1-1) = 0;

                imwrite(final1Image,[exportPath imageName(1:end-4) 's.jpg']);
                imwrite(final2Image,[exportPath imageName(1:end-4) 'm.jpg']);
                imwrite(final3Image,[exportPath imageName(1:end-4) 'l.jpg']);
                imwrite(final4Image,[exportPath imageName(1:end-4) 'xl.jpg']);
                region_size4 = mySize4+2*region_buffer;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds target location closest to the value type inputted that satisfies
% exclusion criteria (not in center, not too close to the edges). There are 
% several output variables for debugging purposes, but fxCoor and fyCoor are
% the important ones as they're the final x/y positions of the target
%
%Input:
%   valueType: target value which can be the min, lower quartile, median,
%   upper quartile or the max within each difference map.
%   fontWidth, fontHeight: target object's width and height
%   csumMap: current summed difference map
%   fid: file identifier (1 if writing images, otherwise, 0)
%   choiceStr: name of the valueType eg: "median"
%   CircleMask: circular section in the centre to exclude
%   res_x, res_y: 600x800
%   mySize75: max[width, height]
%   radiusC: radius of the CircleMask
%
%Output:
%   fxCoor, fyCoor: final x,y positions to use for object insertion
%   xArray, yArray: x,y arrays of all invalid coordinates
%   xCri, yCri: all locations that passed first criterion but not second
%   xallPos,yallPos: all locations of the same value within a matrix
%   trueValue: returns either "True median" or "Close median" depending on
%   the outcome. Or, whichever valueType was chosen.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fxCoor,fyCoor,xArray,yArray,xCri,yCri,xallPos,yallPos,trueValue, rowcolArray] = PlaceWithType(valueType,fontWidth,fontHeight,csumMap,CircleMask,res_x, res_y,letterSize,radiusC)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Instantiate Variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    breakloop = 0; %Break condition for while loop
    xallPos = []; %All x positions if there are duplicate values
    yallPos = []; %All y positions if there are duplicate values
    xArray = []; %Invalid x location
    yArray = []; %Invalid y location
    xCri = []; %All x positions that passed both criterion
    yCri = []; %All y positions that passed both criterion
    fxCoor = []; %Final x coordinate chosen
    fyCoor = []; %Final y coordinate chosen
    size75 = letterSize; %Size of letter (max between width and height)
    letterRad = size75/2;
    trueValue = ''; %Blank string for visualisation to tell user if the value chosen was a true value or a value close to it
    rcidx = 0; %row col index for passed criterion
    cRow = res_y/2; %centre of circleMask
    cCol = res_x/2; %centre of circleMask

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find value type (e.g. median) in Image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [sumedY12, sumedX12] = find(csumMap == valueType);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin Loop
    % Finds value type in summed difference map or a value close to said
    % value.
    % The value (median, max etc) could be located in multiple positions within
    % a given image, and so each one is analysed before 1 being chosen at
    % random to yield our final x,y coordinates.
    % Equation that the loop represents: SUMi(SUMj(1< Xij < WH : Xij ~= C[ij]))
    % x,y coordinates of value is checked agaisnt two criteria:
    % 1. Boundaries of image to prevent letter from truncating off
    % 2. Inner circle to imitate foveal vision as letter should be placed
    %    outside of this region (either parafoveal or peripheral region)
    % If value passes both checks then we choose this value.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while(breakloop == 0)
        % If empty
        if(isempty(sumedX12))
            % Create temporary matrix of absolute values
            tmpsumlq = abs(csumMap - valueType);
            % Reshape into vector
            tmpVector = reshape(tmpsumlq,1,numel(tmpsumlq));
            % Transpose vector and place values in first column
            mat3(:,1) = tmpVector';
            % Place rows in second column
            mat3(:,2) = repmat(1:size(csumMap,1),1,size(csumMap,2))';
            tmp = [];
            % Place columns in third column
            for i=1:size(csumMap,2)
                tmp = [tmp; repmat(ones(1,size(csumMap,1))*i,1)'];
            end
            mat3(:,3) = tmp;
            % Sort matrix rows by the first column to place coordinates in
            % order
            mat3 = sortrows(mat3,1);
            % Find all unique values in the first column
            tmpvals = mat3(:,1);
            tmpvals_unique = unique(tmpvals);
            % Open for loop 1 - number of unique values
            for nextmin = 1:length(tmpvals_unique)
                % Find current unique value in matrix and store how many times
                % it occurs
                idx = find(mat3(:,1)==tmpvals_unique(nextmin));
                for posIndex = 1:length(idx)
                    xallPos(end+1) = mat3(idx(posIndex),3);
                    yallPos(end+1) = mat3(idx(posIndex),2);
                end
                % Open for loop to iterate through these multiple values (if
                % any)
                for rpInd = 1:length(idx)
                    % Store current value and coordinates
                    mypick = idx(rpInd);
                    row = mat3(mypick,2); % X is 2 or 3
                    col = mat3(mypick,3);
                    % Check against criterion (boundary + centre)
                    if(row > (1 + fontHeight) && row < (res_y - fontHeight) && col > (1 + fontWidth) && col < (res_x - fontWidth))
                        d = sqrt((cRow - row).^2 + (cCol - col).^2);
                        if(d>radiusC+letterRad)
                            % add row and column values into array
                            rcidx = rcidx + 1;
                            rowcolArray(rcidx,2) = col;
                            rowcolArray(rcidx,1) = row;
                            if(rpInd == length(idx))
                                breakloop = 1;
                            end
                        else
                            xCri(end+1) = col;
                            yCri(end+1) = row;
                        end
                    else
                        % Else = invalid location
                        xArray(end+1) = col;
                        yArray(end+1) = row;
                    end
                    % If all locations don't pass, loop to the next uniqe value
                    % and store how many times it occurs
                end            
                if(breakloop == 1)
                    % Select one at random and make this our final
                    % position
                    ranPerm = randperm(numel(rowcolArray(1:rcidx,1)));
                    fxCoor = rowcolArray(ranPerm(1),2);
                    fyCoor = rowcolArray(ranPerm(1),1);
                    break;
                end
            end
        else
            % If not empty, check each value against criteria (if there are
            % multiple)
            rcidx = 0;
            row = sumedY12;
            col = sumedX12;
            for tmidx = 1:numel(row)
                if(row(tmidx) > (1 + fontHeight) && row(tmidx) < (res_y - fontHeight) && col(tmidx) > (1 + fontWidth) && col(tmidx) < (res_x - fontWidth))
                    d = sqrt((cRow - row).^2 + (cCol - col).^2);
                    if(d>radiusC+letterRad)
                        % Add row and column values into array
                        rcidx = rcidx + 1;
                        rowcolArray(rcidx,2) = col(tmidx);
                        rowcolArray(rcidx,1) = row(tmidx);
                        % Select one at random and make this our final position
                        ranPerm = randperm(numel(rowcolArray(:,1)));
                        fxCoor = rowcolArray(ranPerm(1),2);
                        fyCoor = rowcolArray(ranPerm(1),1);
                        if(tmidx == numel(row))
                            breakloop = 1;
                        end
                    else
                        xCri(end+1) = col(tmidx);
                        yCri(end+1) = row(tmidx);
                    end
                else
                    xArray(end+1) = col(tmidx);
                    yArray(end+1) = row(tmidx);
                end
            end
            % If all fail, loop back up to top and find next closest
            if(breakloop == 0)
                sumedX12 = [];
            end
        end
        if(breakloop == 1)
            % Select one at random and make this our final position
            ranPerm = randperm(numel(rowcolArray(1:rcidx,1)));
            fxCoor = rowcolArray(ranPerm(1),2);
            fyCoor = rowcolArray(ranPerm(1),1);
            break;
        end
    end
end