%% Upsampling of voxel-based phantom for Phase-contrast X-ray imaging (PCXI).
% Ilian Häggmark, 2020

tic
% Load phantom with coarse voxelation. The code assumes a 16 pixel padding on all sides.
load('originalPhantom')

padSize = 16; 
smoothPadSize = 14; 		% Large padding removed after smoothing operation.
finalPadSize = 2*factor;    % Remaining padding removed after upsampling.

[row, col, layers] = size(originalPhantom);

% Remove padding
row = row - 32;
col = col - 32;
layers = layers - 32;

% Subvolume size
colSize = 10;
rowSize = 10;
layerSize = 670;
  
% material indecies
labels = [0, 29, 1, 2, 95, 125, 151, 200, 88];  
weights = [1, 1.1, 1.1, 1.1, 1.1, 2.2, 1.1, 1.1, 1.5];

NumberOfMaterials = length(labels);
factor = 32; % Upsampling done by this factor
phantomVoxelSize = 50e-6/factor;

regions = ceil(layers./layerSize); % Number of subvolumes in X-ray direction (z) to process
colParts = ceil(col./colSize);	   % Number of subvolumes in horizontal direction to process
rowParts = ceil(row./rowSize);     % Number of subvolumes in vertical direction to process

allProjectedSubvolumes = cell(rowParts,colParts);   % Cell array to keep all upsampled and projected subvolumes.


% Three nested loops to go through all subvolumes.

parfor C = 1:colParts
    cstart = (C-1)*colSize+1;
    cstop = C*colSize;
    if C == colParts
        cstop = col;
    end    
    for R = 1:rowParts
        rstart = (R-1)*rowSize+1;
        rstop = R*rowSize;
        if R == rowParts
            rstop = row;
        end

        % Matrix to keep projected thickness of all materials for one cell.
        % (To reduce size of cell array if many layers are used.)
        projOneCell = zeros(factor*(rstop-rstart+1),factor*(cstop-cstart+1),NumberOfMaterials);

        for L = 1:regions
            lstart = (L-1)*layerSize+1;
            lstop = L*layerSize;
            if L == regions
                lstop = layers;
            end


            fprintf('c %02d/%02d r %02d/%02d l %d/%d ',C,colParts,R,rowParts,L,regions)
            subvolumeTime = tic;

            % Matrix to which all upsampled materials are added to remove any overlap.
            allMaterials = zeros(factor*(rstop-rstart+1),factor*(cstop-cstart+1),factor*(lstop-lstart+1),'uint8');
            
            % Subsets are converted to uint16 [0,65000] to reduce computational time compared to single or double.

            for i = 2:NumberOfMaterials % Ignore air (i=1)
                if i == 2
                    subset = 65000*uint16(originalPhantom(rstart:(rstop+2*padSize),cstart:(cstop+2*padSize),lstart:(lstop+2*padSize)) > 0);
                else
                    subset = 65000*uint16(originalPhantom(rstart:(rstop+2*padSize),cstart:(cstop+2*padSize),lstart:(lstop+2*padSize)) == labels(i));
                end
                if sum(subset(:)) == 0 % No voxels with current label in subvolume
                    fprintf('Found no voxels  ')
                    continue % Go to to next material
                end
                
                % Core section: Smoothing, Upsampling, Edge preservation, and Re-combination.

                subset = imgaussfilt3(subset,2,'FilterSize',9);  % Standad dev. and filter size should typically be adjusted.
                subset = subset((smoothPadSize+1):(end-smoothPadSize),(smoothPadSize+1):(end-smoothPadSize),(smoothPadSize+1):(end-smoothPadSize));         
                subset = imresize3(subset,factor,'cubic');    
                subset = subset(finalPadSize+1:end-finalPadSize,finalPadSize+1:end-finalPadSize,finalPadSize+1:end-finalPadSize);
                subset = logical(round(single(subset)/65000*weights(i)));
                allMaterials(subset) = uint8(labels(i)); % Overwrite to avoid overlap.

                fprintf('Finished material %03d. ',labels(i))
            end
            fprintf('Separate materials ')
            for i = 1:NumberOfMaterials
                materialToProject = allMaterials == labels(i);
                materialToProject = sum(materialToProject,3)*phantomVoxelSize;
                projOneCell(:,:,i) = projOneCell(:,:,i) + materialToProject;
            end
            fprintf(' %0.1f ', toc(subvolumeTime))
            fprintf('Finished subvolume! \n')


        end
        allProjectedSubvolumes{R,C} = projOneCell;
    end
end

clear subset

% Transfer all matrices from the cell array to one matrix.
proj = double(zeros(factor*row,factor*col,NumberOfMaterials));

for C = 1:colParts
    cstart = (C-1)*colSize+1;
    cstop = C*colSize;
    if C == colParts
        cstop = col;
    end    
    for R = 1:rowParts
        rstart = (R-1)*rowSize+1;
        rstop = R*rowSize;
        if R == rowParts
            rstop = row;
        end
        rstartN = factor*(rstart-1)+1;
        rstopN = factor*rstop;
        cstartN = factor*(cstart-1)+1;
        cstopN = factor*cstop;
        proj(rstartN:rstopN,cstartN:cstopN,:) = cell2mat(allProjectedSubvolumes(R,C));
    end
end


% Save material separately as mat-files to limit file size.
for i = 1:NumberOfMaterials
    projmat = proj(:,:,i);
    path = sprintf('projectedMaterial%d',i);
    save(path,'projmat','-v7.3')
end

toc
fprintf('---> DONE! <--- \n')