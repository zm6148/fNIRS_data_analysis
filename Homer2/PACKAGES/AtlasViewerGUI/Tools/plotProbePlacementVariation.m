function plotProbePlacementVariation()

% Selection 1: Inter-subject variability analysis:
%   Plots multiple subject probe geometries on a common atlas volume
%   Plots mean optode locations and their standard deviations
% Selection 2: Probe fabrication error analysis:
%   Plots original source-detector geometries versus digitized probe(s)
%   Plots original optode locations and the digitized mean error distribution


% Find atlasViewer.mat files located one subdirectory deep.
folderContents = dir;
numberOfObjects = length(folderContents);
countSubjects = 0;
for i = 3:1:numberOfObjects
    if (folderContents(i).isdir ~= false)
        tempString = strcat(folderContents(i).name, '/atlasViewer.mat');
        if (exist(tempString, 'file') ~= false)
            countSubjects = countSubjects + 1;
            subFiles{countSubjects} = tempString;
        end
    end
end

% Make sure we have at least 2 subjects
if countSubjects<2
    menu('Warning: need at least 2 subjects with atlasViewer.mat file to calculate probe placement variation.', 'OK');
    return;
end

wd = cd();
pathnm = uigetdir(wd,'Select Group Root Folder containing the Subject Folders');
if pathnm == 0
    return
end
cd(pathnm)

userSelection = menu('Analysis Type','Inter-subject variability','Probe fabrication error');

% Define the distance threshold that isolates short separation channels
ssThresh = 15; % millimeters

if (userSelection == 2)
        
    if ~exist('atlasViewer_SDdesign.mat','file')
        menu('You need the file atlasViewer_SDdesign.mat that contains the original probe design registered to the surface copied into the root of the group folder.','Okay');
        return;
    end
    
    % Setup and Plot Atlas
    atlas = load('atlasViewer_SDdesign.mat', '-mat'); % An atlasViewer.mat 
    % file that contains the original probe design registered to the surface
    % copied into the root of the group folder

    % Head Surface
    f = atlas.headsurf.mesh.faces;
    v = atlas.headsurf.mesh.vertices;
    c = v(:,3); % Creating a vector of equal length to v, to be used as color values
    c(:) = 100; % 100, using clim [100, 200], produces a blue color
    fN = trisurfnorm(v, f); % This function is provided in a separate file

    % Pial Surface
    fp = atlas.pialsurf.mesh.faces;
    vp = atlas.pialsurf.mesh.vertices;
    cp = vp(:,3); % Creating a vector of equal length to vp, to be used as color values
    cp(:) = 150; % 150, using clim [100, 200], produces a green color

    figure('name','SD Layout (Black) and Digitized Probe Positions');
    clf
    hold on
    set(gca,'clim',[100 200])
    h = trisurf(f, v(:,1), v(:,2), v(:,3), c(:)); % Plot head surface
    set(h, 'linestyle', 'none', 'facealpha', 0.25) % Remove lines and set opacity to 25%
    h2 = trisurf(fp, vp(:,1), vp(:,2), vp(:,3), cp(:)); % Plot pial surface
    set(h2, 'linestyle', 'none', 'facealpha', 1) % Remove lines and set opacity to 100%
    lighting phong
    light('position',[-1 0 1])
    light('position',[0 -1 0])
    light('position',[1 0 -1])
    axis('off')
    axis image
    view(0, 0)

    % Original Circle for Patches
    temp = 1:1:21;
    radius = 2.5;
    x = radius * cos(2 * pi * temp(:) / size(temp, 2));
    y = radius * sin(2 * pi * temp(:) / size(temp, 2));
    z(1:size(temp, 2)) = 0;
    circle1 = [x(:), y(:), z(:)];
    circ1N = [0 0 1]; % Normal to the plane of circle 1
    circle2 = circle1; % Create a second circle variable to be modified each loop of s (below)

    % Define the colors to be used for the optodes
    colors1 = colormap;
    for i = 1:1:countSubjects
        patchColor(i,:) = colors1(1+round(63*i/countSubjects),:);
    end

    % Translate each subject's optode positions onto a common atlas and plot
    rotEig(1:3,1:3,1:countSubjects) = 0;
    optN = 100; % Temporary value in case the following loop doesn't run
    for s = 1:1:countSubjects
        % Load Subject Data, Transform Probe, Plot Circles Normal to Surface
        sub = load(subFiles{s}, '-mat');
        p = sub.probe.optpos_reg; % Original registered optode positions
        p(:, 4) = 1;
        T = sub.headvol.T_2mc; % Translation matrix to individual subject space
        Ti(1:4, 1:4) = inv(T(1:4, 1:4)); % Inverse produces translation to common atlas space

        % Search for optodes that are closer together than the ssThresh and exclude them
        useList(1:size(p, 1)) = 1;
        for i = 1:1:size(p, 1)
            for j = 1:1:size(p, 1)
                if (sqrt((p(i, 1) - p(j, 1))^2 + (p(i, 2) - p(j, 2))^2 ...
                    + (p(i, 3) - p(j, 3))^2) < ssThresh)
                    % The following statements assume sources are listed 
                    % before detectors and therefore the optode to exclude 
                    % has the higher number
                    if (i > j)
                        useList(i) = 0;
                    end
                    if (j > i)
                        useList(j) = 0;
                    end
                end
            end
        end
        optN = sum(useList(:));
        p2(1:optN,1:4) = 0;
        j = 1;
        for i=1:1:size(p, 1)
            if (useList(i) ~= 0) % Only entries in the useList are carried to subsequent steps
                p2(j, 1:4) = transpose(Ti * transpose(p(i, 1:4))); % New optode positions in common space
                j = j + 1;
            end
        end

        p2n(1:optN, 1:3) = 0; % Create an array to store the normal vector from the common atlas surface at each optode position
        for i=1:1:optN
            vmin = 1e99; vidx = 0;
            for j = 1:1:size(v, 1)
                foo = sqrt((p2(i,1)-v(j,1))^2 + (p2(i,2)-v(j,2))^2 + (p2(i,3)-v(j,3))^2); % distance from optode to vertex
                if (foo < vmin) % Search for vertex with the minimum distance from the optode
                    vmin = foo;
                    vidx = j;
                end
            end
            [lst, foo] = find(f==vidx); % Using the nearest vertex, find the matching faces
            for j = 1:1:size(lst, 1)
                p2n(i, :) = p2n(i, :) + fN(lst(j), :); % Sum the normal vectors for all matching faces
            end
            p2n(i, :) = p2n(i, :) / size(lst, 1); % Divide by the number of faces to find the average normal vector
            p2T = vrrotvec2mat_copy(vrrotvec_copy(circ1N, p2n(i, :))); % Create a rotation matrix using the original circle normal and the average normal vector for the optode location
            for j = 1:1:size(circle1, 1)
                circle2(j, :) = transpose(p2T * transpose(circle1(j, :))) + p2(i, 1:3); % Rotate and translate the original circle shape into the optode plane and position
            end
            patch(circle2(:, 1), circle2(:, 2), circle2(:, 3), patchColor(s,:)); % Creates a patch using the optode circle on the head volume
        end

        [foo, bar] = eig(sub.headvol.T_2mc(1:3,1:3)); % Calculate the eigenvalues of the rotation matrix from subject space to common atlas space
        rotEig(1:3,1:3,s) = abs(bar); % Store the absolute value of the rotation matrix eigenvalues
        p2all(:,1:4,s) = p2; % Store the subject's translated optode locations for group plot
    end

    % Small Circle for Probe  Patches
    temp = 1:1:21;
    radius = 1;
    x = radius * cos(2 * pi * temp(:) / size(temp, 2));
    y = radius * sin(2 * pi * temp(:) / size(temp, 2));
    z(1:size(temp, 2)) = 0;
    circle1 = [x(:), y(:), z(:)];
    circ1N = [0 0 1]; % Normal to the plane of circle 1
    circle2 = circle1; % Create a second circle variable to be modified each loop of s (below)
    
    % Plot the original probe design (black)
    p2n(1:optN, 1:3) = 0; % Create an array to store the normal vector from the common atlas surface at each optode position
    p2 = atlas.probe.optpos_reg; % Original probe design
    for i=1:1:optN
        vmin = 1e99; vidx = 0;
        for j = 1:1:size(v, 1)
            foo = sqrt((p2(i,1)-v(j,1))^2 + (p2(i,2)-v(j,2))^2 + (p2(i,3)-v(j,3))^2); % distance from optode to vertex
            if (foo < vmin) % Search for vertex with the minimum distance from the optode
                vmin = foo;
                vidx = j;
            end
        end
        [lst, foo] = find(f==vidx); % Using the nearest vertex, find the matching faces
        for j = 1:1:size(lst, 1)
            p2n(i, :) = p2n(i, :) + fN(lst(j), :); % Sum the normal vectors for all matching faces
        end
        p2n(i, :) = p2n(i, :) / size(lst, 1); % Divide by the number of faces to find the average normal vector
        p2T = vrrotvec2mat_copy(vrrotvec_copy(circ1N, p2n(i, :))); % Create a rotation matrix using the original circle normal and the average normal vector for the optode location
        for j = 1:1:size(circle1, 1)
            circle2(j, :) = transpose(p2T * transpose(circle1(j, :))) + p2(i, 1:3); % Rotate and translate the original circle shape into the optode plane and position
        end
        patch(circle2(:, 1), circle2(:, 2), circle2(:, 3), 'k'); % Creates a patch using the optode circle on the head volume
    end

    % Group-Level Analysis
    scale(1) = mean(rotEig(1,1,:)); % Find the average scaling factor in the X-direction
    scale(2) = mean(rotEig(2,2,:)); % Find the average scaling factor in the Y-direction
    scale(3) = mean(rotEig(3,3,:)); % Find the average scaling factor in the Z-direction

    p2probeDesign(1:optN,1:3) = 0;
    p2error(1:optN,1:countSubjects) = 0;
    p2vec(1:optN,1:3,1:countSubjects) = 0;
    p2meanError(1:optN) = 0;
    
    for i = 1:1:optN
        for j = 1:1:3
            p2probeDesign(i,j) = atlas.probe.optpos_reg(i,j); % Original probe design
        end
        for j = 1:1:countSubjects
            p2error(i,j) = sqrt((((p2all(i,1,j)-p2probeDesign(i,1))) / scale(1))^2 + ...
                (((p2all(i,2,j)-p2probeDesign(i,2))) / scale(2))^2 + ...
                (((p2all(i,3,j)-p2probeDesign(i,3))) / scale(3))^2); % Calculate the Euclidean position error and correct for scaling
            p2vec(i,1,j) = (p2all(i,1,j)-p2probeDesign(i,1)) / scale(1);
            p2vec(i,2,j) = (p2all(i,2,j)-p2probeDesign(i,2)) / scale(2);
            p2vec(i,3,j) = (p2all(i,3,j)-p2probeDesign(i,3)) / scale(3);
        end
        p2meanError(i) = mean(p2error(i,:)); % Mean Error, scaled to millimeters
    end

    % Plot a second head and pial surface for group-level results
    figure('name','Mean Geometric Error from Original SD Locations (mm)');
    clf
    hold on
    c = v(:,3);
    c(:) = 100;
    set(gca,'clim',[100 200])
    h = trisurf(f, v(:,1), v(:,2), v(:,3), c(:));
    set(h, 'linestyle', 'none', 'facealpha', 0.25)
    h2 = trisurf(fp, vp(:,1), vp(:,2), vp(:,3), cp(:));
    set(h2, 'linestyle', 'none', 'facealpha', 1)
    lighting phong
    light('position',[-1 0 1])
    light('position',[0 -1 0])
    light('position',[1 0 -1])
    axis('off')
    axis image
    view(0, 0)
    
    % Create a color bar that displays the standard deviation value associated with each color
    colorbar('yticklabel', {0.10 * round(0*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(1*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(2*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(3*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(4*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(5*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(6*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(7*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(8*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(9*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError)), ...
        0.10 * round(10*((max(p2meanError)-min(p2meanError))))+0.10*round(10*min(p2meanError))});

    % Plot the original probe design (black)
    p2n(1:optN, 1:3) = 0; % Create an array to store the normal vector from the common atlas surface at each optode position
    p2 = atlas.probe.optpos_reg; % Original probe design
    for i=1:1:optN
        vmin = 1e99; vidx = 0;
        for j = 1:1:size(v, 1)
            foo = sqrt((p2(i,1)-v(j,1))^2 + (p2(i,2)-v(j,2))^2 + (p2(i,3)-v(j,3))^2); % distance from optode to vertex
            if (foo < vmin) % Search for vertex with the minimum distance from the optode
                vmin = foo;
                vidx = j;
            end
        end
        [lst, foo] = find(f==vidx); % Using the nearest vertex, find the matching faces
        for j = 1:1:size(lst, 1)
            p2n(i, :) = p2n(i, :) + fN(lst(j), :); % Sum the normal vectors for all matching faces
        end
        p2n(i, :) = p2n(i, :) / size(lst, 1); % Divide by the number of faces to find the average normal vector
        p2T = vrrotvec2mat_copy(vrrotvec_copy(circ1N, p2n(i, :))); % Create a rotation matrix using the original circle normal and the average normal vector for the optode location
        for j = 1:1:size(circle1, 1)
            circle2(j, :) = transpose(p2T * transpose(circle1(j, :))) + p2(i, 1:3); % Rotate and translate the original circle shape into the optode plane and position
        end
        patch(circle2(:, 1), circle2(:, 2), circle2(:, 3), 'k'); % Creates a patch using the optode circle on the head volume
    end
    
    p2mean(1:optN,1:3) = 0;
    for i = 1:1:optN
        for j = 1:1:3
            p2mean(i,j) = mean(p2all(i,j,:)); % Calculate the mean optode locations across subjects
        end
    end
    
    % Plot patches at mean of each probe location
    colors = colormap; % Create a matrix using the current color map for pulling color associated with each value
    p2n(1:optN, 1:3) = 0;
    for i=1:1:optN
        % Find the nearest vertex to each group-mean optode position
        vmin = 1e99; vidx = 0;
        for j = 1:1:size(v, 1)
            foo = sqrt((p2mean(i,1)-v(j,1))^2 + (p2mean(i,2)-v(j,2))^2 + (p2mean(i,3)-v(j,3))^2);
            if (foo < vmin)
                vmin = foo;
                vidx = j;
            end
        end
        [lst, foo] = find(f==vidx);
        % Find the normal vector for each group-mean optode position
        for j = 1:1:size(lst, 1)
            p2n(i, :) = p2n(i, :) + fN(lst(j), :);
        end
        p2n(i, :) = p2n(i, :) / size(lst, 1);
        
        p2proj(1:optN,1:3,1:countSubjects) = 0;
        for j=1:1:countSubjects
            p2projT = vrrotvec2mat_copy(vrrotvec_copy(p2n(i, :), [0 0 1]));
            p2proj(i,:,j) = transpose(p2projT * transpose(cross(p2n(i,:), cross(p2vec(i,:,j), p2n(i,:)))));
        end
        
        % Radius-defined ellipse for patches
        temp = 1:1:21;
        radiusX = std(p2proj(i,1,:));
        radiusY = std(p2proj(i,2,:));
        x = radiusX * cos(2 * pi * temp(:) / size(temp, 2));
        y = radiusY * sin(2 * pi * temp(:) / size(temp, 2));
        z(1:size(temp, 2)) = 0;
        circle1 = [x(:), y(:), z(:)];
        circ1N = [0 0 1];
        circle2 = circle1; % Create new variable with same size

        p2T = vrrotvec2mat_copy(vrrotvec_copy(circ1N, p2n(i, :)));
        for j = 1:1:size(circle1, 1)
            circle2(j, :) = transpose(p2T * transpose(circle1(j, :))) + p2mean(i, 1:3);
        end
        color = 1+round(63*(p2meanError(i)-min(p2meanError))/(max(p2meanError)-min(p2meanError))); % This picks the most appropriate of the 64 colormap colors
        patch(circle2(:, 1), circle2(:, 2), circle2(:, 3), ...
            [colors(color, 1) colors(color, 2) colors(color, 3)]); % Creates a patch using the colormap color and group-mean optode circle
    end
    
    % Create UI Table with Error Information
    f = figure('name', 'Probe Creation Error', 'OuterPosition', [300 300 651 500]);
    hold on
    errorTableData(:, 1:3) = p2(1:optN, 1:3);
    errorTableData(:, 4:6) = p2mean(1:optN, 1:3);
    errorTableData(:, 7) = transpose(p2meanError);
    uitable(f, 'Data', errorTableData, 'ColumnName', {'Design X', ...
        'Design Y', 'Design Z', 'Mean X', 'Mean Y', 'Mean Z', ...
        'Mean Geometric Error (mm)'}, 'units','normalized', 'Position',[0 0 1 1]);

else
    
    % Setup and Plot Atlas
    %atlas = load('atlasViewer.mat', '-mat'); % An unmodified atlasViewer.mat file 
    % for the common head volume. If possible, pull this from the Homer2 package.

    atlas = load(subFiles{1}, '-mat');
    T = atlas.headvol.T_2mc; % Translation matrix to individual subject space
    Ti(1:4, 1:4) = inv(T(1:4, 1:4)); % Inverse produces translation to common atlas space

    % Head Surface
    f = atlas.headsurf.mesh.faces;
    v = xform_apply(atlas.headsurf.mesh.vertices, Ti);
    c = v(:,3); % Creating a vector of equal length to v, to be used as color values
    c(:) = 100; % 100, using clim [100, 200], produces a blue color
    fN = trisurfnorm(v, f); % This function is provided in a separate file

    % Pial Surface
    fp = atlas.pialsurf.mesh.faces;
    vp = xform_apply(atlas.pialsurf.mesh.vertices, Ti);
    cp = vp(:,3); % Creating a vector of equal length to vp, to be used as color values
    cp(:) = 150; % 150, using clim [100, 200], produces a green color

    figure('name','Digitized Probe Positions');
    clf
    hold on
    set(gca,'clim',[100 200])
    h = trisurf(f, v(:,1), v(:,2), v(:,3), c(:)); % Plot head surface
    set(h, 'linestyle', 'none', 'facealpha', 0.25) % Remove lines and set opacity to 25%
    h2 = trisurf(fp, vp(:,1), vp(:,2), vp(:,3), cp(:)); % Plot pial surface
    set(h2, 'linestyle', 'none', 'facealpha', 1) % Remove lines and set opacity to 100%
    lighting phong
    light('position',[-1 0 1])
    light('position',[0 -1 0])
    light('position',[1 0 -1])
    axis('off')
    axis image
    view(0, 0)

    % Original Circle for Patches
    temp = 1:1:21;
    radius = 2.5;
    x = radius * cos(2 * pi * temp(:) / size(temp, 2));
    y = radius * sin(2 * pi * temp(:) / size(temp, 2));
    z(1:size(temp, 2)) = 0;
    circle1 = [x(:), y(:), z(:)];
    circ1N = [0 0 1]; % Normal to the plane of circle 1
    circle2 = circle1; % Create a second circle variable to be modified each loop of s (below)

    % Define the colors to be used for the optodes
    colors1 = colormap;
    for i = 1:1:countSubjects
        patchColor(i,:) = colors1(1+round(63*i/countSubjects),:);
    end

    % Translate each subject's optode positions onto a common atlas and plot
    rotEig(1:3,1:3,1:countSubjects) = 0;
    optN = 100; % Temporary value in case the following loop doesn't run
    nSources = 0;
    for s = 1:1:countSubjects
        % Load Subject Data, Transform Probe, Plot Circles Normal to Surface
        sub = load(subFiles{s}, '-mat');
        nSources = sub.probe.nsrc;
        p = sub.probe.optpos_reg; % Original registered optode positions
        p(:, 4) = 1;
        T = sub.headvol.T_2mc; % Translation matrix to individual subject space
        Ti(1:4, 1:4) = inv(T(1:4, 1:4)); % Inverse produces translation to common atlas space

        % Search for optodes that are closer together than the ssThresh and exclude them
        useList(1:size(p, 1)) = 1;
        for i = 1:1:size(p, 1)
            for j = 1:1:size(p, 1)
                if (sqrt((p(i, 1) - p(j, 1))^2 + (p(i, 2) - p(j, 2))^2 ...
                    + (p(i, 3) - p(j, 3))^2) < ssThresh)
                    % The following statements assume sources are listed 
                    % before detectors and therefore the optode to exclude 
                    % has the higher number
                    if (i > j)
                        useList(i) = 0;
                    end
                    if (j > i)
                        useList(j) = 0;
                    end
                end
            end
        end
        optN = sum(useList(:));
        p2(1:optN,1:4) = 0;
        j = 1;
        for i=1:1:size(p, 1)
            if (useList(i) ~= 0) % Only entries in the useList are carried to subsequent steps
                p2(j, 1:4) = transpose(Ti * transpose(p(i, 1:4))); % New optode positions in common space
                j = j + 1;
            end
        end

        p2n(1:optN, 1:3) = 0; % Create an array to store the normal vector from the common atlas surface at each optode position
        for i=1:1:optN
            vmin = 1e99; vidx = 0;
            for j = 1:1:size(v, 1)
                foo = sqrt((p2(i,1)-v(j,1))^2 + (p2(i,2)-v(j,2))^2 + (p2(i,3)-v(j,3))^2); % distance from optode to vertex
                if (foo < vmin) % Search for vertex with the minimum distance from the optode
                    vmin = foo;
                    vidx = j;
                end
            end
            [lst, foo] = find(f==vidx); % Using the nearest vertex, find the matching faces
            for j = 1:1:size(lst, 1)
                p2n(i, :) = p2n(i, :) + fN(lst(j), :); % Sum the normal vectors for all matching faces
            end
            p2n(i, :) = p2n(i, :) / size(lst, 1); % Divide by the number of faces to find the average normal vector
            p2T = vrrotvec2mat_copy(vrrotvec_copy(circ1N, p2n(i, :))); % Create a rotation matrix using the original circle normal and the average normal vector for the optode location
            for j = 1:1:size(circle1, 1)
                circle2(j, :) = transpose(p2T * transpose(circle1(j, :))) + p2(i, 1:3); % Rotate and translate the original circle shape into the optode plane and position
            end
            patch(circle2(:, 1), circle2(:, 2), circle2(:, 3), patchColor(s,:)); % Creates a patch using the optode circle on the head volume
        end

        [foo, bar] = eig(sub.headvol.T_2mc(1:3,1:3)); % Calculate the eigenvalues of the rotation matrix from subject space to common atlas space
        rotEig(1:3,1:3,s) = abs(bar); % Store the absolute value of the rotation matrix eigenvalues
        p2all(:,1:4,s) = p2; % Store the subject's translated optode locations for group plot
    end

    % Group-Level Analysis
    scale(1) = mean(rotEig(1,1,:)); % Find the average scaling factor in the X-direction
    scale(2) = mean(rotEig(2,2,:)); % Find the average scaling factor in the Y-direction
    scale(3) = mean(rotEig(3,3,:)); % Find the average scaling factor in the Z-direction

    p2mean(1:optN,1:3) = 0;
    p2error(1:optN,1:countSubjects) = 0;
    p2sd(1:optN) = 0;
    for i = 1:1:optN
        for j = 1:1:3
            p2mean(i,j) = mean(p2all(i,j,:)); % Calculate the mean optode locations across subjects
        end
        for j = 1:1:countSubjects
            p2error(i,j) = sqrt((((p2all(i,1,j)-p2mean(i,1))) / scale(1))^2 + ...
                (((p2all(i,2,j)-p2mean(i,2))) / scale(2))^2 + ...
                (((p2all(i,3,j)-p2mean(i,3))) / scale(3))^2); % Calculate the Euclidean position error and correct for scaling
        end
        for j = 1:1:countSubjects
            p2sd(i) = p2sd(i) + p2error(i,j)^2; % Accumulate Variance
        end
        p2sd(i) = sqrt(p2sd(i) / countSubjects); % Standard Deviation, scaled to millimeters
    end
    
    filename = fopen('probes_mean.txt', 'wt');
    for i = 1:1:nSources
        fprintf(filename, '%s%s%s %s %s %s\n', 's', num2str(i), ':', num2str(p2mean(i, 1)), num2str(p2mean(i, 2)), num2str(p2mean(i, 3)));
    end
    for i = 1:1:(optN - nSources)
        fprintf(filename, '%s%s%s %s %s %s\n', 'd', num2str(i), ':', num2str(p2mean(i + nSources, 1)), num2str(p2mean(i + nSources, 2)), num2str(p2mean(i + nSources, 3)));
    end
    fprintf(filename, '\nNz: 131.0 188.0 234.0\n');
    fprintf(filename, 'Iz: 130.0 167.0 26.0\n');
    fprintf(filename, 'Ar: 40.0 201.0 124.0\n');
    fprintf(filename, 'Al: 218.0 203.0 124.0\n');
    fprintf(filename, 'Cz: 130.0 47.9 138.0\n');
    fclose(filename);

    % Plot a second head and pial surface for group-level results
    figure('name','Standard Deviation Across Digitized Probes (mm)');
    clf
    hold on
    c = v(:,3);
    c(:) = 100;
    set(gca,'clim',[100 200])
    h = trisurf(f, v(:,1), v(:,2), v(:,3), c(:));
    set(h, 'linestyle', 'none', 'facealpha', 0.25)
    h2 = trisurf(fp, vp(:,1), vp(:,2), vp(:,3), cp(:));
    set(h2, 'linestyle', 'none', 'facealpha', 1)
    lighting phong
    light('position',[-1 0 1])
    light('position',[0 -1 0])
    light('position',[1 0 -1])
    axis('off')
    axis image
    view(0, 0)

    % Create a color bar that displays the standard deviation value associated with each color
    colorbar('yticklabel', {0.10 * round(0*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(1*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(2*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(3*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(4*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(5*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(6*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(7*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(8*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(9*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd)), ...
        0.10 * round(10*((max(p2sd)-min(p2sd))))+0.10*round(10*min(p2sd))});

    % Plot patches at mean of each probe location
    colors = colormap; % Create a matrix using the current color map for pulling color associated with each value
    p2n(1:optN, 1:3) = 0;
    for i=1:1:optN
        % Find the nearest vertex to each group-mean optode position
        vmin = 1e99; vidx = 0;
        for j = 1:1:size(v, 1)
            foo = sqrt((p2mean(i,1)-v(j,1))^2 + (p2mean(i,2)-v(j,2))^2 + (p2mean(i,3)-v(j,3))^2);
            if (foo < vmin)
                vmin = foo;
                vidx = j;
            end
        end
        [lst, foo] = find(f==vidx);
        % Fine the normal vector for each group-mean optode position
        for j = 1:1:size(lst, 1)
            p2n(i, :) = p2n(i, :) + fN(lst(j), :);
        end
        p2n(i, :) = p2n(i, :) / size(lst, 1);

        p2T = vrrotvec2mat_copy(vrrotvec_copy(circ1N, p2n(i, :)));
        for j = 1:1:size(circle1, 1)
            circle2(j, :) = transpose(p2T * transpose(circle1(j, :))) + p2(i, 1:3);
        end
        
        color = 1+round(63*(p2sd(i)-min(p2sd))/(max(p2sd)-min(p2sd))); % This picks the most appropriate of the 64 colormap colors
        patch(circle2(:, 1), circle2(:, 2), circle2(:, 3), ...
            [colors(color, 1) colors(color, 2) colors(color, 3)]); % Creates a patch using the colormap color and group-mean optode circle
    end
    
    % Create UI Table with Mean Probe Information
    f = figure('name', 'Mean Optode Positions', 'OuterPosition', [300 300 412 500]);
    hold on
    tableData = p2mean;
    tableData(:, 4) = transpose(p2sd);
    uitable(f, 'Data', tableData, 'ColumnName', {'X', 'Y', 'Z', 'Standard Deviation (mm)'}, ...
        'units','normalized', 'Position',[0,0,1,1]);
end

cd(wd)
