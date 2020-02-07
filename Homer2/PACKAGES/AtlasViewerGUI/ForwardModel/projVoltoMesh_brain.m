function [mapMesh2Vox, fwmodel] = projVoltoMesh_brain(fwmodel, dirnameOut)

% Project cortical indexed voxels onto pial meshes.
%
% Written by Matteo Caffini
% Modified by Jay Dubb 
%
% usage:
% call function projVoltoMesh in subject directory or registered atlas
% directory.

%  load the segmented Vol

hf = [];

if isempty(fwmodel.projVoltoMesh_brain)
    
    %%%%%%%%%%%%%%%%%%%%
    % load head volume
    %%%%%%%%%%%%%%%%%%%%
    nx = size(fwmodel.headvol.img,1);
    ny = size(fwmodel.headvol.img,2);
    nz = size(fwmodel.headvol.img,3);
    
    [nodeX, nNode, fwmodel.mesh] = reduceMesh(fwmodel, dirnameOut);
    
    if ishandles(hf)
        delete(hf);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % map Vol to LR Mesh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check if there's a segmentation file telling us 
    % what the tissue number of gray matter is in the 
    % seg file. 
    for gm_seg_num=1:length(fwmodel.headvol.tiss_prop)
        if(strcmp(lower(fwmodel.headvol.tiss_prop(gm_seg_num).name), 'gm') | ...
           strcmp(lower(fwmodel.headvol.tiss_prop(gm_seg_num).name), 'gray matter'))
            break;
        end
    end
    i_headvol = uint32(find(fwmodel.headvol.img==gm_seg_num));
    nmiss=0;
    nxy = nx*ny;
    nC = length(i_headvol);
    Amap = single(zeros(nC,1));
    mapMesh2Vox = single(ones(nNode,1000));
    NVoxPerNode = zeros(nNode,1);
    
    % We don't want to slow things down with too frequent updates
    update_interval=ceil(nC/100);  
    hwait = waitbar(0,'Looping over cortical voxels');
    h = 15; % with 15 we miss < 1% of total number of cortex voxels
    for ii=1:nC
        if mod(ii,update_interval)==1
            waitbar(ii/nC,hwait,sprintf('%d of %d',ii,nC));
        end
        [x,y,z] = ind2sub(size(fwmodel.headvol.img),i_headvol(ii));
        x = double(x);
        y = double(y);
        z = double(z);
        
        i_nX = find(nodeX(:,1)>x-h & nodeX(:,1)<x+h & nodeX(:,2)>y-h & nodeX(:,2)<y+h & ...
                    nodeX(:,3)>z-h & nodeX(:,3)<z+h);
        if ~isempty(i_nX)
            % rsep: get the distances from [x y z] to all the points in nodeX(i_nX,:).
            rsep = sum( (nodeX(i_nX,:) - ones(length(i_nX),1)*[x y z]).^2, 2 ).^0.5;
            [foo,imin] = min(rsep);
            Amap(ii) = i_nX(imin);
            NVoxPerNode(Amap(ii)) = NVoxPerNode(Amap(ii))+1; % might be useful?
            mapMesh2Vox(Amap(ii),NVoxPerNode(Amap(ii))) = i_headvol(ii);
        else
            nmiss = nmiss+1; % temporary var, delete when everything works
        end
    end;
    close(hwait);
    ciao = find(mapMesh2Vox(:)==0);
    mapMesh2Vox(ciao) = 1;

    save([dirnameOut, 'projVoltoMesh_brain.mat'], 'mapMesh2Vox');
    fwmodel.projVoltoMesh_brain = [dirnameOut, 'projVoltoMesh_brain.mat'];
    
else
    
    load(fwmodel.projVoltoMesh_brain);
   
end




% --------------------------------------------------------------------------
function [nodeX, nNode, mesh] = reduceMesh(fwmodel, dirnameOut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if number of elements is too large. If greater than 40,000 then
% need to reduce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mesh  = fwmodel.mesh;
nodeX = mesh.vertices;
elem  = mesh.faces;  % Let's put a check in here for number of faces. Good to be less than 70,000. Or should this be at getPialsurf?
nNode = size(nodeX,1);

if size(elem,1)>40000
    idx = menu( sprintf('The pial surface has %d faces. It is recommended that this be less than 40,000.\nShall I reduce this?',size(elem,1)),'Yes','No');
    if idx==1
        hwait = waitbar(0,'Reducing mesh...');
        fv = mesh;
        fvn = reducepatch( fv, 40000/size(elem,1) );
        close(hwait)
        
        hf = figure;
        subplot(1,2,1)
        h=trisurf( fv.faces, fv.vertices(:,1), fv.vertices(:,2), fv.vertices(:,3) );
        title( sprintf(' Original with %d faces', size(fv.faces,1)) )
        set(h,'linestyle','none')
        light
        subplot(1,2,2)
        h=trisurf( fvn.faces, fvn.vertices(:,1), fvn.vertices(:,2), fvn.vertices(:,3) );
        title( sprintf(' New with %d faces', size(fvn.faces,1)) )
        set(h,'linestyle','none')
        light
        
        idx = menu('Accept this?','Yes','No');
        if idx==1
            mesh = fvn;
            nodeX = fvn.vertices;
            elem = fvn.faces;
            nNode = size(nodeX,1);
            write_surf([dirnameOut, 'pialsurf_sensitivity.mesh'], fvn.vertices, fvn.faces);
        else
            idx = menu('Proceed with the original mesh? This could take a long time.','Yes','No');
            if idx==2
                mapMesh2Vox = [];
                mesh = initMesh();
                if ishandles(hf)
                    delete(hf);
                end
                return
            end
        end
    end
end
