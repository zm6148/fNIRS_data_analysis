function fwmodel = getFwmodel(fwmodel, dirname, pialsurf, headsurf, headvol, probe)

if iscell(dirname)
    for ii=1:length(dirname)
        fwmodel = getFwmodel(fwmodel, dirname{ii}, pialsurf, headsurf, headvol, probe);
        if ~fwmodel.isempty(fwmodel)
            return;
        end
        if ~isempty(fwmodel.projVoltoMesh_brain)
            return;
        end
    end
    return;
end

if isempty(dirname)
    return;
end

if dirname(end)~='/' && dirname(end)~='\'
    dirname(end+1)='/';
end
dirnameOut = [dirname, 'fw/'];

if exist([dirnameOut, 'headvol.vox'],'file')
    fwmodel.headvol = load_vox([dirnameOut, 'headvol.vox'], fwmodel.headvol);
    if ~isempty(headvol.center)
        fwmodel.headvol.center = xform_apply(headvol.center, headvol.T_2mc);
    else
        fwmodel.headvol.center = findcenter(fwmodel.headvol.img);
    end
elseif isidentity(headvol.T_2digpts)
    fwmodel.headvol = headvol;
end

if exist([dirnameOut, 'pialsurf_sensitivity.mesh'],'file')
    [v,f] = read_surf([dirnameOut, 'pialsurf_sensitivity.mesh']);
    fwmodel.mesh.vertices = v;
    fwmodel.mesh.faces = f;
else
    fwmodel.mesh = pialsurf.mesh;
end

if exist([dirnameOut, 'headsurf_sensitivity.mesh'],'file')
    [v,f] = read_surf([dirnameOut, 'headsurf_sensitivity.mesh']);
    fwmodel.mesh_scalp.vertices = v;
    fwmodel.mesh_scalp.faces = f;
else
    fwmodel.mesh_scalp = headsurf.mesh;
end

if exist([dirnameOut, 'projVoltoMesh_brain.mat'],'file')
    fwmodel.projVoltoMesh_brain = [dirnameOut, 'projVoltoMesh_brain.mat'];
end

if exist([dirnameOut, 'projVoltoMesh_scalp.mat'],'file')
    fwmodel.projVoltoMesh_scalp = [dirnameOut, 'projVoltoMesh_scalp.mat'];
end


if isempty(fwmodel.fluenceProfFnames)
    
    % Check if there a fluence profile to load in this particular search path
    fluenceProfFnames = dir([dirnameOut, 'fluenceProf*.mat']);
    for ii=1:length(fluenceProfFnames)
        foo = loadFluenceProf([dirnameOut, fluenceProfFnames(ii).name], 'index');
        fwmodel.fluenceProfFnames{foo.index} = [dirnameOut, fluenceProfFnames(ii).name];
    end
    
    if ~isempty(fwmodel.fluenceProfFnames)
        s = load(fwmodel.fluenceProfFnames{1});
        fwmodel.mesh = s.mesh;
        fwmodel.nWavelengths = size(s.intensities,3);
        fwmodel.nphotons = s.nphotons;
        fwmodel.headvol.tiss_prop = s.tiss_prop;
        fwmodel.nFluenceProfPerFile = size(s.srcpos,1);
        fwmodel.voxPerNode = s.voxPerNode;
    end
    
end
    
if exist([dirnameOut, 'Adot.mat'],'file')
    
    load([dirnameOut, 'Adot.mat']);
    
    if exist('tiss_prop','var')
        fwmodel.headvol.tiss_prop = tiss_prop;
    end
    if exist('nphotons','var')
        fwmodel.nphotons = nphotons;
    end
    if ~isempty(fwmodel.mesh)
        if size(fwmodel.mesh.vertices,1) == size(Adot,2)
            fwmodel.Adot = Adot;
            d=getFileDateStruct([dirnameOut, 'Adot.mat']);
            fwmodel.AdotDate = d;
            fwmodel.nWavelengths = size(Adot,3);
        end
    end
    
    if exist([dirnameOut, 'Adot_scalp.mat'],'file')
        load([dirnameOut, 'Adot_scalp.mat']);        
        if ~isempty(fwmodel.mesh_scalp)
            if size(fwmodel.mesh_scalp.vertices,1) == size(Adot_scalp,2)
                fwmodel.Adot_scalp = Adot_scalp;
                d=getFileDateStruct([dirnameOut, 'Adot_scalp.mat']);
            end
        end        
    else
        fwmodel.Adot_scalp     = [];        
    end
    
    set(fwmodel.handles.menuItemImageReconGUI,'enable','on');
    
else
    
    fwmodel.Adot     = [];
    fwmodel.AdotDate = struct('num',0);
    set(fwmodel.handles.menuItemImageReconGUI,'enable','off');
    
    if exist([dirnameOut, 'headvol.vox'], 'file')
        headvol = load_vox([dirnameOut, 'headvol.vox']);
    end
    
    % Wavelength 1
    if exist([dirnameOut, 'fw1.s1.inp'],'file')
        config = read_tMCimg_inp([dirnameOut, 'fw1.s1.inp']);
        fwmodel.nphotons = config.phot_num;
        for ii=1:size(config.tiss_prop,1)
            fwmodel.headvol.tiss_prop(ii).name = headvol.tiss_prop(ii).name; 
            fwmodel.headvol.tiss_prop(ii).scattering(:,1) = config.tiss_prop(ii,1);
            fwmodel.headvol.tiss_prop(ii).anisotropy(:,1) = config.tiss_prop(ii,2);
            fwmodel.headvol.tiss_prop(ii).absorption(:,1) = config.tiss_prop(ii,3);
            fwmodel.headvol.tiss_prop(ii).refraction(:,1) = config.tiss_prop(ii,4);
        end
        fwmodel.nWavelengths = 1;
    end
    
    % Wavelength 2
    if exist([dirnameOut, 'fw2.s1.inp'],'file')
        config = read_tMCimg_inp([dirnameOut, 'fw2.s1.inp']);
        fwmodel.nphotons = config.phot_num;
        for ii=1:size(config.tiss_prop,1)
            fwmodel.headvol.tiss_prop(ii).name = headvol.tiss_prop(ii).name; 
            fwmodel.headvol.tiss_prop(ii).scattering(:,2) = config.tiss_prop(ii,1);
            fwmodel.headvol.tiss_prop(ii).anisotropy(:,2) = config.tiss_prop(ii,2);
            fwmodel.headvol.tiss_prop(ii).absorption(:,2) = config.tiss_prop(ii,3);
            fwmodel.headvol.tiss_prop(ii).refraction(:,2) = config.tiss_prop(ii,4);
        end
        fwmodel.nWavelengths = 2;
    end

end

if ~fwmodel.isempty(fwmodel)
    fwmodel.pathname = dirname;
end


