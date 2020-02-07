function varargout = AtlasViewerGUI(varargin)

% Start AtlasViewerGUI initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AtlasViewerGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AtlasViewerGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
if nargin && ischar(varargin{1}) && ~strcmp(varargin{end},'userargs')
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End AtlasViewerGUI initialization code - DO NOT EDIT

        
        
% ------------------------------------------------------------------
function InitSubj(hObject,handles,argExtern)
global atlasViewer
global DEBUG

DEBUG = 0;

%%% Begin initialization ....

% Create things from scratch

if isempty(argExtern)
    argExtern = {''};
end

handles.ImageRecon = [];
if length(argExtern)>=4
    if length(argExtern{4})>1
        handles.ImageRecon = argExtern{4}(2);
    end
end

dirnameSubj = getSubjDir(argExtern);
dirnameAtlas = getAtlasDir(argExtern);

fprintf('In AtlasViewerGUI_Init:\n');
fprintf('   dirnameAtlas = %s\n', dirnameAtlas);
fprintf('   dirnameSubj = %s\n', dirnameSubj);

cd(dirnameSubj);

atlasViewer.handles.figure = hObject;

% Initialize atlas viewer objects with their respective gui
% handles
objs.axesv       = initAxesv(handles);
objs.headsurf    = initHeadsurf(handles);
objs.pialsurf    = initPialsurf(handles);
objs.labelssurf  = initLabelssurf(handles);
objs.refpts      = initRefpts(handles);
objs.digpts      = initDigpts(handles);
objs.headvol     = initHeadvol();
objs.probe       = initProbe(handles);
objs.fwmodel     = initFwmodel(handles, argExtern);
objs.imgrecon    = initImgRecon(handles);
objs.fs2viewer   = initFs2Viewer(handles,dirnameSubj);

fprintf('   MC application path = %s\n', objs.fwmodel.mc_exepath);
fprintf('   MC application binary = %s\n', objs.fwmodel.mc_exename);

fields = fieldnames(objs);

% Check for a saved viewer state file and restore 
% state if it exists. 
vrnum = [];
if exist([dirnameSubj 'atlasViewer.mat'], 'file')

    load([dirnameSubj 'atlasViewer.mat'],'-mat');
    for ii=1:length(fields)
        if exist(fields{ii},'var')
            % Initialized object exists in saved state. Check its compatibility with current version
            eval(sprintf('b = ~isempty(objs.%s.checkCompatability);', fields{ii}));
            if b==1
                eval(sprintf('%s = objs.%s.checkCompatability(%s);', fields{ii}, fields{ii}, fields{ii}));
            end
            eval(sprintf('atlasViewer.%s = restoreObject(%s, objs.%s);', fields{ii}, fields{ii}, fields{ii}));
        else
            % Initialized object does NOT exist in saved state. Therefore no compatibility issues.  
            eval(sprintf('atlasViewer.%s = restoreObject(objs.%s, objs.%s);', fields{ii}, fields{ii}, fields{ii}));
        end
    end
    if ~isempty(vrnum)
        fprintf('Loading saved viewer state created by AtlasViewerGUI V%s\n', vrnum);
    else
        fprintf('Loading saved viewer state created by a version of AtlasViewerGUI prior to V2.0.1\n');
    end
    
    % Otherwise simply initialize objects from scratch

else
    
    for ii=1:length(fields)
        eval(sprintf('atlasViewer.%s = objs.%s;', fields{ii}, fields{ii}));            
    end
    
end

atlasViewer.dirnameAtlas = dirnameAtlas;
atlasViewer.dirnameSubj  = dirnameSubj;
atlasViewer.dirnameProbe = '';
atlasViewer.groupSubjList = '';
atlasViewer.handles.menuItemRegisterAtlasToDigpts = handles.menuItemRegisterAtlasToDigpts;

% Set the AtlasViewerGUI version number
V = AtlasViewerGUI_version();
if str2num(V{2})==0
    set(hObject,'name', sprintf('AtlasViewerGUI  (v%s) - %s', [V{1}],cd) )
else
    set(hObject,'name', sprintf('AtlasViewerGUI  (v%s) - %s', [V{1} '.' V{2}],cd) )
end
atlasViewer.vrnum = V;




% -----------------------------------------------------------------------
function LoadSubj(hObject, eventdata, handles, argExtern)
global atlasViewer

if isempty(argExtern)
    argExtern = {''};
end

InitSubj(hObject,handles,argExtern);

dirnameAtlas = atlasViewer.dirnameAtlas;
dirnameSubj = atlasViewer.dirnameSubj;
searchPaths = {dirnameSubj; dirnameAtlas};

axesv        = atlasViewer.axesv;
headvol      = atlasViewer.headvol;
headsurf     = atlasViewer.headsurf;
pialsurf     = atlasViewer.pialsurf;
labelssurf   = atlasViewer.labelssurf;
refpts       = atlasViewer.refpts;
digpts       = atlasViewer.digpts;
probe        = atlasViewer.probe;
fwmodel      = atlasViewer.fwmodel;
imgrecon     = atlasViewer.imgrecon;
fs2viewer    = atlasViewer.fs2viewer;
    

if ~exist([dirnameSubj 'atlasViewer.mat'], 'file')

    % Load all objects
    digpts     = getDigpts(digpts, dirnameSubj);
    headvol    = getHeadvol(headvol, searchPaths);
    headsurf   = getHeadsurf(headsurf, headvol.pathname);
    refpts     = getRefpts(refpts, headvol.pathname);
    pialsurf   = getPialsurf(pialsurf, headvol.pathname);
    labelssurf = getLabelssurf(labelssurf, headvol.pathname);
    probe      = getProbe(probe, dirnameSubj, headsurf);
    fwmodel    = getFwmodel(fwmodel, dirnameSubj, pialsurf, headsurf, headvol, probe);
    imgrecon   = getImgRecon(imgrecon, dirnameSubj, fwmodel, pialsurf, probe);
    fs2viewer  = getFs2Viewer(fs2viewer, dirnameSubj, headsurf, headvol, pialsurf);
        
end

% Set orientation and main axes attributes for all objects
if ~refpts.isempty(refpts)
    [headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon] = ...
        setOrientationRefpts(refpts, headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon);
elseif ~headvol.isempty(headvol)
    [refpts, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon] = ...
        setOrientationHeadvol(headvol, refpts, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon);
end

% Display all objects
digpts     = displayDigpts(digpts);
probe      = displayProbe(probe, headsurf);
if isregistered(refpts,digpts) | digpts.isempty(digpts)
    refpts     = displayRefpts(refpts);
    headsurf   = displayHeadsurf(headsurf);
    pialsurf   = displayPialsurf(pialsurf);
    labelssurf = displayLabelssurf(labelssurf);
    fwmodel    = displaySensitivity(fwmodel, pialsurf, labelssurf, probe);
    imgrecon   = displayImgRecon(imgrecon, fwmodel, pialsurf, labelssurf, probe);
    axesv      = displayAxesv(axesv, headsurf, initDigpts());
else
    axesv      = displayAxesv(axesv, initHeadsurf(), digpts);
end

atlasViewer.headsurf    = headsurf;
atlasViewer.pialsurf    = pialsurf;
atlasViewer.labelssurf  = labelssurf;
atlasViewer.refpts      = refpts;
atlasViewer.headvol     = headvol;
atlasViewer.digpts      = digpts;
atlasViewer.probe       = probe;
atlasViewer.fwmodel     = fwmodel;
atlasViewer.imgrecon    = imgrecon;
atlasViewer.axesv       = axesv;
atlasViewer.fs2viewer   = fs2viewer;

% TBD: The workflow should be to go from volume space 
% to dig point to monte carlo space automatically whether
% dig pts are present or not - in which case the transformations 
% are identities and volume space is MC space. 
% Right now calling menuItemRegisterAtlasToDigpts_Callback 
% explicitly that is by a non-graphics event 
% doesn't do anything. It's more like a placeholder. 
menuItemRegisterAtlasToDigpts_Callback();


% Enable menu items 
AtlasViewerGUI_enableDisable();



% ---------------------------------------------------------------------
function [refpts, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon] = ...
     setOrientationHeadvol(headvol, refpts, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon)

if headvol.isempty(headvol)
    return;
end
if isempty(headvol.orientation)
    return;
end

refpts.orientation     = headvol.orientation;
refpts.center          = headvol.center;

headsurf.orientation   = headvol.orientation;
headsurf.center        = headvol.center;

pialsurf.orientation   = headvol.orientation;
pialsurf.center        = headvol.center;

labelssurf.orientation = headvol.orientation;
labelssurf.center      = headvol.center;

if isempty(probe.orientation)
    probe.orientation      = headvol.orientation;
    probe.center           = headvol.center;
end

fwmodel.orientation    = headvol.orientation;
fwmodel.center         = headvol.center;

imgrecon.orientation   = headvol.orientation;
imgrecon.center        = headvol.center;



% ---------------------------------------------------------------------
function [headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon] = ...
    setOrientationRefpts(refpts, headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon)

if refpts.isempty(refpts)
    return;
end

if isempty(refpts.orientation)
    [nz,iz,rpa,lpa,cz] = getLandmarks(refpts);
    [refpts.orientation, refpts.center]  = getOrientation(nz,iz,rpa,lpa,cz);
end 

headvol                = saveHeadvolOrient(headvol, refpts);

headsurf.orientation   = refpts.orientation;
headsurf.center        = refpts.center;

pialsurf.orientation   = refpts.orientation;
pialsurf.center        = refpts.center;

labelssurf.orientation = refpts.orientation;
labelssurf.center      = refpts.center;

if isempty(probe.orientation)
    probe.orientation      = refpts.orientation;
    probe.center           = refpts.center;
end

fwmodel.orientation    = refpts.orientation;
fwmodel.center         = refpts.center;

imgrecon.orientation   = refpts.orientation;
imgrecon.center        = refpts.center;




% ---------------------------------------------------------------------
function [refpts, headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon] = ...
    setOrientationDigpts(digpts, refpts, headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon)

if isempty(digpts.orientation)
    return;
end

refpts.orientation     = digpts.orientation;
refpts.center          = digpts.center;

headvol                = saveHeadvolOrient(headvol, digpts);

headsurf.orientation   = digpts.orientation;
headsurf.center        = digpts.center;

pialsurf.orientation   = digpts.orientation;
pialsurf.center        = digpts.center;

labelssurf.orientation = digpts.orientation;
labelssurf.center      = digpts.center;

probe.orientation      = digpts.orientation;
probe.center           = digpts.center;

fwmodel.orientation    = digpts.orientation;
fwmodel.center         = digpts.center;

imgrecon.orientation   = digpts.orientation;
imgrecon.center        = digpts.center;




% ---------------------------------------------------------------------
function [groupSubjList, dirname] = InitGroup(argExtern)

if isempty(argExtern)
    argExtern = {''};
end

dirname = getSubjDir(argExtern);

groupSubjList = {};

% To find out the subj
[subjDirs, groupDir] = findSubjDirs();
if isempty(groupDir)
    [subjDirs, groupDir] = findSubjDirs('../');
end

if isempty(groupDir)
    return;
end
groupSubjList{1} = groupDir;
for ii=2:length(subjDirs)+1
    groupSubjList{ii} = [groupDir, '/', subjDirs(ii-1).name];
end



% -----------------------------------------------------------------------
function groupSubjList_Callback(hObject,eventdata,handles)
global atlasViewer

if isempty(atlasViewer)
    return;
end

dirnameAtlas = atlasViewer.dirnameAtlas;
axesv = atlasViewer.axesv;
fwmodel = atlasViewer.fwmodel;
imgrecon = atlasViewer.imgrecon;

groupSubjList = getappdata(hObject, 'groupSubjList');
idx = get(hObject,'value');

if isempty(groupSubjList)
    return;
end
if idx>length(groupSubjList)
    return;
end
dirnameSubj = groupSubjList{idx};
fprintf('Loading subject %s ...\n', dirnameSubj);
if dirnameSubj==0
    return;
end

hImageRecon = imgrecon.handles.ImageRecon;
set(hObject,'enable','off');
AtlasViewerGUI(dirnameSubj, dirnameAtlas, fwmodel.mc_exepath, [hObject, hImageRecon], 'userargs');




% -----------------------------------------------------------------------
function hGroupList = displayGroupSubjList(groupSubjList0, hGroupList, hGui)

if isempty(groupSubjList0)
    if ishandles(hGroupList)
        delete(hGroupList);
    end
    return;
end

groupSubjList = {};
for ii=1:length(groupSubjList0)
    [pp,fp] = getpathparts(groupSubjList0{ii});
    
    subjListboxStr = pp{end};
    if ii>1
        subjListboxStr = ['  ', pp{end}];
    end
    groupSubjList{ii} = subjListboxStr;
end

if ishandles(hGroupList)
    hFig = get(hGroupList,'parent');
    set(hGroupList, 'string',groupSubjList);
    figure(hFig);
    return;
end
if ~exist('hGui','var')
    hGui = [];
end
if ishandles(hGui)
    set(hGui,'units','normalized','position',[.33,.08,.66,.85]);
end

hFig = figure('numbertitle','off','menubar','none','name','Group Subject List','units','normalized',...
              'position',[.05,.45,.15,.40],'resize','on', 'visible','off');

hGroupList = uicontrol('parent',hFig,'style','listbox','string',groupSubjList,...
                       'fontsize',10,'units','normalized','position',[.10 .25 .80 .70],'value',1,...
                       'callback',{@groupSubjList_Callback,'hObject'});

setappdata(hGroupList, 'groupSubjList', groupSubjList0);

% Initilize listbox selection to current subject folder
[pp,fs] = getpathparts(pwd);
subjname = pp{end};
k =  find(strcmp(strtrim(groupSubjList), subjname));
if ~isempty(k)
    set(hGroupList, 'value', k);
end





% -----------------------------------------------------------------------
function AtlasViewerGUI_OpeningFcn(hObject, eventdata, handles, varargin)
global atlasViewer
atlasViewer = [];

% Choose default command line output for AtlasViewerGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

initAxesv(handles);

[groupSubjList, dirnameSubj] = InitGroup(varargin);
hGroupList=[];
if length(varargin)>3
    if ~isempty(varargin{4})
        hGroupList = varargin{4}(1);
    end
end
handles.hGroupList = displayGroupSubjList(groupSubjList, hGroupList, hObject);

if ~isempty(dirnameSubj) & dirnameSubj ~= 0
    if length(varargin)<2
        varargin{1} = dirnameSubj;
        varargin{2} = 'userargs';
    else
        varargin{1} = dirnameSubj;
    end
end
LoadSubj(hObject, eventdata, handles, varargin);

fprintf('Subject index = %d\n', atlasViewer.imgrecon.iSubj);

if ishandles(handles.hGroupList)
    set(handles.hGroupList, 'enable','on');
    hParent = get(handles.hGroupList,'parent');
    set(hParent, 'visible','on');
    figure(hParent);
end

atlasViewer.handles.hGroupList = handles.hGroupList;
atlasViewer.groupSubjList = groupSubjList;

if ishandles(atlasViewer.imgrecon.handles.ImageRecon)
    ImageRecon();
end



% -------------------------------------------------------------------
function varargout = AtlasViewerGUI_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------------------------------------------------------
function AtlasViewerGUI_DeleteFcn(hObject, eventdata, handles)
global atlasViewer

fclose all;
if isempty(atlasViewer)
    return;
end
axesv = atlasViewer.axesv;

probe = atlasViewer.probe;
imgrecon = atlasViewer.imgrecon;

if ishandles(probe.handles.hSDgui)
    delete(probe.handles.hSDgui);
end

if ishandles(imgrecon.handles.ImageRecon)
    delete(imgrecon.handles.ImageRecon);
end

if length(axesv)>1
    if ishandles(axesv(2).handles.axesSurfDisplay)
        hp = get(axesv(2).handles.axesSurfDisplay,'parent');
        delete(hp);
    end
end
if ishandles(atlasViewer.handles.hGroupList)
    hFig = get(atlasViewer.handles.hGroupList,'parent');
    delete(hFig);
end
atlasViewer=[];
clear atlasViewer;




% ----------------------------------------------------------------
function radiobuttonShowHead_Callback(hObject, eventdata, handles)
global atlasViewer;
hHeadSurf = atlasViewer.headsurf.handles.surf;

val = get(hObject,'value');

if val==1
    set(hHeadSurf,'visible','on');
elseif val==0
    set(hHeadSurf,'visible','off');
end




% ------------------------------------------------------------------
function editHeadTransparency_Callback(hObject, eventdata, handles)
global atlasViewer;

hHeadSurf = atlasViewer.headsurf.handles.surf;

val_old = get(hHeadSurf,'facealpha');
val = str2num(get(hObject,'string'));

% Error checking 
if isempty(val)
    set(hObject,'string',num2str(val_old));
    return;
end
if ~isnumeric(val)
    set(hObject,'string',num2str(val_old));
    return;
end
if ~isscalar(val)
    set(hObject,'string',num2str(val_old));
    return;
end
if val>1 || val<0
    set(hObject,'string',num2str(val_old));
    return;
end
if ~isempty(hHeadSurf)
    set(hHeadSurf,'facealpha',val);
end




% --------------------------------------------------------------------
function editBrainTransparency_Callback(hObject, eventdata, handles)
global atlasViewer;

hPialSurf = atlasViewer.pialsurf.handles.surf;
hLabelsSurf = atlasViewer.labelssurf.handles.surf;
iFaces = atlasViewer.labelssurf.iFaces;
mesh = atlasViewer.labelssurf.mesh;

val_old = get(hPialSurf,'facealpha');
val = str2num(get(hObject,'string'));

% Error checking 
if isempty(val)
    set(hObject,'string',num2str(val_old));
    return;
end
if ~isnumeric(val)
    set(hObject,'string',num2str(val_old));
    return;
end
if ~isscalar(val)
    set(hObject,'string',num2str(val_old));
    return;
end
if val>1 || val<0
    set(hObject,'string',num2str(val_old));
    return;
end

if ~isempty(hPialSurf)
    set(hPialSurf,'facealpha',val);
end
if ~isempty(hLabelsSurf)
    facevertexalphadata = ones(size(mesh.faces,1),1)*val;
    facevertexalphadata(iFaces) = 1;
    set(hLabelsSurf,'facevertexalphadata',facevertexalphadata);
end




% --------------------------------------------------------------------
function fwmodel = updateGuiControls_AfterProbeRegistration(probe, fwmodel, imgrecon, labelssurf, dirnameSubj)


if ishandles(labelssurf.handles.surf)
    set(probe.handles.menuItemProjectOptodesToCortex, 'enable','on');
    set(probe.handles.menuItemProjectChannelsToCortex, 'enable','on');
else
    set(probe.handles.menuItemProjectOptodesToCortex, 'enable','off');
    set(probe.handles.menuItemProjectChannelsToCortex, 'enable','off');
end

if ~isempty(probe.optpos_reg)
    enableMCGenGuiControls(fwmodel, 'on');
else
    enableMCGenGuiControls(fwmodel, 'off');
end

if ~isempty(probe.ml)
    if isempty(fwmodel.Adot) 
        if ~isempty(fwmodel.fluenceProfFnames)
            enableDisableMCoutputGraphics(fwmodel, 'on');
        else
            enableDisableMCoutputGraphics(fwmodel, 'off');
        end    
        enableImgReconGen(imgrecon, 'off');
        enableImgReconDisplay(imgrecon, 'off');
    else
        if size(probe.ml,1) ~= size(fwmodel.Adot,1)
            fwmodel.Adot = [];
            onoff = 'off';
        else
            onoff = 'on';            
        end
        enableDisableMCoutputGraphics(fwmodel, onoff);
        enableImgReconGen(imgrecon, onoff);
        enableImgReconDisplay(imgrecon, onoff);
    end
else
    enableDisableMCoutputGraphics(fwmodel, 'off');
    enableImgReconGen(imgrecon, 'off');
    enableImgReconDisplay(imgrecon, 'off');
end





% --------------------------------------------------------------------
function probe = probeRegisterSpringsMethod(probe,headvol,refpts)

if isempty(probe)
    menu('probe hasn''t been loaded. Use the Make Probe option in the Tools menu','OK');
    return;
end
if isempty(probe.optpos)
    menu('No source/detector positions. Use Make Probe in the Tools menu','OK');
    return;
end
if isempty([probe.al])
    menu('No anchor points positions. Use Make Probe in the Tools menu','OK');
    return;
end
if isempty([probe.sl])
    menu('No springs list. Use Make Probe in the Tools menu','OK');
    return;
end

% Get registered optode positions and then display springs 
probe = registerProbe2Head(probe,headvol,refpts);




% --------------------------------------------------------------------
function [probe, fwmodel, labelssurf] = ...
    clearRegistration(probe, fwmodel, labelssurf, dirname)

probe = resetProbeGui(probe);
fwmodel = resetSensitivity(fwmodel,probe,dirname);
labelssurf  = resetLabelssurf(labelssurf);





% --------------------------------------------------------------------
function pushbuttonRegisterProbeToSurface_Callback(hObject, eventdata, handles)
global atlasViewer

refpts       = atlasViewer.refpts;
probe        = atlasViewer.probe;
headsurf     = atlasViewer.headsurf;
headvol      = atlasViewer.headvol;
dirnameSubj  = atlasViewer.dirnameSubj;
fwmodel      = atlasViewer.fwmodel;
imgrecon     = atlasViewer.imgrecon;
labelssurf   = atlasViewer.labelssurf;

if isempty(probe.optpos)
    menu('No probe has been loaded or created. Use the SDgui to make or load a probe','ok');
    probe = resetProbe(probe);
    return;
end

% Finish registration
if isempty(probe.sl)
    
    % Register probe by simply pulling (or pushing) optodes toward surface
    % toward (or away from) center of head.
    method = 'digpts';
    probe = pullProbeToHeadsurf(probe,headvol);
    probe.hOptodesIdx = 1;
   
else
    
    % Register probe using springs based method
    method = 'springs';
    probe = probeRegisterSpringsMethod(probe,headvol,refpts);
  
end

% Clear old registration from gui after registering probe to avoid 
% lag time between diplay of initial probe and registered probe
[probe, fwmodel, labelssurf] = ...
    clearRegistration(probe, fwmodel, labelssurf, dirnameSubj);

% View registered optodes on the head surface
probe = viewProbe(probe, 'registered');

% Draw measurement list and save handle
probe = findMeasMidPts(probe);

fwmodel = updateGuiControls_AfterProbeRegistration(probe, fwmodel, imgrecon, labelssurf, dirnameSubj);

probe.hOptodesIdx = 1; 
probe = setProbeDisplay(probe, headsurf, method);

atlasViewer.probe       = probe;
atlasViewer.fwmodel = fwmodel;
atlasViewer.labelssurf  = labelssurf;




% --------------------------------------------------------------------
function menuItemMakeProbe_Callback(hObject, eventdata, handles)
global atlasViewer

probe        = atlasViewer.probe;
labelssurf   = atlasViewer.labelssurf;

hSDgui = atlasViewer.probe.handles.hSDgui;
if isempty(which('SDgui'))
    menu('SDgui doesn''t exist in the search path.','OK');
    return;
end
if ishandles(hSDgui)
    menu('SDgui already active.','OK');
    return;
end
atlasViewer.probe = resetProbe(atlasViewer.probe);
atlasViewer.probe.handles.hSDgui = SDgui(atlasViewer.dirnameProbe,'userargs');
set(atlasViewer.probe.handles.pushbuttonRegisterProbeToSurface,'enable','on');

% Clear labels faces associated with probe to cortex projection (we mark 
% the faces red). It's all new for a new probe.
labelssurf = resetLabelssurf(labelssurf);




% --------------------------------------------------------------------
function menuItemExit_Callback(hObject, eventdata, handles)
global atlasViewer
probe = atlasViewer.probe;

if ishandles(probe.handles.hSDgui) 
    delete(probe.handles.hSDgui);
    probe.handles.hSDgui=[];
end
delete(atlasViewer.handles.figure);
atlasViewer=[];




% --------------------------------------------------------------------
function checkboxHideProbe_Callback(hObject, eventdata, handles)
global atlasViewer;
probe    = atlasViewer.probe;
headsurf = atlasViewer.headsurf;

hideProbe = get(hObject,'value');
probe.hideProbe = hideProbe;
sl = probe.sl;

if isempty(sl)
    probe = setProbeDisplay(probe, headsurf, 'digpts');
else
    probe = setProbeDisplay(probe, headsurf, 'springs');
end

atlasViewer.probe = probe;



% --------------------------------------------------------------------
function checkboxHideMeasList_Callback(hObject, eventdata, handles)
global atlasViewer;
probe = atlasViewer.probe;
headsurf = atlasViewer.headsurf;

hideMeasList = get(hObject,'value');
probe.hideMeasList = hideMeasList;
probe = setProbeDisplay(probe, headsurf);

atlasViewer.probe = probe;


% --------------------------------------------------------------------
function checkboxHideSprings_Callback(hObject, eventdata, handles)
global atlasViewer;
probe = atlasViewer.probe;
headsurf = atlasViewer.headsurf;

hideSprings = get(hObject,'value');
probe.hideSprings = hideSprings;

if hideSprings==0
    set(handles.editSpringLenThresh,'visible','on');
    set(handles.textSpringLenThresh,'visible','on');
else
    set(handles.editSpringLenThresh,'visible','off');
    set(handles.textSpringLenThresh,'visible','off');
end

probe = setProbeDisplay(probe,headsurf,'springs');

atlasViewer.probe = probe;


% --------------------------------------------------------------------
function checkboxHideDummyOpts_Callback(hObject, eventdata, handles)
global atlasViewer;
probe = atlasViewer.probe;
headsurf = atlasViewer.headsurf;

hideDummyOpts = get(hObject,'value');
probe.hideDummyOpts = hideDummyOpts;
probe = setProbeDisplay(probe,headsurf,'springs');

atlasViewer.probe = probe;



% --------------------------------------------------------------------
function dirname = getSubjDir(arg)

dirname = -1;
if length(arg) > 1
    dirname = arg{1};
else
    
    % Rules fr determining if current folder is a subject folder    
    % Check for presence of atlasviewer or homer2 files 
    % in the current folder
    
    % 1. Check for presense of ./anatomical/headsurf.mesh
    dirname = [];
    if exist([pwd, '/anatomical'], 'dir')
        files = dir([pwd, '/anatomical/headsurf.mesh']);
        if ~isempty(files)
            dirname = pwd;
        end
    end
    
    % 2. Check for presense of ./anatomical/headsurf.mesh
    if exist([pwd, '/fw'], 'dir')
        files = dir([pwd, '/fw/fw_all.*']);
        if ~isempty(files)
            dirname = pwd;
        end
        files = dir([pwd, '/fw/headvol.vox']);
        if ~isempty(files)
            dirname = pwd;
        end
    end
    
    % 3. Check for presense of digpts.txt
    if exist([pwd, '/digpts.txt'], 'file')
        dirname = pwd;
    end
    
    % 4. Check for presense of atlasViewer.mat
    if exist([pwd, '/atlasViewer.mat'], 'file')
        dirname = pwd;
    end
    
    % 5. Check for presense of groupResults.mat
    if exist([pwd, '/groupResults.mat'], 'file')
        dirname = pwd;
    end
        
    % 6. Check for presense of SD or .nirs files
    files = dir([pwd, '/*.SD']);
    if ~isempty(files)
        dirname = pwd;
    end
    
    % 7. Check for presense of SD or .nirs files
    files = dir([pwd, '/*.nirs']);
    if ~isempty(files)
        dirname = pwd;
    end
    
    % After checking all the above for insications of subject folder
    % see if dirname is etill empty. If it is ask user for subject dir. 
    if isempty(dirname)
        pause(.1);
        dirname = uigetdir(pwd, 'Please select subject folder');
        if dirname==0
            dirname = pwd;
        end
    end
    
end

if isempty(dirname) | dirname==0
    return;
end

cd(dirname);

dirname(dirname=='\') = '/';

if dirname(end) ~= '/'
    dirname(end+1) = '/';
end



% --------------------------------------------------------------------
function dirname = getAtlasDir(arg)

if ~exist('arg','var')
    arg={};
end

dirname = '';
    
% First check argument for existence of atlas dir
if length(arg) > 2
    dirname = arg{2};
end

% No argument supplied or argument supplied but directory doen't exist
% In this case try finding default using the search paths if we are in
% an IDE.
if isempty(dirname) | ~exist(dirname,'file')
    
    % If we weren't able to get atlas path from config file, then try finding 
    % default dir or have user select it.
    if isempty(dirname)
    
        % If we are a standalone executable using search paths won't work.
        % Force user to select atlas directory.
        dirname = ffpath('Colin/anatomical/headsurf.mesh');
        if isempty(dirname)
            dirname = getDesktopDir();
            if exist([dirname, '/Colin/anatomical/headsurf.mesh'],'file') | exist([dirname, '/Colin/anatomical/headvol.vox'],'file')
                dirname = [dirname, '/Colin'];
                return;
            end
            fprintf('Ask user for atlas dirname.\n');
            dirname = selectAtlasDir(dirname);
            
            % Otherwise success! (Found default atlas dir). Strip path of trailing
            % anatomical/headsurf.mesh to get at the directory.
        elseif ~isempty(dirname)
            dirname = [dirname, '/Colin'];
            fprintf('Found atlas dirname: %s\n', dirname);
        end
        
    end
    
% Argument supplied, dir exists but doesn't contain anatomical dir. It's an invitation
% for the user to pick the atlas from a list of atlases. Supposedly dirname contains
% a database of atlases.
elseif exist(dirname,'file') & ~exist([dirname filesep 'anatomical'],'file')

    if ~exist([dirname filesep 'Colin'],'file') | ~exist([dirname filesep 'Colin' filesep 'anatomical'],'file')        
        dirname = selectAtlasDir(dirname);

    % This is the default. If the root atlas directory exists and contain a valid Colin
    % atlas directory
    elseif exist([dirname filesep 'Colin'],'file')
        dirname = [dirname filesep 'Colin'];
    end
    
end

% Check if we still have no atlas dir and warn user if that's the case
if isempty(dirname) | dirname==0
    menu('Warning: Couldn''t find default atlas directory.','OK');
    dirname = '';
    return;
end

dirname(dirname=='\') = '/';

% Add trailing file separator to dirname if there is none
if dirname(end) ~= '/' 
    dirname(end+1) = '/';
end




% --------------------------------------------------------------------
function dirname = selectAtlasDir(dirname)

if ~exist('dirname','var') | ~exist(dirname,'dir')
    dirname=pwd;
end
if exist(dirname,'file') & exist([dirname, '/anatomical'],'dir')
    if exist([dirname, '/anatomical/headsurf.mesh'],'file') | exist([dirname, '/anatomical/headvol.vox'],'file')
        return;
    end
end
    
while 1
    dirname = uigetdir(dirname,'Atlas directory not found. Please choose atlas directory');
    if dirname==0
        break;
    end
    if exist(dirname,'dir') & (exist([dirname, '/anatomical'],'dir') | ...
       (exist([dirname, '/mri'],'dir') & exist([dirname, '/surf'],'dir')))
        break;
    else
        q = menu('Selected directory is not a valid atlas directory. Please select valid atlas directory or cancel to avoid selecting','Select','Cancel');
        if q==1
            continue;
        elseif q==2
            dirname = [];
            break;
        end
    end
end




% --------------------------------------------------------------------
function pushbuttonZoomIn_Callback(hObject, eventdata, handles)
global atlasViewer;
axesv = atlasViewer.axesv;
ax=[];
for ii=1:length(axesv)
    if ishandles(axesv(ii).handles.pushbuttonZoomIn)
        if hObject==axesv(ii).handles.pushbuttonZoomIn
            ax=axesv(ii);
            break;
        end
    end
end

camzoom(ax.handles.axesSurfDisplay, ax.zoomincr);
atlasViewer.axesv(ii) = ax;




% --------------------------------------------------------------------
function pushbuttonZoomOut_Callback(hObject, eventdata, handles)
global atlasViewer;
axesv = atlasViewer.axesv;
ax=[];
for ii=1:length(axesv)
    if ishandles(axesv(ii).handles.pushbuttonZoomOut)
        if hObject==axesv(ii).handles.pushbuttonZoomOut
            ax=axesv(ii);
            break;
        end
    end
end

camzoom(ax.handles.axesSurfDisplay, 1/ax.zoomincr);
atlasViewer.axesv(ii) = ax;




% --------------------------------------------------------------------
function menuItemChangeSubjDir_Callback(hObject, eventdata, handles)
global atlasViewer

dirnameAtlas = atlasViewer.dirnameAtlas;
groupSubjList = atlasViewer.groupSubjList;
axesv = atlasViewer.axesv;
fwmodel = atlasViewer.fwmodel;
hGroupList = atlasViewer.handles.hGroupList;

if ishandles(hGroupList)
    displayGroupSubjList(groupSubjList, hGroupList);
    return;
end  

dirnameSubj = uigetdir('*.*','Change current subject directory');
if dirnameSubj==0
    return;
end
if length(axesv)>1
    if ishandles(axesv(2).handles.axesSurfDisplay)
        hp = get(axesv(2).handles.axesSurfDisplay,'parent');
        delete(hp);
    end
end
AtlasViewerGUI(dirnameSubj, dirnameAtlas, fwmodel.mc_exepath, hGroupList, 'userargs');




% --------------------------------------------------------------------
function menuItemChangeAtlasDir_Callback(hObject, eventdata, handles)
global atlasViewer
dirnameSubj = atlasViewer.dirnameSubj;
dirnameAtlas = atlasViewer.dirnameAtlas;
axesv = atlasViewer.axesv;
fwmodel = atlasViewer.fwmodel;


% Get the directory containing all the atlases and pass it to the 
% function which lets the user select from the list of atlases.
if ~isempty(dirnameAtlas)
    k = find(dirnameAtlas=='\' | dirnameAtlas=='/');
    dirnameAtlas = selectAtlasDir(dirnameAtlas(1:k(end-1)));
else
    dirnameAtlas = getAtlasDir({''});
end
if isempty(dirnameAtlas) 
    return;
end
if length(dirnameAtlas)==1 
    if dirnameAtlas==0
        return;
    end
end
if ~exist([dirnameSubj,'viewer'],'dir')
    mkdir([dirnameSubj,'viewer']);
else
    delete([dirnameSubj 'viewer/headvol*']);
end
fid = fopen([dirnameSubj,'viewer/settings.cfg'],'w');
fprintf(fid,'%s',dirnameAtlas);
fclose(fid);

% Restart AtlasViewerGUI with the new atlas directory.
AtlasViewerGUI(dirnameSubj, dirnameAtlas, fwmodel.mc_exepath, 'userargs');




% --------------------------------------------------------------------
function menuItemSaveRegisteredProbe_Callback(hObject, eventdata, handles)
global atlasViewer

probe      = atlasViewer.probe;
refpts     = atlasViewer.refpts;

optpos_reg = probe.optpos_reg;
nsrc       = probe.nsrc;

q = menu('Saving registered probe in probe_reg.txt - is this OK? Choose ''No'' to save in other filename or format','Yes','No');
if q==2

    [filename pathname] = uiputfile({'*.mat';'*.txt'},'Save registered probe to file');
    if filename==0
        return;
    end

elseif q==1

    filename = 'probe_reg.txt';

end

qq = menu('Do you want to include the 10-20 reference points?','Yes','No');


k = find(filename=='.');
ext = filename(k(end)+1:end);
if strcmpi(ext,'txt') || isempty(ext)
    
    fid = fopen(filename,'w');
    optpos_s = optpos_reg(1:nsrc,:);
    optpos_d = optpos_reg(nsrc+1:end,:);
    for ii=1:size(optpos_s,1)
        fprintf(fid,'s%d: %0.15f %0.15f %0.15f\n',ii,optpos_s(ii,:));
    end
    for ii=1:size(optpos_d,1)
        fprintf(fid,'d%d: %0.15f %0.15f %0.15f\n',ii,optpos_d(ii,:));
    end
    
    if qq==1
        fprintf(fid,'\n\n\n');
        for ii=1:size(refpts.pos,1)
            fprintf(fid,'%s: %.1f %.1f %.1f\n',refpts.labels{ii},refpts.pos(ii,1),refpts.pos(ii,2),refpts.pos(ii,3) );
        end
    end
    
    fclose(fid);

elseif strcmpi(extenstion,'mat')

    save(filename,'-mat','optpos_reg','nsrc');

end





% --------------------------------------------------------------------
function probe = probe2atlasSpace(headsurf,probe,digpts,refpts)

if isempty(probe.optpos)
    return;
end

% Assign the final atlas and subj fields
atlas.head   = headsurf.mesh.vertices;
jj=1;
for ii=1:length(refpts.labels)
    k = find(strcmpi(digpts.refpts.labels, refpts.labels{ii}));
    if ~isempty(k)
        subj.p1020(jj,:)  = digpts.refpts.pos(k,:);
        subj.l1020{jj}    = digpts.refpts.labels{k};
        atlas.p1020(jj,:) = refpts.pos(ii,:);
        atlas.l1020{jj}   = refpts.labels{ii};
        jj=jj+1;
    end
end

% Bring optodes to atlas space
subj.optodes = probe.optpos;
if ~isfield(subj,'anchor')
    menu('Subject data has no anchor point(s). Defaulting to canonical registration','OK');
    method = 'canonical';
else
    method = getProbeRegMethod();
end
[probe.optpos, T] = reg_subj2atlas(method, subj, atlas);





% --------------------------------------------------------------------
function B = probeDigptsRelated(probe,digpts)

B = 0;

p = probe.optpos;
nopt = probe.noptorig;
nsrc = probe.nsrc;
ndet = nopt-nsrc;

d = [digpts.srcpos; digpts.detpos];
ndigsrc = size(digpts.srcpos,1);
ndigdet = size(digpts.detpos,1);
ndigpts = ndigsrc + ndigdet;

if ndigpts==0
    return;
end
if isempty(p)
    return;
end
if ndigsrc~=nsrc
    return;
end
if ndigdet~=ndet
    return;
end

B = 1;





% --------------------------------------------------------------------
function menuItemImportProbe_Callback(hObject, eventdata, handles)
global atlasViewer

dirnameProbe = atlasViewer.dirnameProbe;
dirnameSubj  = atlasViewer.dirnameSubj;
probe        = atlasViewer.probe;
refpts       = atlasViewer.refpts;
headsurf     = atlasViewer.headsurf;
labelssurf   = atlasViewer.labelssurf;
fwmodel      = atlasViewer.fwmodel;
imgrecon     = atlasViewer.imgrecon;
digpts       = atlasViewer.digpts;
axesv        = atlasViewer.axesv;

[filename pathname] = uigetfile([dirnameProbe '*.*'],'Import subject probe');
if filename==0
    return;
end

% New probe means resetting probe, anatomical labels and sensitivity profile
probe       = resetProbe(probe);
fwmodel     = resetFwmodel(fwmodel);
imgrecon    = resetImgRecon(imgrecon);
labelssurf  = resetLabelssurf(labelssurf);

% Get optodes from various file formats
k = find(filename=='.');
fext = filename(k+1:end);
switch lower(fext)
    case 'mat'

        load([pathname filename],'-mat');        
        if ~exist('subj','var')
            menu('Error: Wrong file format for digitized points file. No subj field','ok');
            return;
        end
        digpts.refpts.pos     = subj.p1020;
        digpts.refpts.labels  = subj.l1020;
        probe = probe2atlasSpace(headsurf, probe, digpts, refpts);

        probe.center = headsurf.center;
        probe.orientation = headsurf.orientation;
    case 'txt'
        
        probe = getProbe(probe, [pathname,filename], headsurf);
        digpts = getDigpts(digpts,[pathname,filename]); 
        if isempty(probe.optpos)
            probe.optpos = digpts.pcpos;
            probe.srcpos = digpts.pcpos;
            probe.nsrc = size(digpts.pcpos,1);
            probe.nopt = size(digpts.pcpos,1);
            probe.noptorig = size(digpts.pcpos,1);
            probe.center = digpts.center;
            probe.orientation = digpts.orientation;
        else 
            probe = probe2atlasSpace(headsurf, probe, digpts, refpts);
        end        
        probe.center = headsurf.center;
        probe.orientation = headsurf.orientation;
        
    case {'sd','nirs'}

        filedata = load([pathname filename],'-mat');
        probe = loadSD(probe,filedata.SD);
        
        if (isempty(digpts.srcpos) & isempty(digpts.detpos)) | ~probeDigptsRelated(probe,digpts)
            % Check if the probe is already registered. One indicator is if
            % the probe isn't flat. Another is if all optodes are on or close to the head surface. 
            [foo1, foo2, d] = nearest_point(headsurf.mesh.vertices, probe.optpos);
            if mean(d)>5 | std(d)>1
            	% We are loading a flat probe that needs to be anchored to the
            	% head. Bring flat probe to some point on head surface to make the
            	% imported probe at least somewhat visible. To do this find
            	% a reference point on the head surface (e.g. Cz) to anchor (ap)
            	% the center of the probe to and translate the probe to that
            	% anchor point.
            	k = find(strcmpi(refpts.labels,'Cz'));
            	if isempty(k)
	                ap = refpts.pos(1,:);
    	        else
        	        ap = refpts.pos(k,:);
            	end
            	c = findcenter(probe.optpos);
            	tx = ap(1)-c(1);
            	ty = ap(2)-c(2);
            	tz = ap(3)-c(3);
            	T = [1 0 0 tx; 0 1 0 ty; 0 0 1 tz; 0 0 0 1];
            	probe.optpos = xform_apply(probe.optpos, T);
            end
            digpts = disableDigpts(digpts);
            
            probe.center = headsurf.center;
            probe.orientation = headsurf.orientation;            
        elseif probeDigptsRelated(probe,digpts)
            % We have digitized some or all source and detector optodes.
            % These digitized points now serve as anchors for the
            % corresponding sources and detectors.
            % NOTE that we are not using digpts directly as in this case
            % digpts should be represented in atlasViewer.probe            
            nSrcDig = size(digpts.srcpos,1);
            nSrc = probe.nsrc;
            nDetDig = size(digpts.detpos,1);
            nDet = probe.noptorig - probe.nsrc;
            
            for iS = 1:nSrcDig
                probe.al{iS,1} = digpts.srcmap(iS);
                probe.al{iS,2} = num2str(digpts.srcpos(iS,:));
            end
            for iD = 1:nDetDig
                probe.al{nSrcDig+iD,1} = nSrc+digpts.detmap(iD);
                probe.al{nSrcDig+iD,2} = num2str(digpts.detpos(iD,:));
            end
            probe.center = digpts.center;
            probe.orientation = digpts.orientation;            
        end

    otherwise

        menu('Error: Unknown file format for digitized points file.','ok');
        return;
end

% Check if probe anchor points exist in refpts
if iscell(probe.al)
    nanchor = size(probe.al,1);
    for ii = 1:nanchor
        if sum(strcmpi(refpts.labels, probe.al{ii,2}))==0
            menu('The selected probe uses anchors points not included in the current anatomy.','Ok');
            return
        end
    end
else
    probe.al = {};
end
probe = viewProbe(probe,'unregistered');

% This is done to not display dummy points by default. It does nothing 
% if the method isn't spring registration.
probe = setProbeDisplay(probe,headsurf,'springs');

atlasViewer.probe        = probe;
atlasViewer.dirnameProbe = pathname;
atlasViewer.labelssurf   = labelssurf;
atlasViewer.digpts       = digpts;
atlasViewer.fwmodel      = fwmodel;
atlasViewer.imgrecon     = imgrecon;




% --------------------------------------------------------------------
function hray = drawRayProjection(p1,p2,headsurf)

if leftRightFlipped(headsurf)
    axesOrd=[2 1 3];
else
    axesOrd=[1 2 3];
end    

hray = line([p1(axesOrd(1)),p2(axesOrd(1))],[p1(axesOrd(2)),p2(axesOrd(2))],...
            [p1(axesOrd(3)),p2(axesOrd(3))],'color','m','linewidth',2);       
set(hray,'tag','MNI projection');
drawnow();




% --------------------------------------------------------------------
function menuItemChooseLabelsColormap_Callback(hObject, eventdata, handles)
global atlasViewer

hLabelsSurf     = atlasViewer.labelssurf.handles.surf;
vertices        = atlasViewer.labelssurf.mesh.vertices;
idxL            = atlasViewer.labelssurf.idxL;
namesL          = atlasViewer.labelssurf.names;
colormaps       = atlasViewer.labelssurf.colormaps;
colormapsIdx    = atlasViewer.labelssurf.colormapsIdx;
iFaces          = atlasViewer.labelssurf.iFaces;

if ~ishandles(hLabelsSurf)
    return;
end

n = length(colormaps);
cmLst = cell(n,1);
for ii=1:n
    cmLst{ii} = sprintf('%s',colormaps(ii).name);
end
cmLst{n+1} = 'Cancel';
ch = menu('Choose Labels Colormap',cmLst);
if ch>n
    return;
end
cm = colormaps(ch).col;
faceVertexCData = cm(idxL,:);
faceVertexCData(iFaces,:) = repmat([1 0 0],length(iFaces),1);
set(hLabelsSurf,'faceVertexCData',faceVertexCData);
atlasViewer.labelssurf.colormapsIdx = ch;




% --------------------------------------------------------------------
function menuItemRegisterAtlasToDigpts_Callback(hObject, eventdata, handles)
global atlasViewer
global DEBUG

refpts       = atlasViewer.refpts;
digpts       = atlasViewer.digpts;
headsurf     = atlasViewer.headsurf;
headvol      = atlasViewer.headvol;
pialsurf     = atlasViewer.pialsurf;
probe        = atlasViewer.probe;
labelssurf   = atlasViewer.labelssurf;
dirnameAtlas = atlasViewer.dirnameAtlas;
dirnameSubj  = atlasViewer.dirnameSubj;
axesv        = atlasViewer.axesv;
fwmodel      = atlasViewer.fwmodel;
imgrecon     = atlasViewer.imgrecon; 

% Check conditions which would make us exit early
if digpts.refpts.isempty(digpts.refpts)
    return;
end 
if refpts.isempty(refpts)
    return;
end 
if all(isregistered(refpts,digpts))
    return;
end

if ~exist('hObject','var')
    hObject = [];
end

% If call menuItemRegisterAtlasToDigpts_Callback is not a GUI event, 
% then exit after setting  setting 
if ~ishandles(hObject)
    return;
end

%%%% Move all the volumes, surfaces and points back to a known space 
%%%% volume space. 



%%%% Move all the volumes, surfaces and points to viewer space



% First determine transformation to monte carlo space from volume 

% Generate transformation from head volume to digitized points space
[rp_atlas, rp_subj] = findCorrespondingRefpts(refpts, digpts);

headvol.imgOrig = headvol.img;
headvol.T_2digpts = gen_xform_from_pts(rp_atlas, rp_subj);

% Register headvol to digpts but first check fwmodel if it's volume 
% is already registered to digpts. if it is then set the headvol object 
% to the fwmodel's headvol and reuse it. 
if ~isregisteredFwmodel(fwmodel, headvol)
    
    [headvol.img, digpts.T_2mc] = ...
        xform_apply_vol_smooth(headvol.img, headvol.T_2digpts);
    
    headvol.T_2mc   = digpts.T_2mc * headvol.T_2digpts;
    headvol.center = xform_apply(headvol.center, headvol.T_2mc);
    
    % The MC space volume changed invalidating fwmodel meshes and 
    % vol to surface mesh. We need to recalculate all of this.
    fwmodel  = resetFwmodel(fwmodel, headvol);
    imgrecon = resetImgRecon(imgrecon);
else
    
    % Reusing MC space headvol from fwmodel. 
    headvol = fwmodel.headvol;
    
    % We know that headvol.T_2mc = digpts.T_2mc * headvol.T_2digpts.
    % Here we need to recover digpts.T_2mc. We can do this from 
    % headvol.T_2mc and headvol.T_2digpts with a little matrix algebra
    digpts.T_2mc = headvol.T_2mc / headvol.T_2digpts;
    
end

% Move digitized pts to monte carlo space
digpts.refpts.pos = xform_apply(digpts.refpts.pos, digpts.T_2mc);
digpts.pcpos      = xform_apply(digpts.pcpos, digpts.T_2mc);
digpts.srcpos     = xform_apply(digpts.srcpos, digpts.T_2mc);
digpts.detpos     = xform_apply(digpts.detpos, digpts.T_2mc);
digpts.optpos     = [digpts.srcpos; digpts.detpos];
digpts.center     = digpts.refpts.center;

% Copy digitized optodes to probe object
probe.optpos = digpts.optpos;

% move head surface to monte carlo space 
headsurf.mesh.vertices   = xform_apply(headsurf.mesh.vertices, headvol.T_2mc);
headsurf.center          = xform_apply(headsurf.center, headvol.T_2mc);
headsurf.centerRotation  = xform_apply(headsurf.centerRotation, headvol.T_2mc);

% move pial surface to monte carlo space 
pialsurf.mesh.vertices   = xform_apply(pialsurf.mesh.vertices, headvol.T_2mc);

% move anatomical labels surface to monte carlo space 
labelssurf.mesh.vertices = xform_apply(labelssurf.mesh.vertices, headvol.T_2mc);

% move ref points to monte carlo space 
refpts.pos = xform_apply(refpts.pos, headvol.T_2mc);

% No need to move fwmodel mesh, we already have our pialsurf in mc space 
% simply assign it to fwmodel. 
if isemptyMesh(fwmodel.mesh)
    fwmodel.mesh       = pialsurf.mesh;
    fwmodel.mesh_scalp = headsurf.mesh;
end

% No need to move imgrecon mesh, we already have our pialsurf in mc space 
% simply assign it to imgrecon. 
if isemptyMesh(imgrecon.mesh)
    imgrecon.mesh      = pialsurf.mesh;
end

if DEBUG
    headvol.mesh.vertices = xform_apply(headvol.mesh.vertices, headvol.T_2mc);
end

%%%% Now display all axes objects 

% Bounding box of axes objects might change - allow the axes limits 
% to change with it dynamically.
set(axesv(1).handles.axesSurfDisplay,{'xlimmode','ylimmode','zlimmode'},{'auto','auto','auto'});
axes(axesv(1).handles.axesSurfDisplay);
cla

% Set the orientation for display puposes to the dig pts 
[refpts, headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon] = ...
    setOrientationDigpts(digpts, refpts, headvol, headsurf, pialsurf, labelssurf, probe, fwmodel, imgrecon);

headsurf       = displayHeadsurf(headsurf);
pialsurf       = displayPialsurf(pialsurf);
labelssurf     = displayLabelssurf(labelssurf);
refpts         = displayRefpts(refpts);
digpts         = displayDigpts(digpts);
probe          = displayProbe(probe, headsurf);

axesv(1) = displayAxesv(axesv(1), headsurf, initDigpts());

set(axesv(1).handles.axesSurfDisplay,{'xlimmode','ylimmode','zlimmode'},{'manual','manual','manual'});

if ishandles(labelssurf.handles.surf)
    set(labelssurf.handles.menuItemSelectLabelsColormap,'enable','on');
else
    set(labelssurf.handles.menuItemSelectLabelsColormap,'enable','off');
end

probe = updateProbeGuiControls(probe,headsurf,'digpts');

%%% Save new atlas coordinates in atlasViewer
atlasViewer.headsurf    = headsurf;
atlasViewer.pialsurf    = pialsurf;
atlasViewer.labelssurf  = labelssurf;
atlasViewer.refpts      = refpts;
atlasViewer.headvol     = headvol;
atlasViewer.digpts      = digpts;
atlasViewer.probe       = probe;
atlasViewer.axesv       = axesv;
atlasViewer.fwmodel     = fwmodel;
atlasViewer.imgrecon    = imgrecon;




% --------------------------------------------------------------------
function menuItemGenerateMCInput_Callback(hObject, eventdata, handles)
global atlasViewer

fwmodel       = atlasViewer.fwmodel;
dirnameSubj   = atlasViewer.dirnameSubj;
probe         = atlasViewer.probe;

qAdotExists = 0;

% Check if there's a sensitivity profile which already exists
if exist([dirnameSubj 'fw/Adot.mat'],'file')
    qAdotExists = menu('Do you want to use the existing sensitivity profile in Adot.mat','Yes','No');
    if qAdotExists == 1
        fwmodel = menuItemGenerateLoadSensitivityProfile_Callback(hObject, struct('EventName','Action'), handles);
        if ~isempty(fwmodel.Adot)
            enableDisableMCoutputGraphics(fwmodel, 'on');
        end
    else
        delete([dirnameSubj 'fw/Adot*.mat']);
    end
end
    
% If neither a sensitivity not fluence profiles exist, then generate MC input and output
if qAdotExists~=1

    %%% Figure out wether we can proceed with the next step in the workflow
    %%% check whether we have what we need (i.e., MC output) to enable the
    %%% sensitivity profile menu item
    fwmodel = existMCOutput(fwmodel,probe,dirnameSubj);    
    if all(fwmodel.errMCoutput==3)
        if ~isempty(probe.ml)
            msg1 = sprintf('MC input and output already exist for this probe.\n');
            msg2 = sprintf('Use ''Generate/Load Sensitivity Profile'' under the\n');
            msg3 = sprintf('Forward Model menu to generate the sensitivity profile');
            menu([msg1,msg2,msg3],'OK');
            enableDisableMCoutputGraphics(fwmodel, 'on');
        else
            msg1 = sprintf('MC input and output already exist for this probe, but file with measurement list\n');
            msg2 = sprintf('is missing. NOTE: The .nirs file from an experiment using this probe\n');
            msg3 = sprintf('should contain the measurement list. Copy this file to the subject directory');
            menu([msg1,msg2,msg3],'OK');
            enableDisableMCoutputGraphics(fwmodel, 'off');
        end
    else
        if ismember(-1,fwmodel.errMCoutput)
            q = menu(sprintf('MC input does not match current probe. Generate new input and output?'),...
                'Yes','No');
            if q==1
                fwmodel = genMCinput(fwmodel, probe, dirnameSubj);
                fwmodel = genMCoutput(fwmodel, probe, dirnameSubj);
            end
        elseif ismember(-2,fwmodel.errMCoutput)
            q = menu(sprintf('MC input doesn''t match current MC settings. Generate new input and output?'),...
                'Yes','No');
            if q==1
                fwmodel = genMCinput(fwmodel, probe, dirnameSubj);
                fwmodel = genMCoutput(fwmodel, probe, dirnameSubj);
            end
        elseif all(fwmodel.errMCoutput>1)
            q = menu(sprintf('MC input exists but newer than ouput. Generate new input and output?'),...
                'Yes','No');
            if q==1
                fwmodel = genMCinput(fwmodel, probe, dirnameSubj);
                fwmodel = genMCoutput(fwmodel, probe, dirnameSubj);
            end
        elseif all(fwmodel.errMCoutput>0)
            q = menu(sprintf('MC input exists but no output.'),...
                'Overwrite input and generate output',...
                'Generate new output only',...
                'Cancel');
            if q==1
                fwmodel = genMCinput(fwmodel, probe, dirnameSubj);
                fwmodel = genMCoutput(fwmodel, probe, dirnameSubj);
            elseif q==2
                fwmodel = genMCoutput(fwmodel, probe, dirnameSubj);
            end
        else
            fwmodel = genMCinput(fwmodel, probe, dirnameSubj);
            fwmodel = genMCoutput(fwmodel, probe, dirnameSubj);
        end
    end
end

atlasViewer.fwmodel = fwmodel;




% --------------------------------------------------------------------
function msg = makeMsgOutOfMem(type)

switch(type)
    case 'outofmem'
        msg = 'There''s not enough contiguous memory to generate a sensitivity profile. ';
        msg = [msg, 'This can happen on older 32-bit systems or 32-bit matlab. '];
        msg = [msg, 'AtlasViewerGUI will try to increase the memory. It will require '];
        msg = [msg, 'a restart of the PC for these changes to take effect.'];
    case 'outofmemlinux'
        msg = 'An error occured which terminated the generation of the sensitivity profile. ';
        msg = [msg, 'One reason for this might be that the system ran out of contiguous memory. '];
        msg = [msg, 'Since matlab''s memory command is only implemented for Windows '];
        msg = [msg, 'the user has to monitor the memory themselves to confirm this. '];
        msg = [msg, 'Some possible solutions might be to increase swap space, terminataing other '];
        msg = [msg, 'processes or adding more physical memory.'];
    case 'restart'
        msg = 'Successfully changed setting for the amount of memory allocated to applications ';
        msg = [msg, '(ran dos command ''bcdedit /set IncreaseUserVa 3072''). '];
        msg = [msg, 'Now close all windows, restart the PC and rerun AtlasViewerGUI. '];
        msg = [msg, 'If that still doesn''t work consider running on a 64-bit system and/or '];
        msg = [msg, 'with a 64-bit matlab.'];
    case 'accessdenied'
        msg = 'Attempt to increase contiguous memory setting failed. This may be due to ';
        msg = [msg, 'lack of administrator priviledges for this account. '];
        msg = [msg, 'Consider acquiring administrator priviledges then re-running or '];
        msg = [msg, 'running AtlasViewerGUI on a 64-bit system and/or with a 64-bit matlab.'];
end



% --------------------------------------------------------------------
function fwmodel = menuItemGenerateLoadSensitivityProfile_Callback(hObject, eventdata, handles)
global atlasViewer

fwmodel     = atlasViewer.fwmodel;
imgrecon    = atlasViewer.imgrecon;
probe       = atlasViewer.probe;
headvol     = atlasViewer.headvol;
pialsurf    = atlasViewer.pialsurf;
headsurf    = atlasViewer.headsurf;
dirnameSubj = atlasViewer.dirnameSubj;
axesv       = atlasViewer.axesv;

T_vol2mc    = atlasViewer.headvol.T_2mc;

try 
    if isempty(eventdata) | strcmp(eventdata.EventName,'Action')
        fwmodel = genSensitivityProfile(fwmodel,probe,headvol,pialsurf,headsurf,dirnameSubj);
    end        
catch ME
    if strcmp(ME.identifier,'MATLAB:nomem') | ...
       strcmp(ME.identifier,'MATLAB:pmaxsize')
        if ispc()
            waitfor(msgbox(makeMsgOutOfMem('outofmem')));
            status = dos('bcdedit /set IncreaseUserVa 3072');
            if status>0
                waitfor(msgbox(makeMsgOutOfMem('accessdenied')));
            else
                waitfor(msgbox(makeMsgOutOfMem('restart')));
            end
            closeallwaitbars();
        else
            waitfor(msgbox(makeMsgOutOfMem('outofmemlinux')));
        end
    else
        rethrow(ME);
    end
    return;
end

% Turn off image recon display
imgrecon = showImgReconDisplay(imgrecon, axesv, 'off', 'off', 'off','off');

% Set image popupmenu to sencitivity
set(imgrecon.handles.popupmenuImageDisplay,'value',1);


if 0 
    hf = figure; 
    hAxes = axes;
    fwmodel = displaySensitivity(fwmodel, pialsurf, [], probe, hAxes);
else
    fwmodel = displaySensitivity(fwmodel, pialsurf, [], probe);
    set(pialsurf.handles.radiobuttonShowPial, 'value',0);
    uipanelBrainDisplay_Callback(pialsurf.handles.radiobuttonShowPial, [], handles);    
end
if ~isempty(fwmodel.Adot)
    imgrecon = enableImgReconGen(imgrecon,'on');
    imgrecon.mesh = fwmodel.mesh;
else
    imgrecon = enableImgReconGen(imgrecon,'off');
end

atlasViewer.fwmodel = fwmodel;
atlasViewer.imgrecon = imgrecon;



% --------------------------------------------------------------------
function menuItemGenFluenceProfile_Callback(hObject, eventdata, handles)
global atlasViewer

fwmodel     = atlasViewer.fwmodel;
imgrecon    = atlasViewer.imgrecon;
probe       = atlasViewer.probe;
headvol     = atlasViewer.headvol;
pialsurf    = atlasViewer.pialsurf;
headsurf    = atlasViewer.headsurf;
dirnameSubj = atlasViewer.dirnameSubj;
axesv       = atlasViewer.axesv;

err = 0;

fwmodel.fluenceProfFnames = {};
probe_section = initProbe();
inp = inputdlg_errcheck({'Number of fluence profiles per file'},'Number of fluence profiles per file', ...
                        1, {num2str(fwmodel.nFluenceProfPerFile)});
if ~isempty(inp)
    fwmodel.nFluenceProfPerFile = str2num(inp{1});
end
sectionsize = fwmodel.nFluenceProfPerFile;

fwmodel.fluenceProf(1) = initFluenceProf();
fwmodel.fluenceProf(2) = initFluenceProf();
for ii=1:sectionsize:probe.noptorig;
    fwmodel.fluenceProf(1).index = fwmodel.fluenceProf(1).index+1;

    iFirst = ii;
    iF = fwmodel.fluenceProf(1).index;
    if sectionsize>probe.noptorig
        iLast = probe.noptorig;
    else
        iLast = mod((iF-1)*sectionsize+sectionsize, probe.noptorig);
        if iFirst>iLast
            iLast = probe.noptorig;
            fwmodel.fluenceProf(1).last = true;
        end
    end
    
    fprintf('iFirst=%d, iLast=%d\n', iFirst, iLast);
    
    % Get number of sources in this probe section
    nsrc = probe.nsrc-iFirst+1;
    if nsrc<0
        nsrc=0;
    elseif iLast <= probe.nsrc
        nsrc = iLast-iFirst+1;
    end
    
    % Get number of dets in this probe section
    if iFirst > probe.nsrc
        ndet = iFirst-iLast+1;
    elseif iFirst <= probe.nsrc & iLast <= probe.nsrc
        ndet = 0;
    elseif iFirst <= probe.nsrc & iLast > probe.nsrc
        ndet = iLast-probe.nsrc;
    end
    
    probe_section.optpos_reg = probe.optpos_reg(iFirst:iLast,:);
    probe_section.nsrc     = nsrc;
    probe_section.ndet     = ndet;
    probe_section.nopt     = nsrc+ndet;
    probe_section.noptorig = nsrc+ndet;
    
    fwmodel = existMCOutput(fwmodel,probe_section,dirnameSubj);
    if ~all(fwmodel.errMCoutput==3 | fwmodel.errMCoutput==2)
        fwmodel = genMCinput(fwmodel, probe_section, dirnameSubj);
        [fwmodel, err] = genMCoutput(fwmodel, probe_section, dirnameSubj, 'noninteractive');
        s = iLast-iFirst+1;
        iFirst = 1;
        iLast = s;
    else
        fprintf('MC output already exists for fluence profile %d\n', fwmodel.fluenceProf(1).index);
    end
    
    if err~=0
        return;
    end
    
    fprintf('Finished with MC output for fluence profile %d\n', fwmodel.fluenceProf(1).index);
    fwmodel = genFluenceProfile(fwmodel, probe, iFirst, iLast, headvol, pialsurf, dirnameSubj);
    fprintf('Completed fluence profile %d\n\n', fwmodel.fluenceProf(1).index);
    
end

atlasViewer.fwmodel = fwmodel;





% --------------------------------------------------------------------
function radiobuttonShowDigpts_Callback(hObject, eventdata, handles)
global atlasViewer

digpts = atlasViewer.digpts;
hRefpts = digpts.handles.hRefpts;
hPcpos = digpts.handles.hPcpos;

val = get(hObject,'value');
if val==1
    if ~isempty(digpts.refpts.pos)
        set(hRefpts(:,1),'visible','on');
        set(hRefpts(:,2),'visible','off');
    end
    if ~isempty(digpts.pcpos)
        set(hPcpos(:,1),'visible','on');
        set(hPcpos(:,2),'visible','off');
    end
else
    set(hRefpts,'visible','off');
    set(hPcpos,'visible','off');
end

    


% --------------------------------------------------------------------
function checkboxOptodeSDMode_Callback(hObject, eventdata, handles)
global atlasViewer
probe    = atlasViewer.probe;
headsurf = atlasViewer.headsurf;

probe = setOptodeNumbering(probe);

probe = setProbeDisplay(probe, headsurf, 'springs');

atlasViewer.probe = probe;





% --------------------------------------------------------------------
function menuItemFs2Viewer_Callback(hObject, eventdata, handles)
global atlasViewer

dirnameSubj = atlasViewer.dirnameSubj;
dirnameAtlas = atlasViewer.dirnameAtlas;

[headvol, status1] = fs2headvol(dirnameSubj);
if status1~=0
    return;
end
status2 = headvol2headsurf(dirnameSubj, headvol);
status3 = fs2pialsurf(dirnameSubj);

status = status1+status2+status3;

if status>0
    menu(sprintf('Error in converting FreeSurfer files to Viewer format.'),'OK');
    return;
end

% Reload subject with it's own, newly-generated anatomical files
AtlasViewerGUI(dirnameSubj, dirnameAtlas, 'userargs');



% --------------------------------------------------------------------
function menuItemCalculateRefpts_Callback(hObject, eventdata, handles)
global atlasViewer

refpts        = atlasViewer.refpts;
headvol       = atlasViewer.headvol;
axesv         = atlasViewer.axesv;
dirnameSubj   = atlasViewer.dirnameSubj;
digpts        = atlasViewer.digpts;

switch(get(hObject, 'tag'))
    case 'menuItem10_20'
        eeg_system = '10-20';
    case 'menuItem10_10'
        eeg_system = '10-10';
    case 'menuItem10_5'
        eeg_system = '10-5';
    case 'menuItem10_2_5'
        eeg_system = '10-2.5';
    case 'menuItem10_1'
        eeg_system = '10-1';
    otherwise
        eeg_system = '10-20';
end

[nz,iz,rpa,lpa,cz] = getLandmarks(refpts);
if isempty(nz) | isempty(iz) | isempty(lpa) | isempty(rpa) | isempty(cz)
    menu(sprintf('One or more landmarks are missing - unable to calculate %s pts.', eeg_system),'OK');
    set(hObject,'enable','off');
    return;
end

%if ~isempty(headvol.img) % change the order of these two so it does surface first
%    refpts = calcRefpts(refpts, headvol);
if isfield(atlasViewer,'headsurf')
    [refpts, ~, err]  = calcRefpts(refpts, atlasViewer.headvol, {}, eeg_system);
    % Need to transform refpts from surface space to volume space
    % also calcRefpts.m saves the refpts to a file, but it doesn't
    % update the refpts2vol.txt file
end

if ~err
    refpts = displayRefpts(refpts);    
    saveRefpts(refpts, dirnameSubj, headvol.T_2mc, 'overwrite');
    atlasViewer.refpts = refpts;
else
    menu('Failed to calculate ref pts because of missing landmarks','OK');
end



% --------------------------------------------------------------------
function menuItemEnableSensitivityMatrixVolume_Callback(hObject, eventdata, handles)
global atlasViewer
fwmodel = atlasViewer.fwmodel;

checked = get(hObject,'checked');
if strcmp(checked,'on');
    set(hObject,'checked','off');
    val=0;
elseif strcmp(checked,'off');
    set(hObject,'checked','on');
    val=1;
else
    return;
end
fwmodel = enableSensitivityMatrixVolume(fwmodel,val);
atlasViewer.fwmodel = fwmodel;




% --------------------------------------------------------------------
function pushbuttonCopyFigure_Callback(hObject, eventdata, handles)
global atlasViewer
axesv       = atlasViewer.axesv;
headsurf    = atlasViewer.headsurf;
fwmodel 	= atlasViewer.fwmodel;

cm = colormap;
clim = caxis;

hf = figure;

hAxes = copyobj(axesv(1).handles.axesSurfDisplay, hf);

% colormap is a propery of figure not axes. Since we don't want to 
% copy the whole figure which is the gui but only the axes, we need to 
% set the figure colormap after the copyobj. 
h = colormap(cm);
caxis(clim);

axesv(1).handles.axesSurfDisplay = hAxes;
camzoom(axesv(1).handles.axesSurfDisplay, 1.3*axesv(1).zoomincr);




% --------------------------------------------------------------------
function menuItemConvertEZStoDigPts_Callback(hObject, eventdata, handles)
global atlasViewer

[filenm,pathnm] = uigetfile('*.ezs','Choose .ezs file to convert');
if filenm==0
    return;
end

wd = cd();
cd(pathnm)
ezmapper2txt(filenm);
cd(wd)

atlasViewer.digpts = getDigpts(atlasViewer.digpts, atlasViewer.dirnameSubj);
atlasViewer.probe = getProbe(atlasViewer.probe, atlasViewer.dirnameSubj, atlasViewer.headsurf);
atlasViewer.probe = displayProbe(atlasViewer.probe);


% --------------------------------------------------------------------
function menuItemSaveViewerState_Callback(hObject, eventdata, handles)
global atlasViewer


axesv       = atlasViewer.axesv(1);
headvol     = atlasViewer.headvol;
headsurf    = atlasViewer.headsurf;
pialsurf    = atlasViewer.pialsurf;
labelssurf  = atlasViewer.labelssurf;
refpts      = atlasViewer.refpts;
digpts      = atlasViewer.digpts;
probe       = atlasViewer.probe;
fwmodel 	= atlasViewer.fwmodel;
imgrecon    = atlasViewer.imgrecon;

dirnameSubj = atlasViewer.dirnameSubj;

saveObjects('atlasViewer.mat', ...
    axesv, ...
    headvol, ...
    headsurf, ...
    pialsurf, ...
    labelssurf, ...
    refpts, ...
    labelssurf, ...
    digpts, ...
    probe, ...
    fwmodel, ...
    imgrecon ...
    );





% --------------------------------------------------------------------
function checkboxOptodeCircles_Callback(hObject, eventdata, handles)
global atlasViewer

probe    = atlasViewer.probe;
headsurf = atlasViewer.headsurf;

val = get(hObject,'value');
if val==1
    probe.optViewMode='circles';
else
    probe.optViewMode='numbers';
end
probe = setOptodeNumbering(probe);
probe = setProbeDisplay(probe, headsurf, 'springs');

atlasViewer.probe = probe;




% --------------------------------------------------------------------
function menuItemLighting_Callback(hObject, eventdata, handles)
global atlasViewer

axesv       = atlasViewer.axesv;

if strcmp(get(hObject,'checked'), 'on');
    set(hObject,'checked', 'off');
    val=0;
elseif strcmp(get(hObject,'checked'), 'off');
    set(hObject,'checked', 'on');
    val=1;
end

name = get(hObject,'label');
iLight = str2num(name(end));

if val==1
    set(axesv(1).handles.lighting(iLight),'visible','on');
else
    set(axesv(1).handles.lighting(iLight),'visible','off');
end



% --------------------------------------------------------------------
function pushbuttonSaveRefpts_Callback(hObject, eventdata, handles)
global atlasViewer;
refpts = atlasViewer.refpts;
headsurf = atlasViewer.headsurf;
digpts = atlasViewer.digpts;
dirnameSubj = atlasViewer.dirnameSubj;
T_vol2mc = atlasViewer.headvol.T_2mc;

command = get(hObject,'String');

if strcmp(command, 'Save Ref Pts')
    
    if all(refpts.handles.selected==-1)
        return;
    end
       
    % However the head mesh we found the ref points on was also transformed
    % by headvol.T_2ras to view left/right correctly. If it's not identity 
    % we need to unapply it to the ref points - which is done inside
    % saveRefpts, when we pass it a viewing tranformation in the volume.    
    if AnatomyExists(dirnameSubj)
        mode = 'overwrite';
    else
        mode = 'nosave';
    end
    refpts = saveRefpts(refpts, dirnameSubj, T_vol2mc, mode);
    refpts = getRefpts(refpts, dirnameSubj);
    refpts = displayRefpts(refpts, atlasViewer.axesv(1).handles.axesSurfDisplay );    
    if length(refpts.labels)>=5
        set(refpts.handles.menuItemCalculateRefpts,'enable','on');
        if ~isempty(digpts.refpts.pos)
            set(digpts.handles.menuItemRegisterAtlasToDigpts,'enable','on');
        end
    else
        set(refpts.handles.menuItemCalculateRefpts,'enable','off');
    end
    
else
    
    labels{1} = command(6:end);
    [refpts, currentPt] = updateRefpts(refpts, headsurf, labels, atlasViewer.axesv(2).handles.axesSurfDisplay);
    disp(sprintf('%s: %s', labels{1}, num2str(currentPt)));
    
end
atlasViewer.refpts = refpts;




% --------------------------------------------------------------------
function findRefptsGUI_DeleteFcn(hObject, eventdata, handles)
global atlasViewer
atlasViewer.axesv(2).handles.axesSurfDisplay=-1;
atlasViewer.refpts.handles.selected(:)=-1;



% --------------------------------------------------------------------
function menuItemFindRefpts_Callback(hObject, eventdata, handles)
global atlasViewer

axesv = atlasViewer.axesv;
headsurf = atlasViewer.headsurf;
headvol = atlasViewer.headvol;

if length(axesv)==2
    if ishandles(axesv(2).handles.axesSurfDisplay)
        return;
    end
end

hf = figure('name','Find Ref Points','units','normalized','position',[.2,.1,.5,.7], 'menubar','none',...
            'toolbar','figure','color',[.2,.2,.2], 'deletefcn',@findRefptsGUI_DeleteFcn);
handles2.axesSurfDisplay = axes;
set(handles2.axesSurfDisplay, {'xlimmode','ylimmode','zlimmode'}, {'manual','manual','manual'});

set(handles2.axesSurfDisplay,'xlim')
set(handles2.axesSurfDisplay,'units','normalized','position',[.10 .30 .80 .65 ]);
handles2.panelRotateZoomAxes = uipanel(hf,'units','normalized','position',[.05,.01,.80,.20]);
handles2.pushbuttonZoomIn = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Zoom In',...
              'units','normalized', 'tag','pushbuttonZoomIn', 'position',[.70,.60,.10,.15], ...
              'callback',@pushbuttonZoomIn_Callback);
handles2.pushbuttonZoomOut = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Zoom Out',...
              'units','normalized', 'tag','pushbuttonZoomOut', 'position',[.70,.40,.10,.15], ...
              'callback',@pushbuttonZoomOut_Callback);
handles2.pushbuttonSaveRefpts = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Save Ref Pts',...
              'units','normalized', 'tag','pushbuttonSaveRefpts', 'position',[.70,.10,.10,.15], ...
              'callback',@pushbuttonSaveRefpts_Callback);
handles2.pushbuttonSaveNz = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Save Nz',...
              'units','normalized', 'tag','pushbutton', 'position',[.85,.82,.10,.15], ...
              'callback',@pushbuttonSaveRefpts_Callback);
handles2.pushbuttonSaveIz = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Save Iz',...
              'units','normalized', 'tag','pushbutton', 'position',[.85,.62,.10,.15], ...
              'callback',@pushbuttonSaveRefpts_Callback);
handles2.pushbuttonSaveLPA = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Save LPA',...
              'units','normalized', 'tag','pushbutton', 'position',[.85,.42,.10,.15], ...
              'callback',@pushbuttonSaveRefpts_Callback);
handles2.pushbuttonSaveRPA = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Save RPA',...
              'units','normalized', 'tag','pushbutton', 'position',[.85,.22,.10,.15], ...
              'callback',@pushbuttonSaveRefpts_Callback);
handles2.pushbuttonSaveCz = ...
    uicontrol('parent',handles2.panelRotateZoomAxes, 'style','pushbutton', 'string','Save Cz',...
              'units','normalized', 'tag','pushbutton', 'position',[.85,.02,.10,.15], ...
              'callback',@pushbuttonSaveRefpts_Callback);
handles2.menuItemLight1=[];
handles2.menuItemLight2=[];
handles2.menuItemLight3=[];
handles2.menuItemLight4=[];
handles2.menuItemLight5=[];
handles2.menuItemLight6=[];
handles2.menuItemLight7=[];
handles2.menuItemLight8=[];
handles2.menuItemViewOrigin = handles.menuItemViewOrigin;

axesv(2) = initAxesv(handles2,1);
headsurf.handles.surf=[];
headsurf = displayHeadsurf(headsurf, handles2.axesSurfDisplay);
set(headsurf.handles.surf,'facealpha',.8);
axesv(2) = displayAxesv(axesv(2), headsurf, initDigpts());
if ~verLessThan('matlab','8.3')
    set(headsurf.handles.surf, 'pickableparts','visible','facealpha',1.0);
end

atlasViewer.axesv = axesv;

hd = datacursormode(hf);
set(hd,'UpdateFcn',@headsurfUpdateFcn,'SnapToDataVertex','on');
datacursormode on



% ---------------------------------------------------------------------
function txt = headsurfUpdateFcn(o,e)
global atlasViewer;

% Instead of getting position with "p = get(e,'position')", we get it 
% from with "p = e.Position". The former isn't compatible with matlab
% version beyond 2014a. 
p = e.Position;
atlasViewer.headsurf.currentPt = p;
txt = sprintf('%0.1f, %0.1f, %0.1f', p(1), p(2), p(3));



% --------------------------------------------------------------------
function menuProbePlacementVariation_Callback(hObject, eventdata, handles)

plotProbePlacementVariation();



% --------------------------------------------------------------------
function menuItemRegisterAtlasToHeadSize_Callback(hObject, eventdata, handles)
global atlasViewer

prompt = {'Head Circumference (cm):','Iz to Nz (cm):','RPA to LPA (cm):'};
dlg_title = 'Input Head Size';
num_lines = 1;
def = {'50','31','37'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

if isempty(answer)
    return
end

headSize.HC = str2num(answer{1})*10;
headSize.IzNz = str2num(answer{2})*10;
headSize.LPARPA = str2num(answer{3})*10;

xo = [70 70 70];
x = fminsearch( @(x) ellipse_1020_costfun(x,headSize),xo);

a = x(1); %  LPA to RPA axis
b = x(2); % Nz to Iz axis
c = x(3); % Cz axis

r10p = 18*3.14159/180;

Cz = [0 0 c];
RPA = [a*cos(r10p) 0 -c*sin(r10p)];
LPA = [-a*cos(r10p) 0 -c*sin(r10p)];
Nz = [0 b*cos(r10p) -c*sin(r10p)];
Iz = [0 -b*cos(r10p) -c*sin(r10p)];

atlasViewer.digpts.refpts.pos    = [Nz; Iz; RPA; LPA; Cz];
atlasViewer.digpts.refpts.labels = {'nz', 'iz', 'rpa', 'lpa', 'cz'};

% fid = fopen('digpts.txt','wt');
% fprintf(fid,'nz: %.2f %.2f %.2f\n',Nz);
% fprintf(fid,'iz: %.2f %.2f %.2f\n',Iz);
% fprintf(fid,'rpa: %.2f %.2f %.2f\n',RPA);
% fprintf(fid,'lpa: %.2f %.2f %.2f\n',LPA);
% fprintf(fid,'cz: %.2f %.2f %.2f\n',Cz);
% fclose(fid);

dirnameSubj  = atlasViewer.dirnameSubj;
if exist([dirnameSubj 'viewer/headvol.vox'],'file')
    delete([dirnameSubj 'viewer/headvol.vox']);
end

menuItemRegisterAtlasToDigpts_Callback(hObject, eventdata, handles)



% --------------------------------------------------------------------
function editSpringLenThresh_Callback(hObject, eventdata, handles)
global atlasViewer;

probe = atlasViewer.probe;
sl         = probe.sl;
hSprings = probe.handles.hSprings;

if ~isempty(probe.optpos_reg)
    optpos = probe.optpos_reg;
elseif ~isempty(probe.optpos)
    optpos     = probe.optpos;
else
    return;
end

cm = [0 0 1; 0 1 1; 0 0 0; 1 1 0; 1 0 0];
sLenThresh = probe.springLenThresh;

foo = str2num(get(hObject,'string'));
if length(foo)~=2
    set(hObject,'string',num2str(sLenThresh));
    return;
elseif foo(1)>=foo(2)
    set(hObject,'string',num2str(sLenThresh));
    return;
end
probe.springLenThresh = foo;
sLenThresh = probe.springLenThresh;

for ii=1:size(sl,1) 
    springLenReg(ii) = dist3(optpos(sl(ii,1),:), optpos(sl(ii,2),:));
    springLenErr(ii) = springLenReg(ii)-sl(ii,3);
    if springLenErr(ii)<-sLenThresh(2)
        k = 1;
    elseif springLenErr(ii)<-sLenThresh(1)
        k = 2;
    elseif springLenErr(ii)>sLenThresh(2)
        k = 5;
    elseif springLenErr(ii)>sLenThresh(1)
        k = 4;
    else
        k = 3;
    end
    set(hSprings(ii),'color',cm(k,:));
end

springLenErr




% --------------------------------------------------------------------
function menuCalcOptProp_Callback(hObject, eventdata, handles)

prompt        = {'Wavelengths (nm)','HbT (uM)','SO2 (%)','Reduced Scattering Coefficient at 800 nm (1/mm)','Scattering Slope b'};
name          = 'Calculate absorption and scattering';
numlines      = 1;
defaultanswer = {'690 830','60','65','0.8','1.5'};
   
answer = inputdlg(prompt,name,numlines,defaultanswer,'on');
if isempty(answer)
    return;
end

wv = str2num(answer{1});
hbt = str2num(answer{2}) * 1e-6;
so2 = str2num(answer{3}) / 100;
a = str2num(answer{4});
b = str2num(answer{5});

e = GetExtinctions( wv );

mua = e(:,1)*(hbt*so2) + e(:,2)*(hbt*(1-so2));
musp = a * exp(-b*(wv-800)/800);

ch = menu( sprintf('For wavelengths %s nm:\nmua = %s 1/mm\nmusp = %s 1/mm\n', num2str(wv),num2str(mua',4),num2str(musp,4)), 'okay');




% --------------------------------------------------------------------
function menuItemProjectProbeToCortex_Callback(hObject, eventdata, handles)
global atlasViewer

if isempty(atlasViewer.labelssurf.colormaps)
    menu( 'No cortical anatomical labels provided for this anatomy.','Ok');
    return;
end

% Get inputs
probe              = atlasViewer.probe;
probe              = initProbeProjection(probe);
optpos_reg         = probe.optpos_reg;
hOptodes           = probe.handles.hOptodes;
hMeasCortex        = probe.handles.hMeasCortex;
hMeasToLabelsProjTbl  = probe.handles.hMeasToLabelsProjTbl;
hRays              = probe.handles.hRays;
nopt               = probe.noptorig;
ml                 = probe.ml;

attractPt      = atlasViewer.headvol.center;

labelssurf     = atlasViewer.labelssurf;
labelssurf     = initLabelssurfProbeProjection(labelssurf);
hLabelsSurf    = labelssurf.handles.surf;
T_labelssurf2vol = labelssurf.T_2vol;
mesh           = labelssurf.mesh;
vertices       = labelssurf.mesh.vertices;
idxL           = labelssurf.idxL;
namesL         = labelssurf.names;

headsurf       = atlasViewer.headsurf;

T_headvol2mc   = atlasViewer.headvol.T_2mc;

hAtlasViewerGUI = atlasViewer.handles.figure;

tag = get(hObject,'tag');

% Project optodes to labeled cortex
if strcmp(tag, 'menuItemProjectOptodesToCortex') | isempty(ml)
    ptsProj = optpos_reg(1:nopt,:);

    % If projecting optodes rather than meas channels, display optodes 
    % in their original registered positions rather than the 
    % ones which were lifted off the head surface for easier viewing.
    probe.hOptodesIdx = 2;
    probe = setProbeDisplay(probe,headsurf);
    figname = 'Optode Projection to Cortex Labels';
else
    ptsProj = probe.mlmp;
    figname = 'Channel Projection to Cortex Labels';
end

% ptsProj_cortex is in viewer space. To get back to MNI coordinates take the
% inverse of the tranformation from mni to viewer space.
ptsProj_cortex = ProjectionBI(ptsProj, vertices);
[ptsClosest iP] = nearest_point(vertices, ptsProj_cortex);
ptsProj_cortex_mni = xform_apply(ptsProj_cortex,inv(T_headvol2mc*T_labelssurf2vol));

% Display optodes on labeled cortex
hMeasCortex = [];
hMeasToLabelsProjTbl = [];
pts = prepPtsStructForViewing(ptsProj_cortex, size(ptsProj_cortex,1), 'numbers','k',11);
hMeasCortex = viewPts(pts, attractPt,  0);
set(hMeasCortex,'visible','off');

% Generate rays from optodes on head to optodes on cortex and
% highlight the faces on the label that they pass through.
faceVertexCData = get(hLabelsSurf,'faceVertexCData');
faceVertexAlphaData = get(hLabelsSurf,'faceVertexAlphaData');
iFaces = [];
for ii=1:size(ptsProj,1)
    p1 = ptsProj(ii,:);
    p2 = ptsProj_cortex(ii,:);
    hRays(ii) = drawRayProjection(p1, p2, headsurf);
    v=p1-p2;
    [t,u,v,iFace] = raytrace(p1,v, mesh.vertices, mesh.faces);
    
    % Find closest face
    [foo,iFaceMin] = min(abs(t(iFace)));
    faceVertexCData(iFace(iFaceMin),:) = repmat([1 0 0],length(iFace(iFaceMin)),1);
    set(hLabelsSurf,'FaceVertexCData',faceVertexCData);
    faceVertexAlphaData(iFace(iFaceMin)) = ones(length(iFace(iFaceMin)),1);
    set(hLabelsSurf,'FaceVertexAlphaData',faceVertexAlphaData);
    iFaces = [iFaces iFace(iFaceMin)];
end
if all(ishandles(hRays))
    set(hRays,'color','k');
else
    return;
end

% Create table associating projected cortex optodes with brain labels
hMeasToLabelsProjTbl = figure('name',figname,'toolbar','none',...
                              'menubar','none','numbertitle','off', ...
                              'units','normalized', 'position',[.1,.02,.35,.8]);
                          
if strcmp(tag, 'menuItemProjectOptodesToCortex') | isempty(ml)
    optlabelsTbl = repmat({'','',''},length(iP),1);
    columnname = {'opt #','opt coord (Monte Carlo)','opt coord (MNI)','label name'};
    columnwidth = {40,160,120,140};
    for k=1:nopt
        optlabelsTbl{k,1} = num2str(k);
        optlabelsTbl{k,2} = num2str(round(ptsProj_cortex(k,:)));
        optlabelsTbl{k,3} = num2str(round(ptsProj_cortex_mni(k,:)));
        j = find(mesh.faces(:,1)==iP(k) | mesh.faces(:,2)==iP(k) | mesh.faces(:,3)==iP(k));
        optlabelsTbl{k,4} = namesL{idxL(j(1))};
        faceVertexCData(j,:) = repmat([1 0 0],length(j),1);
    end
else
    optlabelsTbl = repmat({'','','',''},length(iP),1);
    columnname = {'src','det','ch coord (Monte Carlo) ','ch coord (MNI)','label name'};
    columnwidth = {30,30,160,120,140};
    for k=1:size(ptsProj,1)
        optlabelsTbl{k,1} = num2str(ml(k,1));
        optlabelsTbl{k,2} = num2str(ml(k,2));
        optlabelsTbl{k,3} = num2str(round(ptsProj_cortex(k,:)));
        optlabelsTbl{k,4} = num2str(round(ptsProj_cortex_mni(k,:)));
        j = find(mesh.faces(:,1)==iP(k) | mesh.faces(:,2)==iP(k) | mesh.faces(:,3)==iP(k));
        optlabelsTbl{k,5} = namesL{idxL(j(1))};
        faceVertexCData(j,:) = repmat([1 0 0],length(j),1);
    end
end
ht = uitable('parent',hMeasToLabelsProjTbl,'columnname',columnname,...
             'units','normalized','position',[.1 .1 .8 .8],'columnwidth',columnwidth,...
             'data',optlabelsTbl);

% Save outputs
atlasViewer.probe.handles.hMeasCortex = hMeasCortex;
atlasViewer.probe.handles.hMeasToLabelsProjTbl  = hMeasToLabelsProjTbl;
atlasViewer.probe.handles.hRays = hRays;
atlasViewer.probe.ptsProj_cortex = ptsProj_cortex;
atlasViewer.probe.ptsProj_cortex_mni = ptsProj_cortex_mni;
atlasViewer.labelssurf.iFaces = iFaces;



% --------------------------------------------------------------------
function menuItemViewOrigin_Callback(hObject, eventdata, handles)
global atlasViewer

axesv = atlasViewer.axesv; 
hAxes = axesv.handles.axesSurfDisplay;

onoff = get(hObject,'checked');
if strcmp(onoff, 'on')
    set(hObject, 'checked','off');
elseif strcmp(onoff, 'off')
    set(hObject, 'checked','on');
end

viewOrigin(hAxes, 'donotredraw');




% --------------------------------------------------------------------
function menuItemGetSensitivityatMNICoordinates_Callback(hObject, eventdata, handles)
global atlasViewer
fwmodel    = atlasViewer.fwmodel;

% get user input
prompt={'MNI Coordinates (x, y, z) one or more separated by semicolumns ','radius (mm)','absorption change (mm-1)'};
name='Get brain activation sensitivity at MNI coordinates';
numlines=1;
defaultanswer={'[0 0 0; 0 0 0]','3','0.001'};

answer=inputdlg(prompt,name,numlines,defaultanswer,'on');

if isempty(answer)
    return
end

coordinate_mni = str2num(answer{1});
radius = str2num(answer{2});
absorption_change = str2num(answer{3});

no_mni = size(coordinate_mni,1);

% Load the sensitivity volume
if  exist([atlasViewer.dirnameSubj  'fw' filesep 'AdotVolSum.3pt'],'file')
    cd([atlasViewer.dirnameSubj '/fw'])
    fid = fopen('AdotVolSum.3pt','rb');
    v = single(fread(fid, 'float32'));
    cd(atlasViewer.dirnameSubj);
else
    msg = sprintf('You need to first generate a sensitivity profile with sensitivity matrix volume option enabled.');
    menu(msg,'OK');
    return;
end

% get the dimensions of the volume
dims = ones(1,4);
dims(1) = size(atlasViewer.headvol.img,1);
dims(2) = size(atlasViewer.headvol.img,2);
dims(3) = size(atlasViewer.headvol.img,3);

% reshape to get volume in 3D
f = reshape(v, dims);

% go from Monte Carlo space to MNI space
T_headvol2mc         = atlasViewer.headvol.T_2mc; % colin to Monte Carlo space
T_labelssurf2vol = atlasViewer.labelssurf.T_2vol; % mni to colin

h = waitbar(0,'Please wait while loading volume and calculating...');

% loop around MNI coordinates to get dOD change for each.
for i = 1: no_mni;
    
    coordinate_MC(i,:) = xform_apply(coordinate_mni(i,:), (T_headvol2mc*T_labelssurf2vol));
    
    % generate a sphere with MNI at the center
    vol = gen_blob([radius radius radius],coordinate_MC(i,:),ones(size(f)),0);
    
    % sum the sensitivity within the sphere (this sensitivity is the sum of
    % all optode sensitivities.)
    blob_ind = find(vol>1);
    sum_sensitivity = sum(f(blob_ind));
    vox_vol = 1; %mm^3   Talk to Jay to get the vox vol from the data structure as it is not necessarily 1 mm^3
    
    % estimate optical density change
    deltaMa = absorption_change * vox_vol;
    dOD_blob(i) = sum_sensitivity * deltaMa;
    
    clear vol blob_ind sum_sensitivity;
    waitbar(i/no_mni);
    
end
close(h);% close waitbar
% Save here for later plotting the projection on Atlas
atlasViewer.fwmodel.MNI = coordinate_mni;
atlasViewer.fwmodel.MNI_inMCspace = coordinate_MC;


% Create a table associating projected cortex optodes with brain labels
figname = 'Integrated Sensitivity of a Sphere with MNI Coordinates at the Center';
MNI_2_sensitivity_Table = figure('name',figname,'toolbar','none',...
    'menubar','none','numbertitle','off', ...
    'units','normalized', 'position',[.4 .1 .2 .4]);

cnames = {'x','y','z','radius', 'deltaMa', 'dOD'};
rnames = {'1','2','3'};
columnwidth = {40, 40, 40, 40, 80, 80};
columnformat = {'char','char','char','char','char','char'}; % to allign the entries to left

d = zeros(no_mni,6);
d(1:no_mni,1:3) = coordinate_mni;
d(1:no_mni,4) = radius;
d(1:no_mni,5) = deltaMa;
d(1:no_mni,6) = dOD_blob;

t = uitable(MNI_2_sensitivity_Table ,'Data',d,...
    'ColumnName',cnames,...
    'RowName',rnames,'position',[5 1.5 335 400],'columnwidth',columnwidth,'ColumnFormat', columnformat);
%[left bottom width height]

% Enable display_mni_projection
set(handles.checkbox_Display_MNI_Projection, 'enable','on');
set(handles.checkbox_Display_MNI_Projection, 'value',0);



% --------------------------------------------------------------------
function checkbox_Display_MNI_Projection_Callback(hObject, eventdata, handles)
global atlasViewer

fwmodel    = atlasViewer.fwmodel;
headvol    = atlasViewer.headvol; 
labelssurf = atlasViewer.labelssurf; 
headvol    = atlasViewer.headvol; 
headsurf   = atlasViewer.headsurf;
 
if get(hObject,'value')==1 % if checkbox is checked
    
    prompt={'MNI Coordinates (x, y, z) one or more separated by semicolumns'};
    name='Input MNI Coordinate(s)';
    numlines=1;
    
    % if user already ran brain activation MNI to get the sensitivy
    % have those MNIs as default
    if ~isfield(atlasViewer.fwmodel,'MNI') == 0
        coordinate_mni = atlasViewer.fwmodel.MNI;
        no_mni = size(coordinate_mni,1);
        foo = num2str(coordinate_mni);
        fooi = sprintf(foo(1,:));
        if no_mni == 1;
            defaultanswer = {[fooi]};
        else
            for i = 2:no_mni;
                foon = ['; ' sprintf(foo(i,:))];
                defaultanswer = {[fooi foon]};
                fooi = [fooi foon];
            end %defaultanswer = {[sprintf(foo(1,:)) ';' sprintf(foo(2,:))]}
        end
    else
        defaultanswer={'[0 0 0; 0 0 0]'};
    end
    answer=inputdlg(prompt,name,numlines,defaultanswer,'on');
    if isempty(answer)
        return
    end
    coordinate_mni = str2num(answer{1});
    no_mni = size(coordinate_mni,1);
    
    % get coordinate in MC space
    T_headvol2mc     = atlasViewer.headvol.T_2mc; % colin to Monte Carlo space
    T_labelssurf2vol = atlasViewer.labelssurf.T_2vol; % mni to colin
    for i = 1: no_mni;
        coordinate_MC(i,:) = xform_apply(coordinate_mni(i,:), (T_headvol2mc*T_labelssurf2vol));
    end
    
    % PLOT
    % clean if there is rays or dots left corresponding to previous points
    h = findobj('type','line','tag','MNI projection','color','m');
    set(h,'Visible','off');
    h2 = findobj('Marker','o');
    set(h2,'Visible','off');
    
    headvol   = atlasViewer.headvol;
    vertices  = atlasViewer.labelssurf.mesh.vertices;
    
    % Project MNI in MC space to head surface and pial surface
    for i = 1:no_mni;
        
        coordinate_MC_headsurf = ProjectionBI(coordinate_MC(i,:), atlasViewer.headsurf.mesh.vertices);
        coordinate_MC_cortex = ProjectionBI(coordinate_MC_headsurf, vertices);
        
        % Generate rays from cortex to head surface
        
        p1 = coordinate_MC_headsurf;
        p2 = coordinate_MC_cortex;
        hRays = drawRayProjection(p1,p2,headsurf);
        hold on;
        plot3(coordinate_MC_cortex(2),coordinate_MC_cortex(1),coordinate_MC_cortex(3),'bo','MarkerSize',15,'MarkerFaceColor','b');
        hold off;
    end
    
else % cleans ray(s) and blue spheres from the display.
    
    h1 = findobj('MarkerFaceColor','b','MarkerSize',15); %
    h2 = findobj('type','line','tag','MNI projection');
    set(h1,'Visible','off')
    set(h2,'Visible','off')
    
end




% --------------------------------------------------------------------
function menuItemSetMCApp_Callback(hObject, eventdata, handles)
global atlasViewer
fwmodel = atlasViewer.fwmodel;

% Last resort: If none of the above locate MC app then ask user where it is. 
while 1
    pause(.1)
    [filenm, pathnm] = uigetfile({'*'; '*.exe'}, ['Monte Carlo executable not found. Please select Monte Carlo executable.']);
    if filenm==0
        return;
    end
    
    % Do a few basic error checks
    if istextfile(filenm)
        q = menu('Selected file not an executable. Try again', 'OK', 'Cancel');
        if q==2
            return;
        else
            continue;
        end
    end    
    break;
end
[mc_exename, mc_appname, ext] = searchDirForMCApp(pathnm);
if ~isempty(mc_exename)
    fwmodel.mc_rootpath = pathnm;
    fwmodel.mc_exepath = pathnm;
    fwmodel.mc_exename = mc_exename;
    fwmodel.mc_appname = mc_appname;
    fwmodel.mc_exename_ext = ext;
end

% Set MC options based on app type
fwmodel = setMCoptions(fwmodel);

atlasViewer.fwmodel = fwmodel;



% --------------------------------------------------------------------
function menuItemSaveAnatomy_Callback(hObject, eventdata, handles)
global atlasViewer

headvol      = atlasViewer.headvol;
headsurf     = atlasViewer.headsurf;
pialsurf     = atlasViewer.pialsurf;
labelssurf   = atlasViewer.labelssurf;
refpts       = atlasViewer.refpts;
dirnameSubj  = atlasViewer.dirnameSubj;

saveHeadvol(headvol, dirnameSubj);
saveHeadsurf(headsurf, dirnameSubj, headvol.T_2mc);
savePialsurf(pialsurf, dirnameSubj, headvol.T_2mc);
saveLabelssurf(labelssurf, dirnameSubj, headvol.T_2mc);
saveRefpts(refpts, dirnameSubj, headvol.T_2mc);



% --------------------------------------------------------------------
function menuItemLoadPrecalculatedProfile_Callback(hObject, eventdata, handles)
global atlasViewer

fwmodel = atlasViewer.fwmodel;
headvol = atlasViewer.headvol;
pialsurf = atlasViewer.pialsurf;
probe = atlasViewer.probe;
T_vol2mc = headvol.T_2mc;

dirnameSubj = atlasViewer.dirnameSubj;

pathnm = [headvol.pathname, 'fw/'];

% Check if there's a fluence profile which already exists

for ii=1:length(fwmodel.fluenceProfFnames)
    if ~exist(fwmodel.fluenceProfFnames{ii}, 'file')
        fwmodel.fluenceProfFnames={};
        break;
    end
end
if isempty(fwmodel.fluenceProfFnames)
    pathnm = getDesktopDir();
    fluenceProfFnames = dir([pathnm, 'fluenceProf*.mat']);
    for ii=1:length(fluenceProfFnames)
        foo = loadFluenceProf([pathnm, fluenceProfFnames(ii).name], 'index');
        fwmodel.fluenceProfFnames{foo.index} = [pathnm, fluenceProfFnames(ii).name];
    end
end
if isempty(fwmodel.fluenceProfFnames)
    fluenceProfFnames = dir([pathnm, 'fluenceProf*.mat']);
    while isempty(fluenceProfFnames)
        q = menu('No profile associated with this anatomy. Do you want to find the folder that contains profiles?','Yes','No');
        if q==2
            return;
        else
            pause(.1);
            pathnm = uigetdir(dirnameSubj, 'Search for the fluence profile folder.');
        end
        if pathnm(end)~='/' && pathnm(end)~='\'
            pathnm(end+1)='/';
        end
        fluenceProfFnames = dir([pathnm, 'fluenceProf*.mat']);
    end
    
    for ii=1:length(fluenceProfFnames)
        foo = loadFluenceProf([pathnm, fluenceProfFnames(ii).name], 'index');
        fwmodel.fluenceProfFnames{foo.index} = [pathnm, fluenceProfFnames(ii).name];
    end
end

fwmodel = genSensitivityProfileFromFluenceProf(fwmodel, probe, T_vol2mc, dirnameSubj);

atlasViewer.fwmodel = fwmodel;
menuItemGenerateLoadSensitivityProfile_Callback([], struct('EventName','profile'), handles);




% --------------------------------------------------------------------
function popupmenuImageDisplay_Callback(hObject, eventdata, handles)
global atlasViewer

fwmodel  = atlasViewer.fwmodel;
imgrecon = atlasViewer.imgrecon;
pialsurf = atlasViewer.pialsurf;
axesv    = atlasViewer.axesv;

strs = get(hObject,'string');
val = get(hObject,'value');

if val==1
    % First turn off metrics display
    imgrecon = showImgReconDisplay(imgrecon, axesv(1).handles.axesSurfDisplay, 'off', 'off', 'off', 'off');
    
    % Now display sensitivity
    fwmodel = showFwmodelDisplay(fwmodel, axesv(1).handles.axesSurfDisplay, 'on');
elseif val==2
    % Turn sensitivity display off
    fwmodel = showFwmodelDisplay(fwmodel, axesv(1).handles.axesSurfDisplay, 'off');
    
    % Turn localization error on and resolution display off
    imgrecon = showImgReconDisplay(imgrecon, axesv(1).handles.axesSurfDisplay, 'on', 'off', 'off', 'off');
elseif val==3
    % Turn sensitivity display off
    fwmodel = showFwmodelDisplay(fwmodel, axesv(1).handles.axesSurfDisplay, 'off');

    % Turn resolution on and localization error display off
    imgrecon = showImgReconDisplay(imgrecon, axesv(1).handles.axesSurfDisplay, 'off', 'on', 'off', 'off');
elseif val==4
    % Turn sensitivity display off
    fwmodel = showFwmodelDisplay(fwmodel, axesv(1).handles.axesSurfDisplay, 'off');

    % Turn resolution on and localization error display off
    imgrecon = showImgReconDisplay(imgrecon, axesv(1).handles.axesSurfDisplay, 'off', 'off', 'on', 'off');
elseif val==5
    % Turn sensitivity display off
    fwmodel = showFwmodelDisplay(fwmodel, axesv(1).handles.axesSurfDisplay, 'off');

    % Turn resolution on and localization error display off
    imgrecon = showImgReconDisplay(imgrecon, axesv(1).handles.axesSurfDisplay, 'off', 'off', 'off', 'on');
elseif val==6
    % Turn sensitivity display off
    fwmodel = showFwmodelDisplay(fwmodel, axesv(1).handles.axesSurfDisplay, 'off');

    % Turn resolution on and localization error display off
    imgrecon = showImgReconDisplay(imgrecon, axesv(1).handles.axesSurfDisplay, 'off', 'off', 'off', 'off');
end

set(pialsurf.handles.radiobuttonShowPial, 'value',0);
uipanelBrainDisplay_Callback(pialsurf.handles.radiobuttonShowPial, [], handles);




% --------------------------------------------------------------------
function popupmenuImageDisplay_CreateFcn(hObject, eventdata, handles)
set(hObject, 'string',{'Sensitivity','Localization Error','Resolution','HbO','HbR','None'});
set(hObject, 'value',1);




% --------------------------------------------------------------------
function editSelectChannelSensitivity_new_Callback(hObject, eventdata, handles)
global atlasViewer

fwmodel      = atlasViewer.fwmodel;
probe        = atlasViewer.probe;
mesh         = atlasViewer.pialsurf.mesh;
headsurf     = atlasViewer.headsurf;
pialsurf     = atlasViewer.pialsurf;


Ch = str2num(get(hObject,'string'));
if isempty(Ch)
    set(hObject,'string',num2str(fwmodel.Ch));
    return;
end
if length(Ch)~=2
    set(hObject,'string',num2str(fwmodel.Ch));
    return;
end
iCh = find(probe.ml(:,1)==Ch(1) & probe.ml(:,2)==Ch(2), 1);
if Ch(1)==0 & Ch(2)==0
    iCh = 0;
end
if isempty(iCh)
    set(hObject,'string',num2str(fwmodel.Ch));
    return;
end
fwmodel.Ch = Ch;
fwmodel = displaySensitivity(fwmodel,pialsurf,[],probe);
atlasViewer.fwmodel = fwmodel;




% --------------------------------------------------------------------
function editSensitivityColormapThreshold_new_Callback(hObject, eventdata, handles)
global atlasViewer
fwmodel = atlasViewer.fwmodel;
axesv = atlasViewer.axesv; 
probe = atlasViewer.probe; 

fwmodel = setSensitivityColormap(fwmodel, axesv(1).handles.axesSurfDisplay);
atlasViewer.fwmodel = fwmodel;



% --------------------------------------------------------------------
function editMetricsColormapThreshold_new_Callback(hObject, eventdata, handles)
global atlasViewer
imgrecon = atlasViewer.imgrecon;
axesv = atlasViewer.axesv; 

imgrecon = setImgReconMetricsColormap(imgrecon, axesv(1).handles.axesSurfDisplay);
atlasViewer.imgrecon = imgrecon;



% --------------------------------------------------------------------
function editImageReconColormapThreshold_Callback(hObject, eventdata, handles)
global atlasViewer
imgrecon = atlasViewer.imgrecon;
axesv = atlasViewer.axesv; 

imgrecon = setImgReconColormap(imgrecon, axesv(1).handles.axesSurfDisplay);
atlasViewer.imgrecon = imgrecon;



% --------------------------------------------------------------------
function pushbuttonOpticalPropertiesSet_new_Callback(hObject, eventdata, handles)
global atlasViewer;
headvol = atlasViewer.headvol;
fwmodel = atlasViewer.fwmodel;

% Set tissue properties
name='Input Head Optical Properties';
numlines=1;
if ~isempty(fwmodel.headvol.tiss_prop)
    tiss_prop = fwmodel.headvol.tiss_prop;
else
    tiss_prop = headvol.tiss_prop;
end

prompt = {};
defaultanswer = {};
outstr = {};
for ii=1:length(tiss_prop)
    switch lower(tiss_prop(ii).name)
        case {'skin', 'scalp'}
            prompt{end+1} = 'Scalp Scattering (1/mm)';
            prompt{end+1} = 'Scalp Absroption (1/mm)';
            outstr{end+1} = 'Scalp:';
        case {'skull', 'bone'}
            prompt{end+1} = 'Skull Scattering (1/mm)';
            prompt{end+1} = 'Skull Absroption (1/mm)';
            outstr{end+1} = 'Skull:';
        case {'dm' , 'dura mater'}
            prompt{end+1} = 'Dura Scattering (1/mm)';
            prompt{end+1} = 'Dura Absroption (1/mm)';
            outstr{end+1} = 'Dura:';
        case {'csf', 'cerebral spinal fluid'}
            prompt{end+1} = 'CSF Scattering (1/mm)';
            prompt{end+1} = 'CSF Absroption (1/mm)';
            outstr{end+1} = 'CSF:';
        case {'gm', 'gray matter'}
            prompt{end+1} = 'Gray Scattering (1/mm)';
            prompt{end+1} = 'Gray Absroption (1/mm)';
            outstr{end+1} = 'Gray:';
        case {'wm', 'white matter'}
            prompt{end+1} = 'White Scattering (1/mm)';
            prompt{end+1} = 'White Absroption (1/mm)';
            outstr{end+1} = 'White:';
        case 'other'
            prompt{end+1} = 'Other Scattering (1/mm)';
            prompt{end+1} = 'Other Absroption (1/mm)';
            outstr{end+1} = 'Other:';
        otherwise
            prompt{end+1} = 'Huh Scattering (1/mm)';
            prompt{end+1} = 'Huh Absroption (1/mm)';
            outstr{end+1} = 'Huh:';
    end
       
    defaultanswer{end+1} = num2str(tiss_prop(ii).scattering);
    defaultanswer{end+1} = num2str(tiss_prop(ii).absorption);
end

answer = inputdlg(prompt,name,numlines,defaultanswer,'on');

if ~isempty(answer)
    jj = 0;
    nw = [];
    foos = [];
    for ii=1:length(tiss_prop)
        jj=jj+1;
        foo = str2num(answer{jj});
        nw(ii,1) = length(foo);
        tiss_prop(ii).scattering = foo;
        outstr{ii} = [outstr{ii} ' ' answer{jj}];
        
        jj=jj+1;
        foo = str2num(answer{jj});
        nw(ii,2) = length(foo);
        tiss_prop(ii).absorption = foo;
        outstr{ii} = [outstr{ii} '; ' answer{jj}];
        
        foos = sprintf('%s%s\n',foos,outstr{ii});
    end
    if length(unique(nw))>1
        errordlg('You must enter the same number of wavelengths for each property')
        return;
    end
    
    %set(handles.textOpticalProperties,'string',foos);
    
    headvol.tiss_prop = tiss_prop;
    fwmodel.headvol.tiss_prop = tiss_prop;
    fwmodel.nWavelengths = nw(1);

end

% Set number of photons
answer = inputdlg_errcheck({'Number of photons'},'Number of Photons', 1, {num2str(fwmodel.nphotons)});
if ~isempty(answer)
    fwmodel.nphotons = str2num(answer{1});
end

atlasViewer.headvol = headvol;
atlasViewer.fwmodel = fwmodel;




% --------------------------------------------------------------------
function pushbuttonCalcMetrics_new_Callback(hObject, eventdata, handles)
global atlasViewer

imgrecon     = atlasViewer.imgrecon;
fwmodel 	 = atlasViewer.fwmodel;
pialsurf     = atlasViewer.pialsurf;
probe        = atlasViewer.probe;
dirnameSubj  = atlasViewer.dirnameSubj;

imgrecon = inputParamsImgRecon(imgrecon);
if isempty(imgrecon)
    return;
end

% Set image popupmenu to resolution 
set(imgrecon.handles.popupmenuImageDisplay,'value',3);

% Turn off sensitivity display
set(fwmodel.handles.surf,'visible','off');
setSensitivityColormap(fwmodel, []);

set(imgrecon.handles.hHbO,'visible','off');
set(imgrecon.handles.hHbR,'visible','off');
setImgReconColormap(imgrecon, []);


imgrecon = genImgReconMetrics(imgrecon, fwmodel, dirnameSubj);
imgrecon = displayImgRecon(imgrecon, fwmodel, pialsurf, [], probe);

atlasViewer.imgrecon = imgrecon;



% --------------------------------------------------------------------
function menuItemImageReconGUI_Callback(hObject, eventdata, handles)
global atlasViewer

ImageRecon();




% --------------------------------------------------------------------
function radiobuttonShowRefpts_Callback(hObject, eventdata, handles)

radiobuttonShowRefpts(hObject, eventdata, handles)



% --------------------------------------------------------------
function uipanelBrainDisplay_Callback(hObject, eventdata, handles)

uipanelBrainDisplay(hObject, eventdata, handles);


