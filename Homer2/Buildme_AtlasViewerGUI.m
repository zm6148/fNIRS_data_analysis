cd ./PACKAGES/AtlasViewerGUI
Buildme('AtlasViewerGUI', {'../SDgui'});
if exist('./AtlasViewerGUI.exe','file')
    movefile AtlasViewerGUI.exe ../..
end
cd ../..



