WavePath

load 'D:\Philippe\SwitchDrive\AstroQUT_matlab\ASTROQUT\data_and_results\FEOchandra.mat'

filename='ChandraFEO.mat'

options=astro_setparams();
options.ps=1;
options.nitr=2000;
options.showplot=0;

tic
sol= astrosolveWS(Freal,256,Ereal,Oreal,[],options);
toc

save(filename,'sol')