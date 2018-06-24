function ui_plant
f = ne_dialg('Plant Chamber Instrument',1);
f = ne_dialg(f, 'add', 0, 1, 'gplantsws', 'SW Status' );
f = ne_dialg(f, 'add', 1, 0, 'pplantswssws', 'SW Stat' );
f = ne_dialg(f, 'add', 1, 0, 'pplantswsf', 'Flag' );
f = ne_dialg(f, 'add', 0, 1, 'gplantfcc', 'FCC' );
f = ne_dialg(f, 'add', 1, 0, 'pplantfccf', 'Flow 0' );
f = ne_dialg(f, 'add', 1, 0, 'pplantfccflow1', 'Flow 1' );
f = ne_dialg(f, 'add', 1, 0, 'pplantfccflow2', 'Flow 2' );
f = ne_dialg(f, 'add', 1, 0, 'pplantfcct', 'Temps' );
f = ne_dialg(f, 'add', 1, 0, 'pplantfccs', 'Status' );
f = ne_dialg(f, 'add', 0, 1, 'gplanttm', 'T Mbase' );
f = ne_dialg(f, 'add', 1, 0, 'pplanttmtd', 'T Drift' );
f = ne_dialg(f, 'add', 1, 0, 'pplanttmcpu', 'CPU' );
f = ne_dialg(f, 'add', 1, 0, 'pplanttmram', 'RAM' );
f = ne_dialg(f, 'add', 1, 0, 'pplanttmd', 'Disk' );
f = ne_listdirs(f, 'PLANT_DATA_DIR', 15);
f = ne_dialg(f, 'newcol');
ne_dialg(f, 'resize');
