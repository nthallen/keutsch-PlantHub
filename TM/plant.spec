tmcbase = base.tmc
cmdbase = /usr/local/share/huarp/phrtg.cmd
genuibase = plant.genui

tmcbase = SWStat.tmc
genuibase = SWStat.genui
swsbase = plant.sws

#tmcbase = TS.tmc TS_T_t.tmc
#colbase = TS_col.tmc
#extbase = TS_conv.tmc
#genuibase = TS.genui

tmcbase = fcc.tmc
colbase = fcc_col.tmc
cmdbase = fcc.cmd
genuibase = fcc.genui

Module TMbase

SCRIPT = interact
TGTDIR = $(TGTNODE)/home/plant
IGNORE = Makefile

plantcol : -lsubbuspp
plantsrvr : -lsubbuspp
plantdisp : fcc_conv.tmc Plant.tbl
plantrtgext : rtg.tmc /usr/local/share/oui/cic.oui
doit : plant.doit
plantalgo : fcc_conv.tmc plant.tma plant.sws
