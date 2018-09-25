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
colbase = fcc_col.tmc SB.cc SB.h
cmdbase = fcc.cmd
genuibase = fcc.genui

Module TMbase
Module LICOR

SCRIPT = interact
TGTDIR = $(TGTNODE)/home/plant
IGNORE = Makefile

plantcol : -lsubbuspp
plantsrvr : -lsubbuspp SB.cc SB.h SB.oui
plantdisp : fcc_conv.tmc LICOR_conv.tmc Plant.tbl
plantrtgext : rtg.tmc /usr/local/share/oui/cic.oui
doit : plant.doit
plantalgo : fcc_conv.tmc plant.tma plant.sws
