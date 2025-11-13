tmcbase = base.tmc
genuibase = plant.genui

tmcbase = SWStat.tmc
genuibase = SWStat.genui
swsbase = plant.sws

tmcbase = fcc.tmc
colbase = fcc_col.tmc SB.cc SB.h
cmdbase = fcc.cmd
genuibase = fcc.genui

Module TMbase mode=ignore SWSnot=
# Module LICOR
# Module savelog

SCRIPT = interact
TGTDIR = /home/plant
IGNORE = Makefile '*.exe'

plantcol : -lsubbuspp
plantsrvr : -lsubbuspp SB.cc SB.h SB.oui
# plantdisp : fcc_conv.tmc LICOR_conv.tmc Plant.tbl
plantdisp : fcc_conv.tmc Plant.tbl
# plantrtgext : rtg.tmc /usr/local/share/oui/cic.oui
doit : plant.doit
plantalgo : fcc_conv.tmc plant.tma plant.sws
%%
CXXFLAGS=-g
