LDFLAGS=-Wl,-rpath -Wl,/usr/local/lib -L /usr/local/lib
CPPFLAGS= -I /usr/local/include
CXXFLAGS= -Wall -g
LIBS=-ltmpp -ltm -lnort

SRC = LICOR.cc LICORoui.c
OBJ = ${addsuffix .o, ${basename $(SRC) } } 
TGT = LICOR

$(TGT) : $(OBJ)
	$(LINK.cc) -o $(TGT) $(OBJ) $(LIBS)
	usemsg $(TGT) LICORoui.c
LICORoui.c : LICOR.oui
	oui -o LICORoui.c LICOR.oui
.PHONY : clean
clean :
	rm -f *.o LICORoui.c LICOR
