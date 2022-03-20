all: model6.so model6_eth.so model6_mwi_rwa.so 
	# Comment here preserves the prior tab
model6.so : model6.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model6.cpp', '-O1 -g',DLLFLAGS='')"
model6_eth.so: model6_eth.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model6_eth.cpp', '-O1 -g',DLLFLAGS='')"
model6_mwi_rwa.so: model6_mwi_rwa.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model6_mwi_rwa.cpp', '-O1 -g',DLLFLAGS='')"


clean:
	rm -rf *o
