all: model7.so 
# model7_eth.so model7_mwi.so	model7_zwe.so model7_rwa.so model7_ssd_caf.so ## Comment here preserves the prior tab
model7.so: model7.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7.cpp', '-O1 -g',DLLFLAGS='')"
# model7_eth.so: model7_eth.cpp
# 	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_eth.cpp', '-O1 -g',DLLFLAGS='')"
# model7_mwi.so: model7_mwi.cpp
# 	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_mwi.cpp', '-O1 -g',DLLFLAGS='')"
# model7_zwe.so: model7_zwe.cpp
# 	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_zwe.cpp', '-O1 -g',DLLFLAGS='')"
# model7_rwa.so: model7_rwa.cpp
# 	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_rwa.cpp', '-O1 -g',DLLFLAGS='')"
# model7_ssd_caf.so: model7_ssd_caf.cpp
# 	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_ssd_caf.cpp', '-O1 -g',DLLFLAGS='')"


clean:
	rm -rf *o
