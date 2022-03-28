all: ar1_trend.so ar1.so arima_trend.so arima.so rw.so rw1_trend.so model7.so model7_eth.so model7_mwi.so	model7_zwe.so model7_rwa.so ## Comment here preserves the prior tab
ar1_trend.so: ar1_trend.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('ar1_trend.cpp', '-O1 -g',DLLFLAGS='')"
ar1.so: ar1.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('ar1.cpp', '-O1 -g',DLLFLAGS='')"
arima_trend.so: arima_trend.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('arima_trend.cpp', '-O1 -g',DLLFLAGS='')"
arima.so: arima.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('arima.cpp', '-O1 -g',DLLFLAGS='')"
rw.so: rw.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('rw.cpp', '-O1 -g',DLLFLAGS='')"
rw1_trend.so: rw1_trend.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('rw1_trend.cpp', '-O1 -g',DLLFLAGS='')"
model7.so: model7.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7.cpp', '-O1 -g',DLLFLAGS='')"
model7_eth.so: model7_eth.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_eth.cpp', '-O1 -g',DLLFLAGS='')"
model7_mwi.so: model7_mwi.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_mwi_rwa.cpp', '-O1 -g',DLLFLAGS='')"
model7_zwe.so: model7_zwe.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_zwe.cpp', '-O1 -g',DLLFLAGS='')"
model7_rwa.so: model7_rwa.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('model7_rwa.cpp', '-O1 -g',DLLFLAGS='')"


clean:
	rm -rf *o
