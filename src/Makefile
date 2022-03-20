all: ar1_trend.so ar1.so arima_trend.so arima.so rw.so rw1_trend.so
	# Comment here preserves the prior tab
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


clean:
	rm -rf *o
