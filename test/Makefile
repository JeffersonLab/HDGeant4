INC = ../src
SRC = ../src

options_t: options_t.cc $(SRC)/GlueXUserOptions.cc
	g++ -g -I $(INC) -o $@ $^	

t: t.C
	g++ -O4 -g -std=c++11 -DHDF5_SUPPORT -o $@ \
		-I ${HDF5ROOT}/include \
		-I ${HALLD_RECON_HOME}/${BMS_OSNAME}/include $^ \
		-L ${HALLD_RECON_HOME}/${BMS_OSNAME}/lib -lHDDM -lxstream \
		-L $(HDF5ROOT)/lib -lhdf5_cpp -lhdf5_hl -lhdf5 -lsz -lz -lbz2 -ldl 

hddm_x.hpp: x.xml
	hddm-cpp $^

hddm_x.cpp: x.xml
	hddm-cpp $^

x: x.C hddm_x.cpp
	g++ -O0 -g -std=c++11 -DHDF5_SUPPORT -o $@ \
		-I ${HDF5ROOT}/include -I. \
		-I ${HALLD_RECON_HOME}/${BMS_OSNAME}/include $^ \
		-L ${HALLD_RECON_HOME}/${BMS_OSNAME}/lib -lHDDM -lxstream \
		-L $(HDF5ROOT)/lib -lhdf5_cpp -lhdf5_hl -lhdf5 -lsz -lz -lbz2 -ldl 
