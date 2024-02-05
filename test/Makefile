#
# Makefile for general tests of the core GlueX simulation functionality,
# spanning the full length of the production pipeline
#
# event generation -> simulation -> smearing -> reconstruction -> analysis
#
# author: richard.t.jones at uconn.edu
# version: november 22, 2023

##### sample bggen_example.spec
# LIST
# C
# C ===    INPUT file for BGGEN
# C
# TRIG  10000   number of events to simulate
# C             We expect 395kHz of hadronic rate at high luminosity
# C -- writing out events
# C        HDDM  simple  ntuple
# WROUT    1   0   0   hddm output only
# NPRIEV  5       only print first few events
# EPHLIM  8.4 9.0 energy range in GeV
# RNDMSEQ   9	      random number sequence
# EELEC     12.          electron beam energy
# EPEAK      9.0         coherent peak energy
# ZCOLLIM   7600.        istance to the collimator in cm
# DCOLLIM  0.005         collimator diameter in m
# EPYTHMIN     3.        minimal energy for PYTHIA simulation
# RUNNO    50141         run number for output events
# STOP
##### end of sample bggen_example.spec

##### sample genr8_example.spec
# % run 50141      ! this line is mandatory
# % events 10000   ! this line is mandatory
# % seed 477023945 ! this line is mandatory
# % beamp.x beamp.y beamp.z beamMass
# 0 0 8.5 0
# % targetp.x targetp.y targetp.z targetMass
# 0 0 0 0.938
# % t-channelSlope
#       5.0
#       % number of particles needed to describe the isobar decay of X
#       8
#       %
#       % particle# 0&1 are always the X&Y
# %part#  chld1#  chld2#  parent# Id     nchild   mass    width   charge  flag
# % baryon (Y) decay
#  0       *      *       *       14       0      0.9382   0       +1   11
# % meson (X) decay
#  1       2      5       *      ??       2       2.0      5.0      0     00
#  2       3      4       1      17       2       0.548    0.001    0     00
#  3       *      *       2       1       0       0        0        0     10
#  4       *      *       2       1       0       0        0        0     10
#  5       *      *       1       7       2       0.135    0        0     00
#  6       *      *       5       1       0       0        0        0     10
#  7       *      *       5       1       0       0        0        0     10
# !EOI
##### end of sample genr8_example.spec

##### sample reaction_example.spec
# PLUGINS ReactionFilter
# #Ncores if you want to use all of the hyperthreads available (only use if you have the whole node!)
# #NTHREADS Ncores
# THREAD_TIMEOUT 600
# REST:DATAVERSIONSTRING recon_RunPeriod-2019-11_ver01
# # This was the original timestamp. It had to be updated on 
# # 8/20/2020 due to new TOF calibrations. Runs 71350-71405
# # and 71589-71591 at NERSC were processed with this timestamp.
# #JANA_CALIB_CONTEXT calibtime=2020-07-24-00-00-01
# #calibtime=YYYY-MM-DD-hh-mm-ss (hr is 24-hr format)
# # This is used for runs 71592-71620 (and hopefully all others) at PSC
# #JANA_CALIB_CONTEXT calibtime=2020-09-28-00-00-01
# Reaction1 1_14__17_7_14
# Reaction1:Flags B4_T1_S1_M7_M17_F4
##### end of sample reaction_example.spec

comma = ,

ifdef BGGEN_URL
    bggen_url = $(BGGEN_URL)
else
    bggen_url = https://zeus.phys.uconn.edu/halld/bggen-3-2023/config
endif

ifdef NEVENTS
    nevents = $(NEVENTS)
else
    nevents = 999999999999
endif

ifdef NTHREADS
    nthreads = $(NTHREADS)
else
    nthreads = 6
endif

ifdef BATCHMODE
    batchmode = $(BATCHMODE)
else
    batchmode = 0
endif

ifeq ($(VALGRIND), "gdb")
    exec = valgrind --vgdb-error=0
else ifdef VALGRIND
    exec = valgrind --leak-check=full
else ifdef GDB
    exec = @function run_with_gdb {\
	       executable=$$1;\
               shift;\
	       echo $$executable $$*;\
	       gdb $$executable;\
	   };\
	   run_with_gdb
else
    exec = 
endif
ifdef VERBOSE
    g4options = -V$(VERBOSE) -t$(nthreads)
else
    g4options = -t$(nthreads)
endif

all:
	@echo "GENERAL TESTS OF CORE GLUEX SIMULATION+ANALYSIS FUNCTIONALITY"
	@echo "A list of valid test targets:"
	@echo " bggen_%.hddm (requires user input file bggen_%.spec to be present)"
	@echo " genr8_%.hddm (requires user input file genr8_%.spec to be present)"
	@echo " hdgeant4_%.hddm (uses control file control_%.in, if present)"
	@echo " hdgeant31_%.hddm (uses control file  control_%.in, if present)"
	@echo " hdgeant34_%.hddm (uses control file  control_%.in, if present)"
	@echo " smeared_%.hddm (assumes % starts with either hdgeant3*_ or hdgeant4_)"
	@echo " rest_%.hddm (reads from smeared_%.hddm)"
	@echo " reaction_<reaction>_%.root (assumes reaction_<reaction>.spec is present)"
	@echo "  where <reaction> can be any of the following:"
	@echo "   etaS2gS,pi0,p"
	@echo "   ... to be extended ..."

bggen_%.hddm: bggen_%.spec
	@test -f pythia-geant.map || wget $(bggen_url)/pythia-geant.map || exit 1
	@test -f pythia.dat || wget $(bggen_url)/pythia.dat || exit 1
	@test -f particle.dat || wget $(bggen_url)/particle.dat || exit 1
	@cp $< fort.15 || exit 1
	$(exec) bggen || exit 1
	@mv bggen.hddm $@ || exit 1
	@rm fort.15

genr8_%.hddm: genr8_%.spec
	@echo "genr8 -Agenr8_output.txt \
	 -M$$(awk '/^. events /{print $$3}' $<) \
	 -r$$(awk '/^. run /{print $$3}' $<) \
	 -s$$(awk '/^. seed /{print $$3}' $<) < $<"
	$(exec) genr8 -Agenr8_output.txt \
	 -M$$(awk '/^. events /{print $$3}' $<) \
	 -r$$(awk '/^. run /{print $$3}' $<) \
	 -s$$(awk '/^. seed /{print $$3}' $<) < $< || exit 1
	genr8_2_hddm -V"0 0 50 80" genr8_output.txt || exit 1
	@mv genr8_output.hddm $@
	@rm -f genr8_output.txt

control_%.in: control.in_test_template
	@sed "s/^INFI.*/INFILE '$*.hddm'/" $< > $@ || exit 1

hdgeant4_%.hddm: control_%.in %.hddm
	@sed "s/^OUTF.*/OUTFILE '$@'/" < $< > control.in
	$(exec) hdgeant4 $(g4options) || exit 1

hdgeant31_%.hddm: control_%.in %.hddm
	@sed "s/^OUTF.*/OUTFILE '$@'/" $< > control.in || exit 1
	@sed -i "s/^HADR .*/HADR 1/" control.in || exit 1
	@echo hdgeant -xml=ccdb://GEOMETRY/main_HDDS.xml,run=$$( \
	 awk 'BEGIN{run=0}/^RUNNO /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 awk 'BEGIN{run=0}/^RUNG /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 hddm-xml -n 1 $*.hddm | awk -F'"' '/runNo/{print $$4}')
	$(exec) hdgeant -xml=ccdb://GEOMETRY/main_HDDS.xml,run=$$( \
	 awk 'BEGIN{run=0}/^RUNNO /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 awk 'BEGIN{run=0}/^RUNG /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 hddm-xml -n 1 $*.hddm | awk -F'"' '/runNo/{print $$4}') || \
	 exit 1

hdgeant34_%.hddm: control_%.in %.hddm
	@sed "s/^OUTF.*/OUTFILE '$@'/" $< > control.in || exit 1
	@sed -i "s/^HADR .*/HADR 4/" control.in || exit 1
	@echo hdgeant -xml=ccdb://GEOMETRY/main_HDDS.xml,run=$$( \
	 awk 'BEGIN{run=0}/^RUNNO /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 awk 'BEGIN{run=0}/^RUNG /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 hddm-xml -n 1 $*.hddm | awk -F'"' '/runNo/{print $$4}')
	$(exec) hdgeant -xml=ccdb://GEOMETRY/main_HDDS.xml,run=$$( \
	 awk 'BEGIN{run=0}/^RUNNO /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 awk 'BEGIN{run=0}/^RUNG /{run=$$2}END{if(run>0){print run}else{exit 1}}' $< || \
	 hddm-xml -n 1 $*.hddm | awk -F'"' '/runNo/{print $$4}') || \
	 exit 1

smeared_%.hddm: %.hddm
	$(exec) mcsmear \
        -PNTHREADS=$(nthreads) \
        -PEVENTS_TO_KEEP=$(nevents) \
	-Pprint -PJANA:BATCH_MODE=$(batchmode) \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        $< || exit 1
	@mv $*_smeared.hddm $@

rest_%.hddm: hd_recon.config smeared_%.hddm
	$(exec) hd_root \
	--config=$< \
        --nthreads=$(nthreads) \
        -PNTHREADS=$(nthreads) \
        -PEVENTS_TO_KEEP=$(nevents) \
        -PJANA:BATCH_MODE=$(batchmode) \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        -PTRK:SAVE_TRUNCATED_DEDX=1 \
        smeared_$*.hddm || exit 1
	@mv dana_rest.hddm $@

pi0eta2g = etaS2gS$(comma)pi0$(comma)p
pi0eta3pi0 = etaS3pi0S$(comma)pi0$(comma)p
pi0eta3pi = etaSpi+$(comma)pi-$(comma)pi0S$(comma)pi0$(comma)p

reaction_$(pi0eta2g)_%.root: reaction_$(pi0eta2g).spec rest_%.hddm
	$(exec) hd_root \
	--config=$< \
        --nthreads=$(nthreads) \
        -PNTHREADS=$(nthreads) \
        -PEVENTS_TO_KEEP=$(nevents) \
        -PJANA:BATCH_MODE=$(batchmode) \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        rest_$*.hddm || exit 1
	@mv hd_root.root $@

reaction_$(pi0eta3pi0)_%.root: reaction_$(pi0eta3pi0).spec rest_%.hddm
	$(exec) hd_root \
	--config=$< \
        --nthreads=$(nthreads) \
        -PNTHREADS=$(nthreads) \
        -PEVENTS_TO_KEEP=$(nevents) \
        -PJANA:BATCH_MODE=$(batchmode) \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        rest_$*.hddm || exit 1
	@mv hd_root.root $@

reaction_$(pi0eta3pi)_%.root: reaction_$(pi0eta3pi).spec rest_%.hddm
	$(exec) hd_root \
	--config=$< \
        --nthreads=$(nthreads) \
        -PNTHREADS=$(nthreads) \
        -PEVENTS_TO_KEEP=$(nevents) \
        -PJANA:BATCH_MODE=$(batchmode) \
        -PTHREAD_TIMEOUT_FIRST_EVENT=3600 \
        -PTHREAD_TIMEOUT=600 \
        rest_$*.hddm || exit 1
	@mv hd_root.root $@

test:
	@make reaction_$(pi0eta2g)_hdgeant4_bggen_10k.root
	@make reaction_$(pi0eta3pi0)_hdgeant4_bggen_10k.root
	@make reaction_$(pi0eta3pi)_hdgeant4_bggen_10k.root
	@make reaction_$(pi0eta2g)_hdgeant34_bggen_10k.root
	@make reaction_$(pi0eta3pi0)_hdgeant34_bggen_10k.root
	@make reaction_$(pi0eta3pi)_hdgeant34_bggen_100k.root
	@make reaction_$(pi0eta2g)_hdgeant4_genr8_eta,pi0,p_10k.root
	@make reaction_$(pi0eta3pi0)_hdgeant4_genr8_eta,pi0,p_10k.root
	@make reaction_$(pi0eta3pi)_hdgeant4_genr8_eta,pi0,p_10k.root
	@make reaction_$(pi0eta2g)_hdgeant31_genr8_eta,pi0,p_10k.root
	@make reaction_$(pi0eta3pi0)_hdgeant31_genr8_eta,pi0,p_10k.root
	@make reaction_$(pi0eta3pi)_hdgeant31_genr8_eta,pi0,p_10k.root

%.config:
	@test -f ../$@ || (echo "$@ not found, cannot continue" && exit 1)
	@cp -s ../$@ .

%.spec:
	@test -f ../$@ || (echo "$@ not found, cannot continue" && exit 1)
	@cp -s ../$@ .

control.in_test_template:
	@test -f ../$@ || (echo "$@ not found, cannot continue" && exit 1)
	@cp -s ../$@ .

clean:
	@rm -rf pythia-geant.map pythia.dat particle.dat bggen.his fort.15 *.astate control.in *.his *.hbook *.rz smear.root flukaerr.dat

.SECONDARY: $(%.hddm) $(control_%.in)
