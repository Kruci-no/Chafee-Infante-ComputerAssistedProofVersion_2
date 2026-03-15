# a list of all the programs in your package
PROGS = ChafeeInfante/sampleDyn ChafeeInfante/findPeriodicPoint ChafeeInfante/CAProof ChafeeInfante/C1AroundUstableOrbit
# a list of all your units to be linked with your programs
OTHERS = DissipativePDE/Algebra/algebra DissipativePDE/Set/set DissipativePDE/VectorField/vectorField DissipativePDE/SolverPDE/solverPDE DissipativePDE/VectorFieldMaker/vectorFieldMaker ChafeeInfante/GallerkinProjections/gallerkinProjections ChafeeInfante/InOut/InOut Utils/SampleDyn/sampleDyn Utils/FinderAttractingFixedPoint/finderAttractingFixedPoint ChafeeInfante/ChafeeInfanteVecField/chafeeInfanteVecField  #ChafeeInfante/ChafeeInfanteVecFieldC1/chafeeInfanteVecFieldC1
#
# directory where capd scripts are (e.g. capd-config)
CAPDBINDIR = ../CAPD/build/bin/

# setting compiler and linker flags
CAPDFLAGS = `${CAPDBINDIR}capd-config --cflags`
CAPDLIBS = `${CAPDBINDIR}capd-config --libs`
CXXFLAGS += ${CAPDFLAGS} -O2 

# directory where object and dependancy files will be created
OBJDIR = .obj/

#============ the following should not be changed =========

OTHERS_STRIPPED = $(notdir ${OTHERS}) 
OTHERS_OBJ = ${OTHERS:%=${OBJDIR}%.o}
OBJ_FILES = ${OTHERS_OBJ} ${PROGS:%=${OBJDIR}%.o}
$(info $$var is [${OBJ_FILES}])


.PHONY: all
all: ${PROGS}

# rule to link executables
${PROGS}: % : ${OBJDIR}%.o ${OTHERS_OBJ}
	${CXX} -o $@ $< ${OTHERS_OBJ} ${CAPDLIBS}

# include files with dependencies
-include ${OBJ_FILES:%=%.d}

#rule to compile .cpp files and generate corresponding files with dependencies
${OBJ_FILES}: ${OBJDIR}%.o : %.cpp
	@mkdir -p ${OBJDIR}
	@mkdir -p ${OBJDIR}$(dir $<)
	$(CXX) ${CXXFLAGS} -MT $@ -MD -MP -MF ${@:%=%.d} -c -o $@ $<

# rule to clean all object files, dependencies and executables
.PHONY: clean
clean:
	rm -f ${OBJDIR}*.o ${OBJDIR}*.o.d ${PROGS}cd