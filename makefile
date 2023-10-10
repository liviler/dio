EXE_NAME = dio
EXE_DIR = ./exe
SRC_DIR = ./src
MOD_DIR = ./src/mod
OBJ_DIR = ./src/obj

# SOURCES = $(wildcard $(SRC_DIR)/*.f90)
SOURCES = 	./src/dio.f90\
		./src/module_constants.f90\
		./src/module_globals.f90 \
		./src/module_inout.f90\
		./src/module_force.f90\
		./src/module_basis.f90\
		./src/module_mathmethods.f90\
		./src/module_nucleus.f90\
		./src/module_preparation.f90\
		./src/module_wavefunction.f90\
		./src/module_field.f90\
		./src/module_matrix.f90\

OBJECTS = $(patsubst $(SRC_DIR)/%.f90, ${OBJ_DIR}/%.o, $(SOURCES))


#default: gfortran
default: gfortran #ifort

# compiled by gfortran
gfortran: FC = gfortran
gfortran: FFLAGS = -O2 -J ${MOD_DIR} #-Wall
gfortran: printConfiguration ${EXE_NAME} printEndInformation

# compiled by ifort
ifort: FC = ifort
ifort: FFLAGS = -O2 -module:${MOD_DIR}
ifort: printConfiguration ${EXE_NAME} printEndInformation

printConfiguration:
	@echo "=============Compiling with ${FC}====================="
	@echo "src path: ${SRC_DIR}/"
	@echo "mod path: ${MOD_DIR}/"
	@echo "obj path: ${OBJ_DIR}/"
	@echo "exe path: ${EXE_DIR}/${EXE_NAME}"
	@echo "------------------------------------------------------"
	
printEndInformation:
	@echo "------------------------------------------------------"
	@echo -e "\033[32mCompilation Finished ! \033[0m "
	@echo -e "src path: \033[32m ${SRC_DIR}/ \033[0m"
	@echo -e "exe path: \033[32m ${EXE_DIR}/${EXE_NAME} \033[0m"
	@ENDPrintEndInformation=]]]]]] # make brackets can be matched

${EXE_NAME}:${OBJECTS} | ${EXE_DIR}
	@echo "compiling $@ ......"
	${FC} ${FFLAGS} -o ${EXE_DIR}/${EXE_NAME} $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | ${OBJ_DIR}  ${MOD_DIR}
	@echo "compiling  $@ ......"
	$(FC) $(FFLAGS) -c $< -o $@ 

${OBJ_DIR}:
	mkdir -p ${OBJ_DIR}

${MOD_DIR}:
	mkdir -p ${MOD_DIR}

${EXE_DIR}:
	mkdir -p ${EXE_DIR}

# Dependencies
${OBJ_DIR}/dio.o:$(filter-out ${OBJ_DIR}/dio.o, ${OBJECTS})

${OBJ_DIR}/module_globals.o: $(OBJ_DIR)/module_constants.o

${OBJ_DIR}/module_inout.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o

${OBJ_DIR}/module_force.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o

${OBJ_DIR}/module_basis.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o

${OBJ_DIR}/module_mathmethods.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o

${OBJ_DIR}/module_nucleus.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o

${OBJ_DIR}/module_preparation.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o \
								 ${OBJ_DIR}/module_inout.o ${OBJ_DIR}/module_force.o \
								 ${OBJ_DIR}/module_basis.o ${OBJ_DIR}/module_mathmethods.o \
								 ${OBJ_DIR}/module_nucleus.o

${OBJ_DIR}/module_wavefunction.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o

${OBJ_DIR}/module_field.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o \
						   ${OBJ_DIR}/module_inout.o

${OBJ_DIR}/module_matrix.o: ${OBJ_DIR}/module_constants.o ${OBJ_DIR}/module_globals.o

debug:
	@echo "src path: ${SRC_DIR}/"
	@echo "mod path: ${MOD_DIR}/"
	@echo "obj path: ${OBJ_DIR}/"
	@echo "exe path: ${EXE_DIR}/${EXE_NAME}"
	@echo "SOURCES : ${SOURCES}"
	@echo "OBJECTS : ${OBJECTS}"

clean:
	rm -f ${EXE_DIR}/${EXE_NAME} $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

deepclean:
	rm -rf ${EXE_DIR} $(OBJ_DIR) $(MOD_DIR) ${SRC_DIR}/*.mod
