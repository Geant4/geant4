SRCS =  pvalue.cpp        

OBJ_DIR = ${PWD}
OBJS    = ${addprefix ${OBJ_DIR}/,${SRCS:.cpp=.o}}
EXES    = ${addprefix ${OBJ_DIR}/,${SRCS:.cpp=}}

C_INCLUDES += -I${PI_DIR}/include
C_INCLUDES += -I${DIR_STAT}/include/StatisticsTesting -I${DIR_STAT}/include/UserLayer-aida 
C_INCLUDES += -I${GSL_DIR}/include/gsl 
C_INCLUDES += `aida-config --includes`

LIBS += -L${DIR_STAT}/lib -lStatisticsTesting 
LIBS += `aida-config --libs` 
LIBS += -L${GSL_DIR}/lib -lgsl -lgslcblas 

C_FLAGS += --pic 

ifdef DEBUG_PRINT
C_FLAGS += -DDEBUG_PRINT
endif

ifdef FULLY_STATIC_EXECUTABLE
  LDFLAGS  += -static
endif 

all: examples

examples: ${OBJ_DIR} ${OBJS} ${EXES}

clean:
	@echo "cleaning up"
	@rm -f ${OBJS} ${EXES}

${OBJ_DIR}:
	mkdir -p ${OBJ_DIR}

#
# rules:
#

${OBJ_DIR}/%.o: %.cpp
	@echo compiling $*
	g++ -c -o ${OBJ_DIR}/$*.o ${C_FLAGS} ${C_INCLUDES} $*.cpp

${OBJ_DIR}/%: ${OBJ_DIR}/%.o
	@echo linking $*
	g++ ${LDFLAGS} ${C_FLAGS} ${C_INCLUDES} $< ${LIBS} -o $@

