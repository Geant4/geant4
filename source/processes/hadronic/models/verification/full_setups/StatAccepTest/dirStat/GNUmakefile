SRCS =  pvalue.cpp        

OBJ_DIR = ${PWD}
OBJS    = ${addprefix ${OBJ_DIR}/,${SRCS:.cpp=.o}}
EXES    = ${addprefix ${OBJ_DIR}/,${SRCS:.cpp=}}

###C_INCLUDES += -I${STAT_DIR} -I${GSL_INCLUDES}
###C_INCLUDES += -I${STAT_DIR}/Comparison/include 
###C_INCLUDES += -I${STAT_DIR}/StatisticsTesting 
C_INCLUDES += -I${PI_DIR}/include
C_INCLUDES += -I${PI_DIR}/include/StatisticsTesting
C_INCLUDES += -I${GSL_DIR}/include/gsl 
C_INCLUDES += `aida-config --includes`

GSL_LIBS=-L${GSL_DIR}/lib -lgsl -lgslcblas

LIBS += -llcg_StatisticsTesting `aida-config --libs` ${GSL_LIBS} 

C_FLAGS += --pic 

ifdef DEBUG_PRINT
C_FLAGS += -DDEBUG_PRINT
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
	g++ ${C_FLAGS} ${C_INCLUDES} $< ${LIBS} -o $@

