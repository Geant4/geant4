CXXFLAGS=-Wall -O2
OBJS=                              \
	G4Unit.o                   \
	test/unitSelftest.o       \
	test/fake1.o               \
	test/fake2.o               \
	test/utility.o             \
	test/template.o		\
	../unit.o
RUNNER=unitTest

REMOVE=rm -f

.cc.o: ${CXX} -c -o $@ ${CXXFLAGS} $<

.PHONY: all list clean

all: ${RUNNER}
	./${RUNNER}

list: ${RUNNER}
	./${RUNNER} --list

${RUNNER}: ${OBJS}
	${CXX} ${CXXFLAGS} -o$@ $^

clean:
	${REMOVE} ${RUNNER} ${OBJS} 
ca:	clean