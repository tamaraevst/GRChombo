# temporary makefile for testing this BosonStar initial data generator

# compiler setting
CC=g++
CFLAGS = -isystem /usr/include/boost_1_68_0/
CFLAGS += -I ..
CFLAGS += -I ${GRCHOMBO_SOURCE}/../Examples/BosonStar/
CFLAGS += -I ${GRCHOMBO_SOURCE}/utils/
CFLAGS += -Wall -g

test: test.o
		${CC} $< -o $@

%.o: %.cpp
		${CC} ${CFLAGS} -c $< -o $@

clean:
		rm *.o test
