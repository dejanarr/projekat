# ------------------ Compilation options ------------------------

# Loads math library.
LIBS = -lm `pkg-config --libs opencv4` -lsystemc
GET = get
# Compiler flags for release
CFLAGS = -Wall -O0 -march=x86-64 -mmmx -msse -msse2 -msse3 `pkg-config --cflags opencv4`
# Compiler flags for debugging
# CFLAGS = -Wall -g -pg `pkg-config --cflags opencv4`
# Compiler flags for profiling
# CFLAGS = -Wall -g -O4 -pg
DEPEND= makedepend $(CFLAGS)

CC = g++

# --------------------- Code modules ----------------------------

# Source files
SRCS = main.cpp image.cpp fasthessian.cpp imload.cpp
# Object files
OBJ = main.o image.o fasthessian.o imload.o
# Object files for shared library
SOBJ = image.o fasthessian.o surf.o
# Definitions
DEFS = image.h fasthessian.h ipoint.h imload.h

# ------------------------ Rules --------------------------------
#$(SRCS):
#        $(GET) $@

surf.ln: ${OBJ}
	${CC} -o $@ ${CFLAGS} ${OBJ} ${LIBS}

# Implicit rule used by Gnu Make: $(CC) -c $(CPPFLAGS) $(CFLAGS)

imload.o: imload.h imload.cpp image.h
	${CC} -c ${CFLAGS} imload.cpp

surf.o : surf.h surf.cpp
	${CC} -c ${CFLAGS} -fPIC surf.cpp

image.o: image.h image.cpp
	${CC} -c ${CFLAGS} -fPIC image.cpp

fasthessian.o: ${DEFS} fasthessian.cpp
	${CC} -c ${CFLAGS} -fPIC fasthessian.cpp

main.o: ${DEFS} ${SRCS}
	${CC} -c  ${CFLAGS} main.cpp

clean:
	-rm *.o surf.ln ../data/*.surf

#depend: $(SRCS)
#        $(DEPEND) $(SRCS)
