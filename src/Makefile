#Disable building with rendering by passing RENDER= to make
RENDER = yes

DEFINES=-D_GNU_SOURCE
#MARCH=native
MARCH=core2

OBJECTS = system.o math.o task.o measure.o samplers.o physics.o extraInteractions.o integrator.o monteCarlo.o world.o spgrid.o tinymt/tinymt64.o render.o octave.o
EXTRA_RENDER_OBJECTS = font.o mathlib/vector.o mathlib/quaternion.o mathlib/matrix.o

LIBS = -lm
EXTRA_RENDER_LIBS = -lfreetype -lSDL -lGL

ifeq ($(RENDER), yes)
	OBJECTS += $(EXTRA_RENDER_OBJECTS)
	LIBS += $(EXTRA_RENDER_LIBS)
else
	DEFINES += -DNO_RENDER
endif

WARNINGS = -pedantic -Wextra -Wall -Wwrite-strings -Wshadow -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes
PROFILE = 
DEBUG = -DNDEBUG
GCC_OPTIM = -O4 -flto -fwhole-program -mmmx -msse -msse2 -msse3 -fexcess-precision=fast -ffast-math -finline-limit=2000 -fmerge-all-constants -fmodulo-sched -fmodulo-sched-allow-regmoves -fgcse-sm -fgcse-las -fgcse-after-reload -funsafe-loop-optimizations
GCC_CFLAGS = $(WARNINGS) -std=gnu99 -pipe -march=$(MARCH) -ggdb $(GCC_OPTIM) $(DEBUG) $(DEFINES) -I/usr/include/freetype2

CLANG_OPTIM = -O4 -march=$(MARCH)
CLANG_CFLAGS = $(CLANG_OPTIM) $(DEBUG) $(DEFINES) -ggdb -I/usr/include/freetype2 


#use 'make target CLANG=' to build with clang instead of gcc
CLANG=no
ifneq ($(CLANG), no)
	CC = clang
	CFLAGS = $(CLANG_CFLAGS)
else
	CC = gcc
	CFLAGS = $(GCC_CFLAGS)
endif

#use 'make target STATIC=' to generate a statically linked binary
STATIC=no
ifneq ($(STATIC), no)
	STATIC_FLAGS = -static
endif

#use 'make target LANGEVIN=' to build with alterative Langevin integrator
LANGEVIN=regular
ifneq ($(LANGEVIN), regular)
	DEFINES += -DALTERNATIVE_LANGEVIN
endif


all: hairpin

profile:
	make -B PROFILE="-pg"

#Debug (asserts etc) *without* optimizations
debug:
	make -B GCC_OPTIM="-O0" CLANG_OPTIM="-O0" DEBUG="-DDEBUG"

#Debug (asserts etc) *with* optimizations
fastdebug:
	make -B DEBUG="-DDEBUG"



doubleStrand: $(OBJECTS) doubleStrand.o
	@echo "LD	doubleStrand"
	@$(CC) -o doubleStrand $(CFLAGS) $(OBJECTS) doubleStrand.o $(LIBS) -ggdb $(PROFILE)


hairpin: $(OBJECTS) hairpin.o
	@echo "LD	hairpin"
	@$(CC) $(STATIC_FLAGS) -o hairpin $(CFLAGS) $(OBJECTS) hairpin.o $(LIBS) -ggdb $(PROFILE)
hairpindebug:
	make -B GCC_OPTIM="-O0" CLANG_OPTIM="-O0" DEBUG="-DDEBUG" hairpin
hairpinCombined: hairpinCombined.o
	@echo "LD	hairpinCombined"
	@$(CC) $(STATIC_FLAGS) -o hairpinCombined $(CFLAGS) hairpinCombined.o $(LIBS) -ggdb $(PROFILE)
hairpinCombinedProfile:
	make -B PROFILE="-pg" hairpinCombined


diff: $(OBJECTS) diffusion.o
	@echo "LD	diffusion"
	@$(CC) $(STATIC_FLAGS) -o diffusion $(CFLAGS) $(OBJECTS) diffusion.o $(LIBS) -ggdb $(PROFILE)
diffdebug:
	make -B GCC_OPTIM="-O0" CLANG_OPTIM="-O0" DEBUG="-DDEBUG" diff
difffastdebug:
	make -B DEBUG="-DDEBUG" diff


tests: $(OBJECTS) tests.o
	@echo "LD	hairpin"
	@$(CC) $(STATIC_FLAGS) -o tests $(CFLAGS) $(OBJECTS) tests.o $(LIBS) -ggdb $(PROFILE)
testsdebug:
	make -B GCC_OPTIM="-O0" CLANG_OPTIM="-O0" DEBUG="-DDEBUG" tests


.c.o:
	@echo "CC	$@"
	@$(CC) -o $@ $(CFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(EXTRA_RENDER_OBJECTS) doubleStrand.o hairpin.o diffusion.o tests.o

