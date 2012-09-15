#Disable building with rendering by passing RENDER= to make
RENDER = yes

DEFINES=-D_GNU_SOURCE

OBJECTS = system.o math.o task.o measure.o samplers.o physics.o integrator.o world.o spgrid.o tinymt/tinymt64.o render.o octave.o
EXTRA_RENDER_OBJECTS = font.o mathlib/vector.o mathlib/quaternion.o mathlib/matrix.o

LIBS = -lm
EXTRA_RENDER_LIBS = -lfreetype -lSDL -lGL

ifeq ($(RENDER), yes)
	OBJECTS += $(EXTRA_RENDER_OBJECTS)
	LIBS += $(EXTRA_RENDER_LIBS)
else
	DEFINES += -DNO_RENDER
endif

WARNINGS = -pedantic -Wextra -Wall -Wwrite-strings -Wshadow -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wunsafe-loop-optimizations
PROFILE = 
DEBUG = -DNDEBUG
GCC_OPTIM = -O4 -flto -fexcess-precision=fast -ffast-math -finline-limit=2000 -fmerge-all-constants -fmodulo-sched -fmodulo-sched-allow-regmoves -fgcse-sm -fgcse-las -fgcse-after-reload -funsafe-loop-optimizations
GCC_CFLAGS = $(WARNINGS) -std=c99 -pipe -march=native -ggdb $(GCC_OPTIM) $(DEBUG) $(DEFINES) $(PROFILE) -I/usr/include/freetype2

CLANG_OPTIM = -O4 -march=core2
CLANG_CFLAGS = $(CLANG_OPTIM) $(DEBUG) $(DEFINES) -I/usr/include/freetype2 

#use 'make target CLANG=' to build with clang instead of gcc
CLANG=no
ifneq ($(CLANG), no)
	CC = clang
	CFLAGS = $(CLANG_CFLAGS)
else
	CC = gcc
	CFLAGS = $(GCC_CFLAGS)
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

main: $(OBJECTS) main.o
	@echo "LD	main"
	@$(CC) -o main $(CFLAGS) $(OBJECTS) main.o $(LIBS) -ggdb $(PROFILE)

doubleStrand: $(OBJECTS) doubleStrand.o
	@echo "LD	doubleStrand"
	@$(CC) -o doubleStrand $(CFLAGS) $(OBJECTS) doubleStrand.o $(LIBS) -ggdb $(PROFILE)

hairpin: $(OBJECTS) hairpin.o
	@echo "LD	hairpin"
	@$(CC) -o hairpin $(CFLAGS) $(OBJECTS) hairpin.o $(LIBS) -ggdb $(PROFILE)
hairpindebug:
	make -B GCC_OPTIM="-O0" CLANG_OPTIM="-O0" DEBUG="-DDEBUG" hairpin

diff: $(OBJECTS) diffusion.o
	@echo "LD	diffusion"
	@$(CC) -o diffusion $(CFLAGS) $(OBJECTS) diffusion.o $(LIBS) -ggdb $(PROFILE)

diffdebug:
	make -B GCC_OPTIM="-O0" CLANG_OPTIM="-O0" DEBUG="-DDEBUG" diff

difffastdebug:
	make -B DEBUG="-DDEBUG" diff

.c.o:
	@echo "CC	$@"
	@$(CC) -o $@ $(CFLAGS) -c $<

clean:
	rm -f *.o mathlib/*.o tinymt/*.o

