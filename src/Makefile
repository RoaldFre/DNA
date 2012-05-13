#Disable building with rendering by passing RENDER=no to make
RENDER = yes

OBJECTS = system.o math.o task.o measure.o samplers.o physics.o world.o spgrid.o tinymt/tinymt64.o render.o
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
OPTIM = -O4 -flto -fexcess-precision=fast -ffast-math -finline-limit=2000 -fmerge-all-constants -fmodulo-sched -fmodulo-sched-allow-regmoves -fgcse-sm -fgcse-las -fgcse-after-reload -funsafe-loop-optimizations
CFLAGS = $(WARNINGS) -std=c99 -pipe -march=native -ggdb $(OPTIM) $(DEBUG) $(DEFINES) $(PROFILE) -I/usr/include/freetype2

all: main

profile:
	make -B PROFILE="-pg"

#Debug (asserts etc) *without* optimizations
debug:
	make -B OPTIM="-O0" DEBUG="-DDEBUG"

#Debug (asserts etc) *with* optimizations
fastdebug:
	make -B DEBUG="-DDEBUG"

main: $(OBJECTS) main.o
	@echo "LD	main"
	@gcc -o main $(CFLAGS) $(OBJECTS) main.o $(LIBS) -ggdb $(PROFILE)

diff: $(OBJECTS) diffusion.o
	@echo "LD	diffusion"
	@gcc -o diffusion $(CFLAGS) $(OBJECTS) diffusion.o $(LIBS) -ggdb $(PROFILE)

diffdebug:
	make -B OPTIM="-O0" DEBUG="-DDEBUG" diff

difffastdebug:
	make -B DEBUG="-DDEBUG" diff

.c.o:
	@echo "CC	$@"
	@gcc -o $@ $(CFLAGS) -c $<

clean:
	rm -f *.o mathlib/*.o tinymt/*.o

