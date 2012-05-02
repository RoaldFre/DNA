OBJECTS = font.o system.o task.o measure.o samplers.o physics.o world.o spgrid.o render.o mathlib/vector.o mathlib/quaternion.o mathlib/matrix.o main.o
WARNINGS = -pedantic -Wextra -Wall -Wwrite-strings -Wshadow -Wcast-qual -Wstrict-prototypes -Wmissing-prototypes -Wunsafe-loop-optimizations
PROFILE = 
DEBUG = -DNDEBUG
OPTIM = -O4 -flto -fexcess-precision=fast -ffast-math -finline-limit=2000 -fmerge-all-constants -fmodulo-sched -fmodulo-sched-allow-regmoves -fgcse-sm -fgcse-las -fgcse-after-reload -funsafe-loop-optimizations  ${DEBUG}
CFLAGS = $(WARNINGS) -std=c99 -pipe -march=native -ggdb $(OPTIM) $(PROFILE) -I/usr/include/freetype2

all: main

profile:
	make -B PROFILE="-pg"

#Debug (asserts etc) *without* optimizations
debug:
	make -B OPTIM="-O0" DEBUG="-DDEBUG"

#Debug (asserts etc) *with* optimizations
fastdebug:
	make -B DEBUG="-DDEBUG"

main: $(OBJECTS)
	@echo "LD	main"
	@gcc -o main $(CFLAGS) $(OBJECTS) -lfreetype -lm -lSDL -lGL -ggdb $(PROFILE)

.c.o:
	@echo "CC	$@"
	@gcc -o $@ $(CFLAGS) -c $<

clean:
	rm -f *.o mathlib/*.o

