FLAGS:=										\
-I include									\
-no-pie 									\
-D _DEBUG									\
-ggdb3										\
-std=c++17									\
-O0											\
-Wall										\
-Wextra										\
-Weffc++									\
-Waggressive-loop-optimizations				\
-Wc++14-compat								\
-Wmissing-declarations						\
-Wcast-align								\
-Wcast-qual									\
-Wchar-subscripts							\
-Wconditionally-supported					\
-Wconversion								\
-Wctor-dtor-privacy							\
-Wempty-body								\
-Wfloat-equal								\
-Wformat-nonliteral							\
-Wformat-security							\
-Wformat-signedness							\
-Wformat=2									\
-Winline									\
-Wlogical-op								\
-Wnon-virtual-dtor							\
-Wopenmp-simd								\
-Woverloaded-virtual						\
-Wpacked									\
-Wpointer-arith								\
-Winit-self									\
-Wredundant-decls							\
-Wshadow									\
-Wsign-conversion							\
-Wsign-promo								\
-Wstrict-null-sentinel						\
-Wstrict-overflow=2							\
-Wsuggest-attribute=noreturn				\
-Wsuggest-final-methods						\
-Wsuggest-final-types						\
-Wsuggest-override							\
-Wswitch-default							\
-Wswitch-enum								\
-Wsync-nand -Wundef							\
-Wunreachable-code							\
-Wunused									\
-Wuseless-cast								\
-Wvariadic-macros							\
-Wno-literal-suffix							\
-Wno-missing-field-initializers				\
-Wno-narrowing								\
-Wno-old-style-cast							\
-Wno-varargs								\
-Wstack-protector							\
-fcheck-new									\
-fsized-deallocation						\
-fstack-protector							\
-fstrict-overflow							\
-flto-odr-type-merging						\
-fno-omit-frame-pointer						\
-Wlarger-than=8192							\
-Wstack-usage=8192							\
-Werror=vla									\
-fsanitize=address,alignment,bool,bounds,enum,float-cast-overflow,float-divide-by-zero,integer-divide-by-zero,leak,nonnull-attribute,null,object-size,return,returns-nonnull-attribute,shift,signed-integer-overflow,undefined,unreachable,vla-bound,vptr

BINDIR:=bin
OUTPUT:=Mandelbrot
SRCDIR:=source
SOURCE:=$(wildcard ${SRCDIR}/*.cpp)
OBJECTS:=$(addsuffix .o,$(addprefix ${BINDIR}/,$(basename $(notdir ${SOURCE}))))

all: ${OUTPUT}

${OUTPUT}:${OBJECTS}
	g++ ${FLAGS} ${OBJECTS} -o ${BINDIR}/${OUTPUT} -lsfml-system -lsfml-graphics -lsfml-window
${OBJECTS}: ${SOURCE} ${BINDIR}
	$(foreach SRC,${SOURCE},$(shell g++ -c ${SRC} ${FLAGS} -o $(addsuffix .o,$(addprefix ${BINDIR}/,$(basename $(notdir ${SRC}))))))
clean:
	rm -rf ${BINDIR}
${SOURCE}:

${BINDIR}:
	mkdir ${BINDIR}
