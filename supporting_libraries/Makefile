#CXX = icpc
#CC = icc
#CXXFLAGS = -O2 -I . -I ../include -Wall -ffast-math -DLINUX -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

CXX = g++
CC = gcc
CXXFLAGS = -O2 -I . -I ../include -Wall -ffast-math -DLINUX -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

# Some new gcc compilers need the std= directive

CXXFLAGS = -O2 -I . -I ../include -Wall -ffast-math -DLINUX -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 
# Auto check gcc version from Raj Guha
GCCVERSIONGTEQ4 := $(shell expr `gcc -dumpversion | cut -f1,2 -d.` \>= 4.3)
ifeq "$(GCCVERSIONGTEQ4)" "1"
    CXXFLAGS += -std=gnu++0x
endif


AR = ar
CP = cp
RM = rm


CFLAGS = -O2 -I . -I ../include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

LIBRARY_OBJECTS = IWString_class.o cmdline.o iwstring_data_source.o logical_expression.o msi_object.o new_int.o\
	iwzlib.o iwgrep-2.5.o grep-2.5.regex.o iwstring_and_file_descriptor.o iwrandom.o mtrand.o iwwrite.o iwstring.o iwbits.o\
	du_bin2ascii.o bits_in_common.o iwstrncasecmp.o int_comparator.o write_space_suppressed_string.o\
	dash_f.o iw_stl_hash_map.o  iwwrite_block.o KahanSum.o

all:libiwsupport.a

libiwsupport.a: $(LIBRARY_OBJECTS)
	$(AR) cr $@ $(LIBRARY_OBJECTS)

install: libiwsupport.a
	$(CP) $< ../lib

clean:
	-$(RM) $(LIBRARY_OBJECTS) libiwsupport.a

uninstall:
	-$(RM) ../lib/libiwsupport.a
