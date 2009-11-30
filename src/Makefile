CC = g++
CXXFLAGS = -I. -I../minmaxer/source.triangulation/ `pkg-config --cflags gtkglext-1.0 libglade-2.0`
LDFLAGS = -L../minmaxer/minmaxer.bin
GTKLIBS = `pkg-config --libs gtkglext-1.0 libglade-2.0`
OSGLIBS = -losg -losgViewer
LDLIBS = -lGL -lGLU -lglut $(GTKLIBS) $(OSGLIBS)

all: triangulate vis vis1 vis2

triangulate: triangulate.o ../minmaxer/source.triangulation/libtriangulation.a

triangulate.o: triangulate.cc

vis: vis.o

vis.o: vis.cc

vis1: vis1.o

vis1.o: vis1.cc

vis2: vis2.o

vis2.o: vis2.cc
