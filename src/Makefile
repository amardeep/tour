CC = g++
CXXFLAGS = -I. -I../minmaxer/source.triangulation/ `pkg-config --cflags gtkglext-1.0 libglade-2.0` -I../OpenSceneGraph-2.8.2/include `python2.5-config --includes`
LDFLAGS = -L../minmaxer/minmaxer.bin -L../OpenSceneGraph-2.8.2/lib
GTKLIBS = `pkg-config --libs gtkglext-1.0 libglade-2.0`
OSGLIBS = -losg -losgViewer -losgSim
PYLIBS = `python2.5-config --libs`
LDLIBS = -lGL -lGLU -lglut $(GTKLIBS) $(OSGLIBS) $(PYLIBS)

all: triangulate vis vis1 vis2 vis3 vis_points gui

triangulate: triangulate.o ../minmaxer/source.triangulation/libtriangulation.a

triangulate.o: triangulate.cc

vis: vis.o

vis.o: vis.cc

vis1: vis1.o

vis1.o: vis1.cc

vis2: vis2.o

vis2.o: vis2.cc

vis3: vis3.o

vis3.o: vis3.cc

vis_points: vis_points.cc

gui: gui.o

gui.o: gui.cc