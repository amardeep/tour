#include <osgViewer/Viewer>
#include <osgText/Text>
#include <osg/Geometry>
#include <osg/ShapeDrawable>

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;
using namespace osg;

using osg::Geometry;
using osg::Vec3;
using osg::Vec4;
using osg::Vec3Array;
using osg::Vec4Array;
using osg::DrawElementsUInt;
using osg::PrimitiveSet;
using osg::Geode;
using osg::StateSet;

void ReadTour(Vec3Array* v) {
  ifstream fin("hw4.tour");
  double x, y, z;
  while (fin >> x >> y >> z) {
    v->push_back(Vec3(x, y, z));
  }
  cout << "Num tour vertices = " << v->size() << endl;
}

void ReadVertices(Vec3Array* v, Vec3 *minv, Vec3 *maxv) {
  ifstream fin("hw4.vertices");
  float x, y, z;
  int index = 0;
  while (fin >> x >> y >> z) {
    v->push_back(Vec3(x, y, z));
    if (index == 0) {
      minv->set(x, y, z);
      maxv->set(x, y, z);
    } else {
      minv->set(min(x, minv->x()),
                min(y, minv->y()),
                min(z, minv->z()));
      maxv->set(max(x, maxv->x()),
                max(y, maxv->y()),
                max(z, maxv->z()));
    }
    index++;
  }
  cout << "min: "
       << (*minv)[0] << " "
       << (*minv)[1] << " "
       << (*minv)[2] << endl;
  cout << "max: "
       << (*maxv)[0] << " "
       << (*maxv)[1] << " "
       << (*maxv)[2] << endl;
}

void ReadTriangles(Geometry* geo) {
  ifstream fin("hw4.edges");
  int v1, v2, v3;
  while (fin >> v1 >> v2 >> v3) {
    DrawElementsUInt* t = new DrawElementsUInt(PrimitiveSet::TRIANGLES, 0);
    t->push_back(v1);
    t->push_back(v2);
    t->push_back(v3);
    geo->addPrimitiveSet(t);
  }
}

osg::Node* CreateScene() {
  // create drawable geometry object
  Geometry* geo = new Geometry;

  // add 3 vertices creating a triangle
  Vec3Array* v = new Vec3Array;
  Vec3 minv, maxv;
  ReadVertices(v, &minv, &maxv);
  geo->setVertexArray(v);

  // create a primitive set (add index numbers)
  ReadTriangles(geo);

  // create color array data for each corner of our triangle
  Vec4Array* colors = new Vec4Array;
  for (int i = 0; i < v->size(); ++i) {
    colors->push_back(Vec4(drand48(), drand48(), drand48(), 1.0 ) );
  }
  geo->setColorArray(colors);

  // make sure that our geometry is using one color per vertex
  geo->setColorBinding(Geometry::BIND_PER_VERTEX);

  // create geometry node that will contain all our drawables
  Geode* geode = new Geode;
  StateSet* stateSet = geode->getOrCreateStateSet();
  stateSet->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
  geode->addDrawable(geo);

  Vec3Array* tour = new Vec3Array;
  ReadTour(tour);

  for (int i = 0; i < tour->size(); ++i) {
    geode->addDrawable(new ShapeDrawable(new Box((*tour)[i], 500.0)));
  }

  return geode;
}


osg::Node* CreateTour() {
  Geometry* geo = new Geometry;

  Vec3Array* v = new Vec3Array;
  ReadTour(v);
  geo->setVertexArray(v);


  // create color array data for each corner of our triangle
  Vec4Array* colors = new Vec4Array;
  for (int i = 0; i < v->size(); ++i) {
    colors->push_back(Vec4(drand48(), drand48(), drand48(), 1.0 ) );
  }
  geo->setColorArray(colors);

  // make sure that our geometry is using one color per vertex
  geo->setColorBinding(Geometry::BIND_PER_VERTEX);

  // create geometry node that will contain all our drawables
  Geode* geode = new Geode;
  StateSet* stateSet = geode->getOrCreateStateSet();
  stateSet->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
  geode->addDrawable(geo);

  return geode;
}


int main(int argc, char **argv) {
  osg::ref_ptr<osgViewer::Viewer> viewer = new osgViewer::Viewer;

  // make the viewer create a 512x512 window and position it at 32, 32
  viewer->setUpViewInWindow(32, 32, 512, 512);

  // set the scene-graph data the viewer will render
  viewer->setSceneData(CreateScene());

  // execute main loop
  return viewer->run();
}
