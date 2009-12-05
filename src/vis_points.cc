#include <osgViewer/Viewer>
#include <osgText/Text>
#include <osg/Geometry>
#include <osg/ShapeDrawable>

#include <osgDB/ReadFile>
#include <osg/io_utils>
#include <osgUtil/IntersectionVisitor>
#include <osgUtil/LineSegmentIntersector>

#include <osgSim/LineOfSight>
#include <osgSim/HeightAboveTerrain>
#include <osgSim/ElevationSlice>

#include <osgUtil/Optimizer>
#include <osgUtil/SmoothingVisitor>

#include <osgViewer/CompositeViewer>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#define STEP 100
#define DRAWSTEP 15

using namespace std;
using namespace osg;


// Return a transition value between 0 & 1 following a cubic
// spline to accentuate changes in beginning and end, and flatten
// the middle.
float getTransitionValue(float t) {
  if (t < 0) return 0;
  if (t > 1) return 1;

  float a = .8;
  float b = 1 - a;

  float val =  3*(1-t)*(1-t)*t*a +
               3*(1-t)*t*t*b +
               t*t*t;
  return val;
}


void ReadPoints(Vec3Array* v, string filename) {
  ifstream fin(filename.c_str());
  double x, y, z;
  while (fin >> x >> y >> z) {
    v->push_back(Vec3(x, y, z));
  }
  cout << "Num vertices from " << filename << " = " << v->size() << endl;
}


void ReadTour(Vec3Array* v) {
  return ReadPoints(v, "hw4.tour");
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


// Reads bezier control points from file
void ReadControl(Vec3Array* v) {
  ifstream fin("tour3.txt");
  double x, y, z;
  char buffer[100];
  fin >> buffer;
  fin >> buffer;
  fin >> buffer;

  while (fin >> x >> y) {
    v->push_back(Vec3(x, y, 0));
  }
  cout << "Num control vertices = " << v->size() << endl;
}


// Return height of terrain at a poing
Vec3d getHeight(osg::Node* scene, double x, double y){
  BoundingSphere bs = scene->getBound();
  Vec3d start = Vec3d(x, y, 8 * bs.radius());
  Vec3d end = Vec3d(x, y, -8 * bs.radius());

  using osgSim::LineOfSight;
  osgSim::LineOfSight los;
  los.addLOS(start, end);;
  los.computeIntersections(scene);
  const LineOfSight::Intersections& intersections = los.getIntersections(0);
  for (LineOfSight::Intersections::const_iterator itr = intersections.begin();
       itr != intersections.end(); ++itr) {
      return *itr;
  }
}


double L2Norm(Vec3 x){
 return (sqrt(pow(x.x(), 2) +
              pow(x.y(), 2) +
              pow(x.z(), 2)));
}


bool isLinear(Vec3 x, Vec3 y, Vec3 z){
  Vec3 slope1, slope2;
  double diff;
  slope1= y - x,
  slope2= z - y;
  slope1= slope1 / (L2Norm(slope1));
  slope2= slope2 / (L2Norm(slope2));
  diff=L2Norm(slope1-slope2);
  return (diff < .000001);
}

//Returns a Vec3 with (slope1, slope2, height), Slope 1 is the max slope 
//from control point P1 to a terrain point, slope2 is the max slope from
//P3, height is the maximum height along the segment.  
Vec3 getMaxSlopeHeight(Vec3 P1, Vec3 P2, Vec3 P3, osg::Node* terrain){
  Vec3 slopeHeight;
  slopeHeight.set(0.0,0.0,0.0);
  double slope1=0;
  double slope2=0;
  double height=0;
  Vec3 point;
  double sample, slope;
  double t;
  for(int j=1; j<(STEP); j++){
      t=(double)j/STEP;
      point = P1*((1-t)*(1-t)) + P2 * (2*(1-t)*t) + P3*(t*t);
      
      //get Height of current point
      sample = getHeight(terrain, point.x(), point.y()).z();
      if(sample>height)
	height=sample;
	
      //Calculate slope from P1
      slope=(sample-P1.z())/t;
      if(slope>slope1)
	slope1=slope;
      
       //Calculate slope from P3
      slope=(sample-P3.z())/(1-t);
      if(slope>slope2)
	slope2=slope;
  }
  slopeHeight.set(slope1, slope2, height);
  return slopeHeight;
}

osg::Node* CreateScene() {
  Geometry* geo = new Geometry;

  Vec3Array* v = new Vec3Array;
  Vec3 minv, maxv;
  ReadVertices(v, &minv, &maxv);
  geo->setVertexArray(v);

  // create a primitive set (add index numbers)
  ReadTriangles(geo);

  // create color array data for each corner of our triangle
  Vec4Array* colors = new Vec4Array;
  float span_z = maxv.z()- minv.z();
  for (int i = 0; i < v->size(); ++i) {
    float t = ((*v)[i].z() - minv.z()) / span_z;
    float color = getTransitionValue(t);
    colors->push_back(Vec4(.1, color, .5, 1.0 ));
  }
  geo->setColorArray(colors);

  // make sure that our geometry is using one color per vertex
  geo->setColorBinding(Geometry::BIND_PER_VERTEX);

  // create geometry node that will contain all our drawables
  Geode* geode = new Geode;
  StateSet* stateSet = geode->getOrCreateStateSet();
  stateSet->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
  //osgUtil::SmoothingVisitor::smooth(*geo);
  geode->addDrawable(geo);

  Vec3Array* tour = new Vec3Array;
  ReadTour(tour);
  for (int i = 0; i < tour->size(); ++i) {
    //geode->addDrawable(new ShapeDrawable(new Box((*tour)[i], 150.0)));
  }

  Vec3Array* points = new Vec3Array;
  ReadPoints(points, "points.txt");
  for (int i = 0; i < points->size(); ++i) {
    ShapeDrawable* s = new ShapeDrawable(new Sphere((*points)[i], 30.0));
    s->setColor(Vec4(1, .2, .2, 1));
    geode->addDrawable(s);
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
  ref_ptr<osgViewer::Viewer> viewer = new osgViewer::Viewer;
  viewer->setUpViewInWindow(32, 32, 800, 800);
  osg::ref_ptr<osg::Node> scene = CreateScene();
  viewer->setSceneData(scene);
  return viewer->run();
}
