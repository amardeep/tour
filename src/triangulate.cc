#include "triangulation.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
int GLOBAL_INTERRUPT_ALGORITHM = 0;
extern "C" {
  void copyCoordinatesToGraph(int nofSites, int *x, int *y, int *z,
                              int eliminateDuplicates, char **gExternal);
  void copyGraphToListOfTriangles (char *gExternal, triangleList **tl);


  void planeSweep(char* gExternal);
  void delaunay1(char* gExternal);
}

void readSites(string filename, int **x, int **y, int **z, int *numSites) {
  ifstream fin(filename.c_str());

  string ignore;
  int rx, ry, rz;
  vector<int> vx, vy, vz;
  while (fin >> ignore >> rx >> ry >> rz) {
    vx.push_back(rx);
    vy.push_back(ry);
    vz.push_back(rz);
  }

  *numSites = vx.size();

  *x = new int[*numSites];
  *y = new int[*numSites];
  *z = new int[*numSites];

  for (int i = 0; i < *numSites; ++i) {
    (*x)[i] = vx[i];
    (*y)[i] = vy[i];
    (*z)[i] = vz[i];
  }
}


void writeSites(string filename, int *x, int *y, int *z, int numSites) {
  ofstream fout(filename.c_str());

  for (int i = 0; i < numSites; ++i) {
    fout << x[i] << " " << y[i] << " " << z[i] << endl;
  }
  fout.close();
}

void writeTriangles(string filename, triangleList *list) {
  ofstream fout(filename.c_str());

  for (int i = 0; i < list->nofTriangles; ++i) {
    fout << list->v[i][0] << " " << list->v[i][1] << " "
         << list->v[i][2] << endl;
  }
  fout.close();
}


int main() {
  int *x, *y, *z, n;
  char *g;
  triangleList *list;

  cout << "Reading file hw4.heights" << endl;
  readSites("hw4.heights", &x, &y, &z, &n);
  cout << "No. of sites = " << n << endl;

  copyCoordinatesToGraph (n, x, y, z, 1, &g);
  planeSweep(g);
  delaunay1(g);

  copyGraphToListOfTriangles(g, &list);
  cout << "num triangles: " << list->nofTriangles << endl;
  cout << "max triangles: " << list->maxTriangles << endl;

  writeSites("hw4.vertices", x, y, z, n);
  writeTriangles("hw4.edges", list);
}
