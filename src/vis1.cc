#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <glade/glade.h>

#include <GL/gl.h>
#include <GL/glu.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

GLuint theTorus;
GLuint theSphere;
GLuint theVertices;
GLuint theTriangles;

GLint *vertices;
int num_vertices;

int minx, miny, minz;
int maxx, maxy, maxz;

GLint *triangles;
int num_triangles;

/* Draw a torus */
void torus(int numc, int numt) {
  int i, j, k;
  double s, t, x, y, z, twopi;

  twopi = 2 * (double)M_PI;
  for (i = 0; i < numc; i++) {
    glBegin(GL_QUAD_STRIP);
    for (j = 0; j <= numt; j++) {
      for (k = 1; k >= 0; k--) {
        s = (i + k) % numc + 0.5;
        t = j % numt;

        x = (1+.1*cos(s*twopi/numc))*cos(t*twopi/numt);
        y = (1+.1*cos(s*twopi/numc))*sin(t*twopi/numt);
        z = .1 * sin(s * twopi / numc);
        glVertex3f(x, y, z);
      }
    }
    glEnd();
  }
}

/* Create display list with Torus and initialize state*/
void init(void) {
  theTorus = glGenLists(1);
  glNewList(theTorus, GL_COMPILE);
  torus(8, 25);
  glEndList();

  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(3, GL_INT, 0, vertices);

  theVertices = glGenLists(1);
  glNewList(theVertices, GL_COMPILE);
  glPointSize(3);
  glDrawArrays(GL_POINTS, 0, num_vertices);
  glEndList();

  theTriangles = glGenLists(1);
  glNewList(theTriangles, GL_COMPILE);
  glDrawElements(GL_TRIANGLES, num_triangles*3, GL_UNSIGNED_INT, triangles);
  glEndList();

}

void read_vertices() {
  vector<int> x, y, z;
  int sx, sy, sz;

  ifstream fin("hw4.vertices");
  while (fin >> sx >> sy >> sz) {
    x.push_back(sx);
    y.push_back(sy);
    z.push_back(sz);
  }
  cout << "Num vertices = " << x.size() << endl;
  num_vertices = x.size();

  minx = maxx = x[0];
  miny = maxy = y[0];
  minz = maxz = z[0];

  vertices = new GLint[x.size() * 3];
  for (int i = 0; i < x.size(); ++i) {
    vertices[i*3] = x[i];
    vertices[i*3 + 1] = y[i];
    vertices[i*3 + 2] = z[i];
    minx = min(minx, x[i]);
    miny = min(miny, y[i]);
    minz = min(minz, z[i]);
    maxx = max(maxx, x[i]);
    maxy = max(maxy, y[i]);
    maxz = max(maxz, z[i]);
  }

  cout << "minmax "
       << minx << " " << maxx << " - "
       << miny << " " << maxy << " - "
       << minz << " " << maxz << endl;
}

void read_edges() {
  vector<int> v;

  ifstream fin("hw4.edges");
  int v1, v2, v3;

  while (fin >> v1 >> v2 >> v3) {
    v.push_back(v1);
    v.push_back(v2);
    v.push_back(v3);
  }
  cout << "Num triangles = " << v.size()/3 << endl;
  num_triangles = v.size() / 3;

  triangles = new GLint[v.size()];
  for (int i = 0; i < v.size(); ++i) {
    triangles[i] = v[i];
  }
}

extern "C" {

gboolean on_window1_delete_event(GtkWidget *widget, GdkEvent *event,
                                 gpointer user_data) {
  gtk_main_quit ();
  return FALSE;
}

void on_drawingarea1_realize(GtkWidget *widget, gpointer user_data) {
  GdkGLContext *glcontext = gtk_widget_get_gl_context(widget);
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);

  GLUquadricObj *qobj;
  static GLfloat light_diffuse[] = {1.0, 0.0, 0.0, 1.0};
  static GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};

  if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return;

  qobj = gluNewQuadric ();
  gluQuadricDrawStyle (qobj, GLU_FILL);
  theSphere = glGenLists(1);
  glNewList (theSphere, GL_COMPILE);
  gluSphere (qobj, 1, 20, 20);
  glEndList ();

  init();

  glLightfv (GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv (GL_LIGHT0, GL_POSITION, light_position);
  glEnable (GL_LIGHTING);
  glEnable (GL_LIGHT0);
  glEnable (GL_DEPTH_TEST);

  glClearColor (1.0, 1.0, 1.0, 1.0);
  glClearDepth (1.0);

  glViewport (0, 0,
              widget->allocation.width, widget->allocation.height);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluPerspective (120.0, 1.0, 1.0, 10000.0);

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
  GLdouble eyez = maxz + (maxz - minz);
  GLdouble eyex = (minx + maxx) / 2.0;
  GLdouble eyey = (miny + maxy) / 2.0;

  gluLookAt (0, 0, 5,
             0, 0, 0,
             0.0, 1.0, 0.0);
  glTranslatef (0.0, 0.0, -5000.0);

  gdk_gl_drawable_gl_end (gldrawable);
}


gboolean on_drawingarea1_configure_event(GtkWidget *widget,
                                         GdkEventConfigure *event,
                                         gpointer user_data) {
  GdkGLContext *glcontext = gtk_widget_get_gl_context(widget);
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);

  if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return FALSE;
  glViewport (0, 0, widget->allocation.width, widget->allocation.height);
  gdk_gl_drawable_gl_end (gldrawable);

  return FALSE;
}


gboolean on_drawingarea1_expose_event(GtkWidget *widget,
                                      GdkEventExpose *event,
                                      gpointer user_data) {
  GdkGLContext *glcontext = gtk_widget_get_gl_context(widget);
  GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(widget);

  if (!gdk_gl_drawable_gl_begin (gldrawable, glcontext)) return FALSE;
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glCallList(theSphere);
  //glCallList(theTorus);
  glCallList(theVertices);
  glCallList(theTriangles);
  gdk_gl_drawable_swap_buffers(gldrawable);
  gdk_gl_drawable_gl_end (gldrawable);

  return FALSE;
}


void on_button1_clicked(GtkButton *button, gpointer user_data) {
  gtk_main_quit ();
}

}  // end extern "C"


int main (int argc, char **argv) {
  GdkGLConfig *glconfig;
  GladeXML *xml;
  GtkWidget *window;
  GtkWidget *drawingarea;

  read_vertices();
  read_edges();

  gtk_init (&argc, &argv);
  gtk_gl_init (&argc, &argv);

  /* Try double-buffered visual */
  int mode = GDK_GL_MODE_RGB | GDK_GL_MODE_DEPTH | GDK_GL_MODE_DOUBLE;
  glconfig = gdk_gl_config_new_by_mode(GdkGLConfigMode(mode));

  if (glconfig == NULL) {
    cerr << "Cannot find the double-buffered visual." << endl;
    cerr << "Trying single-buffered visual" << endl;
    return 1;
  }

  xml = glade_xml_new("simple.glade", NULL, NULL);
  glade_xml_signal_autoconnect(xml);
  window = glade_xml_get_widget (xml, "window1");
  gtk_container_set_reallocate_redraws(GTK_CONTAINER(window), TRUE);

  drawingarea = glade_xml_get_widget(xml, "drawingarea1");
  /* Add OpenGL-capability to drawingarea1. */
  gtk_widget_set_gl_capability(drawingarea,
                               glconfig,
                               NULL,
                               TRUE,
                               GDK_GL_RGBA_TYPE);
  gtk_widget_show (window);
  gtk_main ();
  return 0;
}
