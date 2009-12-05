import numpy as np
import gts

class mystruct:
  pass

class Point:
  def __init__(self, x=0, y=0, z=0):
    self.x = x
    self.y = y
    self.z = z

  def set(self, x, y, z = None):
    self.x = x
    self.y = y
    if not z is None:
      self.z = z

  def __repr__(self):
    return "(%s %s %s)" % (self.x, self.y, self.z)

# Read vertices data from file
def read_vertices(vertices, filename = "hw4.vertices"):
  f = open(filename)
  for line in f:
    s = [float(x) for x in line.split()]
    p = gts.Vertex(s[0], s[1], s[2])
    vertices.append(p)
  print "Read %d points" % len(vertices)
  f.close()

# Read tour sights data from file
def read_tour(points, filename = "hw4.tour"):
  print "Reading tour sites..."
  return read_vertices(points, filename)

# Read triangle info from file
def read_triangles(triangles, filename = "hw4.faces"):
  f = open(filename)
  for line in f:
    s = [int(x) for x in line.split()]
    s.sort()
    triangles.append(s)
  print "Read %d triangles" % len(triangles)
  f.close()

# Make a gts.Triangle array from vertices and triangle indices
def make_gts_triangles(gts_vertices, face_indices, gts_triangles):
  for f in face_indices:
    triangle = gts.Triangle(gts_vertices[f[0]],
                            gts_vertices[f[1]],
                            gts_vertices[f[2]])
    if (triangle.orientation() < 0):
      triangle.revert()
    gts_triangles.append(triangle)

# get span of a set of points
def get_extents(points):
  pmin = None
  pmax = None
  for p in points:
    if pmin is None:
      pmin = Point(p.x, p.y, p.z)
    if pmax is None:
      pmax = Point(p.x, p.y, p.z)
    pmin.x = min(pmin.x, p.x)
    pmin.y = min(pmin.y, p.y)
    pmin.z = min(pmin.z, p.z)
    pmax.x = max(pmax.x, p.x)
    pmax.y = max(pmax.y, p.y)
    pmax.z = max(pmax.z, p.z)

  dx = pmax.x - pmin.x
  dy = pmax.y - pmin.y
  dz = pmax.z - pmin.z
  span = Point(dx, dy, dz)
  return (pmin, pmax, span)

class HeightLookup:
  def __init__(self, vertices, faces):
    self.V = vertices
    self.F = faces
    self.n = 200  # divide into nxn buckets
    self.hmap = None
    (self.vmin, self.vmax, self.vspan) = get_extents(self.V)
    self.create_map()

  def create_map(self):
    self.hmap = [[[] for a in range(self.n)] for b in range(self.n)]

    for f in self.F:
      (fmin, fmax, fspan) = get_extents(f.vertices())
      bx1 = (fmin.x - self.vmin.x)/self.vspan.x * (self.n - 1)
      bx2 = (fmax.x - self.vmin.x)/self.vspan.x * (self.n - 1)
      by1 = (fmin.y - self.vmin.y)/self.vspan.y * (self.n - 1)
      by2 = (fmax.y - self.vmin.y)/self.vspan.y * (self.n - 1)
      for bx in range(int(bx1), int(bx2) + 1):
        for by in range(int(by1), int(by2) + 1):
          self.hmap[bx][by].append(f)

  def get_height(self, point, deb = False):
    bx = (point.x - self.vmin.x)/self.vspan.x * (self.n - 1)
    by = (point.y - self.vmin.y)/self.vspan.y * (self.n - 1)
    for f in self.hmap[int(bx)][int(by)]:
      if deb: print [v.coords() for v in f.vertices()]
      if point.is_in(f) >= 0:
        return f.interpolate_height(point)
    return None


def get_height(point, faces):
  for f in faces:
    if point.is_in(f) >= 0:
      print [v.coords() for v in f.vertices()]
      return f.interpolate_height(point)
  return None


def get_vector_2d(p1, p2):
  x = p2.x - p1.x
  y = p2.y - p1.y
  l = math.sqrt(x*x + y*y)
  return Point(x/l, y/l)

# returns a vector which is average of
# vector from P1 to P2 and P2 to P3 in 2D
def get_tangent(p1, p2, p3):
  v1 = get_vector_2d(p1, p2)
  v2 = get_vector_2d(p2, p3)
  x = (v1.x + v2.x) / 2.0
  y = (v1.y + v2.y) / 2.0
  l = math.sqrt(x*x + y*y)
  return Point(x/l, y/l)


# (v.x, v.y): tangent vector
# l: length to interpolate along tangent vector
# p1, p2: new points
def insert_points(p, v, l, p1, p2):
  p1.x = p.x - v.x * l;
  p1.y = p.y - v.y * l;
  p2.x = p.x + v.x * l;
  p2.y = p.y + v.y * l;


def interpolate_2d(p1, p2, l):
  p = Point()
  v = get_vector_2d(p1, p2)
  p.x = p1.x + l*v.x
  p.y = p1.y + l*v.y
  return p

def mid_point(p1, p2):
  return Point((p1.x + p2.x)/2, (p1.y + p2.y)/2, (p1.z + p2.z)/2)

def bez(p0, p1, p2, t):
  p = Point()
  p.x = (1-t)*(1-t)*p0.x + 2*(1-t)*t*p1.x + t*t*p2.x;
  p.y = (1-t)*(1-t)*p0.y + 2*(1-t)*t*p1.y + t*t*p2.y;
  p.z = (1-t)*(1-t)*p0.z + 2*(1-t)*t*p1.z + t*t*p2.z;
  return p

def dist(p1, p2):
  sqdist = (p1.x - p2.x)**2 + (p1.y - p2.y)**2 + (p1.z - p2.z)**2
  return math.sqrt(sqdist)

# sample bezier between sample points s1 and s2
def sample_bezier(p0, p1, p2, s1, s2, min_dist):
  t = (s1.t + s2.t) / 2.0
  s = bez(p0, p1, p2, t)
  s.t = t

  d1 = dist(s1, s2)
  d2 = dist(s1, s)
  d3 = dist(s, s2)
  if max(d1, d2, d3) < min_dist:
    return []

  V1 = sample_bezier(p0, p1, p2, s1, s, min_dist)
  V2 = sample_bezier(p0, p1, p2, s, s2, min_dist)
  return V1 + [s] + V2



def get_max_height_line(hmap, p0, p1):
  d = dist(p0, p1)
  mh = hmap.get_height(gts.Point(p0.x, p0.y))
  for i in range(0, 201):
    p = interpolate_2d(p0, p1, i*d/200)
    h = hmap.get_height(gts.Point(p.x, p.y))
    mh = max(h, mh)
  return mh


def rotate_point_2d(p, theta):
  np = Point()
  np.x = p.x * math.cos(theta) - p.y * math.sin(theta)
  np.y = p.x * math.sin(theta) + p.y * math.cos(theta)
  return np
