ts_x = []
ts_y = []
ts_z = []
ts = []
ml = 1000

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

minView = Point(-23469, -19540, -2)
maxView = Point(21859, 19774, 1282)

# read tour sites
f = open("hw4.tour")
for line in f:
  s = [float(x) for x in line.split()]
  p = Point(s[0], s[1], s[2])
  ts.append(p)
  ts_x.append(p.x)
  ts_y.append(p.y)
  ts_z.append(p.z)

nts = len(ts)
print "No. of sights:", nts, ts[0]

def get_vector(p1, p2):
  x = p2.x - p1.x
  y = p2.y - p1.y
  l = sqrt(x*x + y*y)
  return Point(x/l, y/l)

# returns a vector which is average of
# vector from P1 to P2 and P2 to P3
def get_tangent(p1, p2, p3):
  v1 = get_vector(p1, p2)
  v2 = get_vector(p2, p3)
  x = (v1.x + v2.x) / 2.0
  y = (v1.y + v2.y) / 2.0
  l = sqrt(x*x + y*y)
  return Point(x/l, y/l)

tang = get_tangent(Point(0, 0), Point(0, 5), Point(5, 5))
print tang

# (v.x, v.y): tangent vector
# l: length to interpolate along tangent vector
# p1, p2: new points
def insert_points(p, v, l, p1, p2):
  p1.x = p.x - v.x * l;
  p1.y = p.y - v.y * l;
  p2.x = p.x + v.x * l;
  p2.y = p.y + v.y * l;

p1 = Point()
p2 = Point()
insert_points(Point(0, 5), tang, 1, p1, p2)
print p1, p2


def interpolate(p1, p2, l):
  p = Point()
  v = get_vector(p1, p2)
  p.x = p1.x + l*v.x
  p.y = p1.y + l*v.y
  return p

print interpolate(Point(0, 0), Point(0, 5), 2)

np = []  # new points to interpolate
cp = []  # bezier control points


plot_s = []  # original sites
plot_i = []  # inserted interpolated
plot_c = []  # control points

np.append(ts[0])
plot_s.append(ts[0])

first_cp = Point()

for i in range(1, len(ts) - 1):
  p1 = ts[i - 1]
  p2 = ts[i]
  p3 = ts[i + 1]
  tang = get_tangent(p1, p2, p3)

  np1 = Point()
  np2 = Point()
  insert_points(p2, tang, ml, np1, np2)
  p2.prev = np1
  p2.next = np2

ts[0].next = Point(ts[0].x, ts[0].y)
ts[-1].prev = Point(ts[-1].x, ts[-1].y)

for i in range (0, len(ts) - 1):
  p1 = ts[i]
  p2 = ts[i + 1]
  tang = get_vector(p1.next, p2.prev)

  p3 = interpolate(p1.next, p2.prev, ml)
  p5 = interpolate(p2.prev, p1.next, ml)
  p4 = Point((p3.x + p5.x)/2, (p3.y + p5.y)/2)

  np.append(p1)
  np.append(p3)
  np.append(p4)
  np.append(p5)

  plot_s.append(p1)
  plot_i.append(p3)
  plot_i.append(p4)
  plot_i.append(p5)
  plot_c.append(p1.next)
  plot_c.append(p2.prev)


header = """<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">

<svg width="1000px" height="1000px" viewBox="-25000 -25000 50000 50000"
     xmlns="http://www.w3.org/2000/svg" version="1.1">

  <title>Example quad01 - quadratic Bezier commands in path data</title>

  <rect x="-23469" y="-19540" width="45328" height="39314"
        fill="none" stroke="blue" stroke-width="15" />
"""

footer = '</svg>'


index = 1
s = []
s.append('<path d="M%d, %d' % (np[0].x, np[0].y))
# calculate first control point as midpoint of first two points
cx = (np[0].x + np[1].x) / 2.0
cy = (np[0].y + np[1].y) / 2.0
s.append('Q%d,%d %d,%d' % (cx, cy, np[1].x, np[1].y))
for i in range(2, len(np)):
  s.append("T%d,%d" % (np[i].x, np[i].y))
s.append('" \n fill="none" stroke="red" stroke-width="40"  /> ')
curve =  " ".join(s)


def svg_points(points, color, width=130):
  s = []
  s.append('<g fill="%s">' % color)
  for p in points:
    s.append('<circle cx="%d" cy="%d" r="%s"/>' % (p.x, p.y, width))
  s.append('</g>')
  return "\n".join(s)

f = open("t2.svg", "w")
print >> f, header
print >> f, curve
print >> f, svg_points(plot_s, "black")
print >> f, svg_points(plot_c, "#9f9", 100)
print >> f, svg_points(plot_i, "gray")
print >> f, footer
f.close()
