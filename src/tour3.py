import math
import util

ts_x = []
ts_y = []
ts_z = []
ts = []
ml = 1000


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


# reflect point p around origin o
def reflect(p, o):
  n = Point()
  n.x = 2*o.x - p.x
  n.y = 2*o.y - p.y
  return n

cp = []  # de-boor control points

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

  if i == 0:
    cp1 = reflect(p2.prev, p5)
    cp2 = reflect(cp1, p1)
    cp.append(cp2)
    cp.append(cp1)
  elif i < len(ts) - 2:
    cp.append(p1.prev)
    cp.append(p1.next)
    cp.append(reflect(p1.next, p3))
    cp.append(reflect(p2.prev, p5))
  else:
    cp.append(p1.prev)
    cp.append(p1.next)
    cp1 = reflect(p1.next, p3)
    cp2 = reflect(cp1, p2)
    cp.append(cp1)
    cp.append(cp2)


# arc length of a quadratic bezier curve
# http://segfaultlabs.com/graphics/qbezierlen/
def arc_length(p0, p1, p2):
  a = Point()
  b = Point()
  a.x = p0.x - 2*p1.x + p2.x
  a.y = p0.y - 2*p1.y + p2.y
  b.x = 2*p1.x - 2*p0.x
  b.y = 2*p1.y - 2*p0.y

  A = 4*(a.x*a.x + a.y*a.y)
  B = 4*(a.x*b.x + a.y*b.y)
  C = b.x*b.x + b.y*b.y

  Sabc = 2 * sqrt(A+B+C)
  A_2 = sqrt(A)
  A_32 = 2 * A * A_2
  C_2 = 2 * sqrt(C)
  BA = B/A_2

  logterm = log((2*A_2+BA+Sabc)/(BA+C_2))
  if logterm != logterm:
    print "Got a nan:", p0, p1, p2
    logterm = 0

  val = (A_32*Sabc + A_2*B*(Sabc-C_2) + (4*C*A-B*B)*logterm) / (4*A_32);
  return val

print arc_length(Point(0, 0), Point(12, 0), Point(20, 0))
print arc_length(Point(0, 0), Point(2, 0), Point(20, 0))


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


finfo = open("tour3.txt", "w")


s = []
print >> finfo, "#Bezier Control points"
x = (cp[0].x + cp[1].x) / 2.0
y = (cp[0].y + cp[1].y) / 2.0
s.append('<path d="M%d, %d' % (x, y))
print >> finfo, x, y
for i in range(1, len(cp) - 1):
  x = (cp[i].x + cp[i+1].x) / 2.0
  y = (cp[i].y + cp[i+1].y) / 2.0
  s.append('Q%d,%d %d,%d' % (cp[i].x, cp[i].y, x, y))
  print >> finfo, cp[i].x, cp[i].y
  print >> finfo, x, y
s.append('" \n fill="none" stroke="red" stroke-width="40"  /> ')
curve = " ".join(s)


s = []
print >> finfo, "\n\n#De-boor control points"
s.append('<path d="')
for i in range(0, len(cp)):
  cmd = "L"
  if i == 0: cmd = "M"
  s.append("%s %s %s" % (cmd, cp[i].x, cp[i].y))
  print >> finfo, cp[i].x, cp[i].y
s.append('" \n fill="none" stroke="#ccc" stroke-width="80"  /> ')
db = " ".join(s)

def svg_points(points, color, width=130):
  s = []
  s.append('<g fill="%s">' % color)
  for p in points:
    s.append('<circle cx="%d" cy="%d" r="%s"/>' % (p.x, p.y, width))
  s.append('</g>')
  return "\n".join(s)

f = open("t3.svg", "w")
print >> f, header
print >> f, db
print >> f, curve
print >> f, svg_points(ts, "black")
print >> f, svg_points(cp, "#9f9", 100)
print >> f, footer
f.close()
finfo.close()
