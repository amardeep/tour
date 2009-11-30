from numpy import *
from scipy import *
from scipy import linalg
import matplotlib.pyplot as plt

ts = []

minView = [-23469, -19540, -2]
maxView = [21859, 19774, 1282]


ml = 1000

# read tour sites
f = open("hw4.tour")
tourSites = []
for line in f:
  s = [float(x) for x in line.split()]
  tourSites.append(s)
ts = mat(tourSites)

p = []
cp = []
x1 = ts[0, 0]
y1 = ts[0, 1]
cp.append((0, 0))
p.append((x1, y1))

x3 = ts[1, 0]
y3 = ts[1, 1]
x2 = (x1 + x3) / 2
y2 = (y1 + y3) / 2
cp.append((x2, y2))
p.append((x3, y3))

for i in range(1, len(ts) - 1):
  l = sqrt((x3-x1)**2 + (y3-y1)**2)
  x2 = x3 + ml*(x3-x1)/l
  y2 = y3 + ml*(y3-y1)/l
  x1 = x3
  y1 = y3
  nx = ts[i + 1, 0]
  ny = ts[i + 1, 1]
  l = sqrt((x2-nx)**2 + (y2-ny)**2)
  x3 = x2 - ml*(x2-nx)/l
  y3 = y2 - ml*(y2-ny)/l
  cp.append((x2, y2))
  p.append((x3, y3))

  x2 = (x3 + nx) / 2
  y2 = (y3 + ny) / 2
  cp.append((x2, y2))
  p.append((nx, ny))
  x1 = x3
  y1 = y3
  x3 = nx
  y3 = ny



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
s.append('<path d="M%d, %d' % (p[0][0], p[0][1]))
for i in range(1, len(p)):
  s.append("Q%d,%d %d,%d" % (cp[i][0], cp[i][1], p[i][0], p[i][1]))
s.append('" \n fill="none" stroke="red" stroke-width="30"  /> ')
curve =  " ".join(s)

s = []
s.append('<g fill="black">')
s.append('<circle cx="0" cy="0" r="200"/>')
for i in ts:
  s.append('<circle cx="%d" cy="%d" r="100"/>' % (i[0,0], i[0, 1]))
s.append('</g>')
cps = "\n".join(s)

f = open("t1.svg", "w")
print >> f, header
print >> f, curve
print >> f, cps
print >> f, footer
f.close()
