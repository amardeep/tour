from numpy import *
from scipy import *
from scipy import linalg
import matplotlib.pyplot as plt

ts = []

minView = [-23469, -19540, -2]
maxView = [21859, 19774, 1282]


# read tour sites
f = open("hw4.tour")
tourSites = []
for line in f:
  s = [float(x) for x in line.split()]
  tourSites.append(s)
ts = mat(tourSites)

plt.plot(ts[1:-1, 0], ts[1:-1, 1], 'ro', ts[:, 0], ts[:, 1], 'b:')
plt.plot([ts[0, 0]], [ts[0, 1]], 'b>')
plt.plot([ts[-1, 0]], [ts[-1, 1]], 'b<')
plt.axis([minView[0], maxView[0] + 5000,
          minView[1] - 5000, maxView[1] + 5000])
plt.title('De-boor control points')
plt.show()


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
s.append("<path d=\"")
for i in ts:
  x = i[0, 0]
  y = i[0, 1]
  if index == 1:
    s.append("M%d,%d" % (x, y))
  elif index == 2:
    s.append("Q0,0 %d,%d" % (x, y))
  else:
    s.append("T%d,%d" % (x, y))
  index = index + 1
s.append('" \n fill="none" stroke="red" stroke-width="30"  /> ')
curve =  " ".join(s)

s = []
s.append('<g fill="black">')
s.append('<circle cx="0" cy="0" r="200"/>')
for i in ts:
  s.append('<circle cx="%d" cy="%d" r="100"/>' % (i[0,0], i[0, 1]))
s.append('</g>')
cps = "\n".join(s)

f = open("t.svg", "w")
print >> f, header
print >> f, curve
print >> f, cps
print >> f, footer
f.close()
