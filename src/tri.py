import math
import gts
import numpy as np
from util import *
import matplotlib.pyplot as plt

V = []  # terrain map vertices but in form of gts.Vertex array
F_indices = []  # triangle vertex indices
F = []  # faces
T = []  # tour

read_vertices(V)
read_triangles(F_indices)
make_gts_triangles(V, F_indices, F)
read_tour(T)

(vmin, vmax, vspan) = get_extents(V)
print "Terrain extent: ", vmin, vmax
print "Terrain span: ", vspan

hmap = HeightLookup(V, F)

class TourPoint:
  def __init__(self, v):
    self.v = v   # tour point

TP = []
for v in T:
  TP.append(TourPoint(v))

ml = 1000   # max turning radius


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


# compute tangent direction and control points for interior points
for i in range(1, len(TP) - 1):
  tp = TP[i]
  p1 = TP[i-1].v
  p2 = TP[i].v
  p3 = TP[i+1].v
  tp.tangent = get_tangent(p1, p2, p3)

  # for theta in np.arange(-math.pi/4, math.pi/4, .05):
  #   tang = rotate_point_2d(tp.tangent, theta)

  # if i == 1: break

  np1 = Point()
  np2 = Point()
  insert_points(p2, tp.tangent, ml, np1, np2)
  tp.prev = np1
  tp.next = np2
  print get_max_height_line(hmap,
                            interpolate_2d(np1, np2, -ml),
                            interpolate_2d(np2, np1, -ml)),
  print tp.v.z

# compute tangent direction for end points of tour
TP[0].tangent = get_vector_2d(TP[0].v, TP[1].prev)
TP[0].prev = interpolate_2d(TP[0].v, TP[1].prev, -ml)
TP[0].next = interpolate_2d(TP[0].v, TP[1].prev, ml)

TP[-1].tangent = get_vector_2d(TP[-1].v, TP[-2].next)
TP[-1].prev = interpolate_2d(TP[-1].v, TP[-2].next, ml)
TP[-1].next = interpolate_2d(TP[-1].v, TP[-2].next, -ml)

# compute control points
for i in range(len(TP)):
  tp = TP[i]
  tp.cp1 = tp.prev
  tp.cp2 = tp.next
  if i == len(TP)-1:
    tp.cp3 = tp.cp4 = tp.cp5 = None
  else:
    tp1 = TP[i+1]
    tp.cp3 = interpolate_2d(tp.next, tp1.prev, 2*ml)
    tp.cp5 = interpolate_2d(tp1.prev, tp.next, 2*ml)
    tp.cp4 = Point((tp.cp3.x + tp.cp5.x)/2, (tp.cp3.y + tp.cp5.y)/2)

v = []
for tp in TP:
  v.append(tp.cp1)
  v.append(tp.cp2)
  if not tp.cp3 is None:
    v.append(tp.cp3)
    #v.append(tp.cp4)
    v.append(tp.cp5)

# x = [a.x for a in v]
# y = [-a.y for a in v]
# plt.plot(x, y, 'r-')
# plt.axis([-20000, 20000, -20000, 20000])
# plt.show()

def sample_bspline(cp, min_dist, samples):
  for i in range(len(cp) - 2):
    p0 = mid_point(cp[i], cp[i+1])
    p1 = Point(cp[i+1].x, cp[i+1].y, cp[i+1].z)
    p2 = mid_point(cp[i+1], cp[i+2])
    samples.append(p0)
    p0.t = 0
    p2.t = 1
    samples.extend(sample_bezier(p0, p1, p2, p0, p2, min_dist))
  samples.append(mid_point(cp[-2], cp[-1]))

def sample_tour_2d(TP, hmap, min_dist):
  # go through all the tour points
  for i in range (len(TP)-1):
    tp = TP[i]
    tp1 = TP[i+1]
    # bspline samples between point i and i+1
    samples = tp.samples_2d = []
    # cumulative arc length at each point
    dists = tp.dists_2d = []
    # heights at each sample point
    heights = tp.heights_2d = []

    # bspline for segment between point i and i+1
    v = [tp.cp1, tp.cp2, tp.cp3, tp.cp5, tp1.cp1, tp1.cp2]
    sample_bspline(v, min_dist, samples)

    d = 0
    for j in range(len(samples)):
      assert samples[j].z == 0
      if j > 0:
        d += dist(samples[j], samples[j-1])
      dists.append(d)
      heights.append(
          hmap.get_height(gts.Point(samples[j].x, samples[j].y)))

sample_tour_2d(TP, hmap, 200)

def plot_bspline(cp, min_dist, pat = 'b-'):
  samples = []
  sample_bspline(cp, min_dist, samples)
  x = [a.x for a in samples]
  y = [-a.y for a in samples]
  plt.plot(x, y, pat)
  plt.axis([-20000, 20000, -20000, 20000])
  plt.show()

def plot_heights(cp, min_dist, pat = 'b-'):
  samples = []
  sample_bspline(cp, min_dist, samples)
  x = []
  y = []
  d = 0
  for i in range(len(samples)):
    assert samples[i].z == 0
    if i > 0:
      d += dist(samples[i], samples[i-1])
    x.append(d)
    y.append(hmap.get_height(gts.Point(samples[i].x, samples[i].y)))
  plt.plot(x, y, pat)
  plt.axis([0, 200000, -10, 1300])
  plt.show()
  return (x, y)


#plot_bspline(v, 500, 'b,')

#plt.plot(TP[0].dists_2d, TP[0].heights_2d, 'r-')

cum_dist = 0
for i in range(len(TP) - 1):
  tp = TP[i]
  plt.plot([cum_dist], [tp.v.z], 'rx')
  x = [cum_dist + d for d in tp.dists_2d]
  cum_dist += tp.dists_2d[-1]
  plt.plot(x, tp.heights_2d, 'r-')
plt.plot([cum_dist], [TP[-1].v.z], 'rx')

(x, y) = plot_heights(v, 200, 'b-')

#v1 = [TP[0].cp1, TP[0].cp2, TP[0].cp3, TP[0].cp5, TP[1].cp1, TP[1].cp2]
#plot_bspline(v1, 100, 'b-')
#plot_heights(v1, 100, 'b-')




fpoints = open("points.txt", "w")
for s in np.arange(0, 1, .01):
  for t in np.arange(0, 1, .01):
    x = vmin.x + vspan.x * s
    y = vmin.y + vspan.y * t
    z = hmap.get_height(gts.Point(x, y))
    if not z is None:
      print >> fpoints, x, y, z
fpoints.close()
