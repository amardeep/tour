import math
import gts
import numpy as np
from util import *
import matplotlib.pyplot as plt
import bisect

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
    # specified tour point
    self.v_orig = v
    # actual tour point through which we pass near the specified one
    self.v = gts.Point(v.x, v.y, v.z)

ml = 200   # max turning radius
Z_MARGIN_MIN = 100.0
MARGIN_D = 40
MARGIN_D = min(MARGIN_D, 50)

use_margind_z = True
use_margind_xy = True


# compute point at which height is the lowest
# by comparing height at the boundary of the 2dx2d square
def best_point(hmap, v, d):
  incr = max(d/100.0, 1)
  minx = miny = minh = None
  for x in np.arange(v.x-d, v.x+d, incr):
    h1 = hmap.get_height(gts.Point(x, v.y-d))
    h2 = hmap.get_height(gts.Point(x, v.y+d))
    if minh is None or h1 < minh:
      minx, miny, minh = x, v.y-d, h1
    if minh is None or h2 < minh:
      minx, miny, minh = x, v.y+d, h2
  for y in np.arange(v.y-d, v.y+d, incr):
    h1 = hmap.get_height(gts.Point(v.x-d, y))
    h2 = hmap.get_height(gts.Point(v.x+d, y))
    if minh is None or h1 < minh:
      minx, miny, minh = v.x-d, y, h1
    if minh is None or h2 < minh:
      minx, miny, minh = v.x+d, y, h2
  return (minx, miny, minh)


TP = []
for v in T:
  tp = TourPoint(v)
  # add MARGIN_D to increase height
  if use_margind_z:
    tp.v.z += MARGIN_D
  if use_margind_xy:
    minx, miny, minh = best_point(hmap, v, MARGIN_D)
    tp.v.x = minx
    tp.v.y = miny
  TP.append(tp)

# compute tangent direction and control points for interior points
for i in range(1, len(TP) - 1):
  tp = TP[i]
  p1 = TP[i-1].v
  p2 = TP[i].v
  p3 = TP[i+1].v
  tp.tangent = get_tangent(p1, p2, p3)

  j = 1
  maxdiff = None
  maxtang = None
  maxtheta = None
  while True:
    theta = int(j/2) * .05
    if theta > math.pi/4: break
    if j%2 == 0: theta *= -1
    j += 1

    tang = rotate_point_2d(tp.tangent, theta)
    np1 = Point()
    np2 = Point()
    insert_points(p2, tang, ml, np1, np2)
    h =  get_max_height_line(hmap,
                             interpolate_2d(np1, np2, -ml),
                             interpolate_2d(np2, np1, -ml))
    diff = tp.v.z - h
    #print tp.v.z, diff
    if diff > Z_MARGIN_MIN:
      maxtang = tang
      maxdiff = diff
      maxtheta = theta
      break

    if maxdiff is None or diff > maxdiff:
      maxdiff = diff
      maxtang = tang
      maxtheta = theta

  if (maxdiff < 0):
    print "no solution"
    #assert(False)
  print maxdiff, maxtheta * 180 / math.pi

  tp.tangent = maxtang
  np1 = Point()
  np2 = Point()
  insert_points(p2, maxtang, ml, np1, np2)
  tp.prev = np1
  tp.next = np2


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

def sample_bspline(cp, min_dist, samples, distfn = dist):
  for i in range(len(cp) - 2):
    p0 = mid_point(cp[i], cp[i+1])
    p1 = Point(cp[i+1].x, cp[i+1].y, cp[i+1].z)
    p2 = mid_point(cp[i+1], cp[i+2])
    samples.append(p0)
    p0.t = 0
    p2.t = 1
    samples.extend(sample_bezier(p0, p1, p2, p0, p2, min_dist, distfn))
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
    v = [tp.cp1, tp.cp2, tp.cp3, tp.cp4, tp.cp5, tp1.cp1, tp1.cp2]
    # make a copy with z = 0
    v = [gts.Point(a.x, a.y) for a in v]
    sample_bspline(v, min_dist, samples)

    d = 0
    for j in range(len(samples)):
      assert samples[j].z == 0
      if j > 0:
        d += dist(samples[j], samples[j-1])
      dists.append(d)
      heights.append(
          hmap.get_height(gts.Point(samples[j].x, samples[j].y)))

sample_tour_2d(TP, hmap, 50)


def get_slope_index(tp, start, h1, h2):
  points_slope = (h2 - h1)/(tp.dists_2d[-1] - tp.dists_2d[start])
  maxslope = None
  maxindex = None
  for i in range(start + 1, len(tp.samples_2d)):
    ds = tp.dists_2d[i] - tp.dists_2d[start]
    dh = tp.heights_2d[i] - h1
    slope = dh/ds
    if maxslope is None or slope > maxslope:
      maxslope = slope
      maxindex = i

  if maxslope > points_slope:
    return maxindex
  else:
    return None

def sample_tour_3d(TP, hmap, min_dist):
  for i in range(len(TP)-1):
    tp = TP[i]
    tp1 = TP[i+1]
    tp.cp1.z = tp.cp2.z = tp.v.z
    tp1.cp1.z = tp1.cp2.z = tp1.v.z

    h1 = tp.v.z
    h2 = tp1.v.z
    next_index = 0
    tp.cp3ds = []
    while True:
      next_index = get_slope_index(tp, next_index, h1, h2)
      if next_index is None:
        break
      p = tp.samples_2d[next_index]
      h1 = tp.heights_2d[next_index] + ml
      tp.cp3ds.append(gts.Point(p.x, p.y, h1))
      print "cp %.2f %.2f %.2f" % (p.x, p.y, h1)

    # compute heights of control points cp3 and cp5 interpolating from nearby
    if len(tp.cp3ds) > 0:
      cp = tp.cp3ds[0]
    else:
      cp = tp1.cp1
    ds = dist_2d(tp.cp2, cp)
    dz = cp.z - tp.v.z
    ds1 = dist_2d(tp.cp2, tp.cp3)
    dz1 = ds1/ds * dz
    tp.cp3.z = tp.v.z + dz1

    if len(tp.cp3ds) > 0:
      cp = tp.cp3ds[-1]
    else:
      cp = tp.cp2
    ds = dist_2d(cp, tp1.cp1)
    dz = cp.z - tp1.v.z
    ds1 = dist_2d(tp.cp5, tp1.cp1)
    dz1 = ds1/ds * dz
    tp.cp5.z = tp1.v.z + dz1

    # bspline for segment between point i and i+1
    v = [tp.cp1, tp.cp2, tp.cp3] + tp.cp3ds + [tp.cp5, tp1.cp1, tp1.cp2]
    s = []
    sample_bspline(v, min_dist, s)
    tp.samples_3d = s

    # update heights and cum distances
    d = 0
    samples = tp.samples_3d
    dists = tp.dists_3d = []
    heights = tp.heights_3d = []

    for j in range(len(samples)):
      if j > 0:
        d += dist(samples[j], samples[j-1])
      dists.append(d)
      heights.append(
          hmap.get_height(gts.Point(samples[j].x, samples[j].y)))


sample_tour_3d(TP, hmap, 50)



def plot_bspline(cp, min_dist, pat = 'b-'):
  samples = []
  sample_bspline(cp, min_dist, samples)
  x = [a.x for a in samples]
  y = [-a.y for a in samples]
  plt.plot(x, y, pat)
  plt.axis([-20000, 20000, -20000, 20000])

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
  return (x, y)


plt.figure(1)
plot_bspline(v, 100, 'b-')

#plt.plot(TP[0].dists_2d, TP[0].heights_2d, 'r-')

plt.figure(2)
cum_dist = 0
for i in range(len(TP) - 1):
  tp = TP[i]
  plt.plot([cum_dist], [tp.v.z], 'rx')
  x = [cum_dist + d for d in tp.dists_3d]
  cum_dist += tp.dists_2d[-1]
  plt.plot(x, tp.heights_3d, 'r-')
  y = [a.z for a in tp.samples_3d]
  plt.plot(x, y, 'g-')
plt.plot([cum_dist], [TP[-1].v.z], 'rx')

#(x, y) = plot_heights(v, 200, 'b-')

#v1 = [TP[0].cp1, TP[0].cp2, TP[0].cp3, TP[0].cp5, TP[1].cp1, TP[1].cp2]
#plot_bspline(v1, 100, 'b-')
#plot_heights(v1, 100, 'b-')


plt.show()


fpoints = open("points.txt", "w")
for s in np.arange(0, 1, .01):
  for t in np.arange(0, 1, .01):
    x = vmin.x + vspan.x * s
    y = vmin.y + vspan.y * t
    z = hmap.get_height(gts.Point(x, y))
    if not z is None:
      print >> fpoints, x, y, z
fpoints.close()

fsamples = open("samples.txt", "w")
for i in range(len(TP) - 1):
  tp = TP[i]
  for s in tp.samples_3d:
    print >> fsamples, s.x, s.y, s.z
fsamples.close()
