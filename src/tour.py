import math
import gts
import numpy as np
from util import *
import matplotlib.pyplot as plt
import bisect


# config params
ml = 200   # max turning radius
Z_MARGIN_MIN = 100.0
MARGIN_D = 40
MARGIN_D = min(MARGIN_D, 50)

# config flags
USE_MARGIN_D_Z = True
USE_MARGIN_D_XY = True


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


def best_direction(hmap, p1, p2, p3, ):
  # get initial direction based on prev and next points
  if p1 is None:
    assert not p3 is None
    best_tang = get_vector_2d(p2, p3)
  elif p3 is None:
    assert not p1 is None
    best_tang = get_vector_2d(p1, p2)
  else:
    best_tang = get_tangent(p1, p2, p3)

  # try to jiggle to a better position
  j = 1
  maxdiff = None
  maxtang = None
  maxtheta = None
  while True:
    theta = int(j/2) * .05
    if theta > math.pi/4: break
    if j%2 == 0: theta *= -1
    j += 1

    tang = rotate_point_2d(best_tang, theta)
    np1 = Point()
    np2 = Point()
    insert_points(p2, tang, ml, np1, np2)
    h =  get_max_height_line(hmap,
                             interpolate_2d(np1, np2, -ml),
                             interpolate_2d(np2, np1, -ml))
    diff = p2.z - h
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
    print "best_direction: no solution"
  print "best_direction:", maxdiff, maxtheta * 180 / math.pi

  return (maxtang, maxdiff)


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


class TourPoint:
  def __init__(self, hmap, v):
    # HeightLookup instance
    self.hmap = hmap
    # specified tour point
    self.v_orig = gts.Point(v.x, v.y, v.z)
    # actual tour point through which we pass near the specified one
    self.v = gts.Point(v.x, v.y, v.z)
    # best direction
    self.tangent = None
    # two 2d control points along the direction
    self.prev = None
    self.next = None

  def setup_point(self):
    self.compute_best_pos_()

  def setup_direction(self, prev, next):
    self.compute_tangent_(prev, next)

  def compute_best_pos_(self):
    if USE_MARGIN_D_Z:
      self.v.z += MARGIN_D

    if USE_MARGIN_D_XY:
      minx, miny, minh = best_point(hmap, self.v, MARGIN_D)
      self.v.x = minx
      self.v.y = miny

  def compute_tangent_(self, prev, next):
    (self.tangent, maxdiff) = best_direction(self.hmap, prev, self.v, next)
    if maxdiff < 0:
      print "Warning: maxdiff < 0 at", str_point(self.v_orig)
    self.prev = gts.Point()
    self.next = gts.Point()
    insert_points(self.v, self.tangent, ml, self.prev, self.next)
    self.prev.z = self.v.z
    self.next.z = self.v.z



class TourSegment:
  def __init__(self, hmap, tp1, tp2):
    self.hmap = hmap
    self.tp1 = tp1
    self.tp2 = tp2

    self.cps_2d = []
    self.samples_2d = []
    self.dists_2d = []
    self.heights_2d = []

  def compute(self, min_dist):
    self.compute_2d_(min_dist)
    self.compute_3d_(min_dist)

  def compute_2d_(self, min_dist):
    # compute 2d control points
    self.compute_2d_cps_()

    # sample spline
    sample_bspline(self.cps_2d, min_dist, self.samples_2d)

    # compute height and cumulative distance at each sample
    d = 0
    samples = self.samples_2d
    for j in range(len(samples)):
      assert samples[j].z == 0
      if j > 0:
        d += dist(samples[j], samples[j-1])
      self.dists_2d.append(d)
      self.heights_2d.append(
          self.hmap.get_height(gts.Point(samples[j].x, samples[j].y)))

  def compute_2d_cps_(self):
    # compute base control points of the segment
    self.cp1 = self.tp1.prev
    self.cp2 = self.tp1.next
    self.cp3 = interpolate_2d(self.tp1.next, self.tp2.prev, 2*ml)
    self.cp4 = interpolate_2d(self.tp2.prev, self.tp1.next, 2*ml)
    self.cp5 = self.tp2.prev
    self.cp6 = self.tp2.next

    # create a copy with z = 0 for 2d control points
    cps = [self.cp1, self.cp2, self.cp3, self.cp4, self.cp5, self.cp6]
    self.cps_2d = [gts.Point(v.x, v.y) for v in cps]

  def get_slope_index_(self, start, h1, h2):
    # slope between sample at start with h=h1 and last point with h=h2
    points_slope = (h2 - h1)/(self.dists_2d[-1] - self.dists_2d[start])

    # maximum slope from start such that it doesn't hit the terrain
    maxslope = None
    maxindex = None
    maxdh1 = None
    for i in range(start + 1, len(self.samples_2d)):
      ds = self.dists_2d[i] - self.dists_2d[start]
      dh = self.heights_2d[i] - h1
      slope = dh/ds
      dh1 = points_slope * ds
      if maxslope is None or slope > maxslope:
        maxslope = slope
        maxindex = i
        maxdh1 = dh1 - dh

    if maxslope > points_slope or maxdh1 < 40:
      return maxindex
    else:
      print ">>", maxdh1
      return None


  def compute_3d_(self, min_dist):
    h1 = self.tp1.v.z
    h2 = self.tp2.v.z
    next_index = 0
    cps = []  # inserted control points to lift the curve

    while True:
      next_index = self.get_slope_index_(next_index, h1, h2)
      if next_index is None:
        break
      p = self.samples_2d[next_index]   # point in 2d to lift to a deboor cp
      h1 = self.heights_2d[next_index] + ml*1.5  # height of that point
      cps.append(gts.Point(p.x, p.y, h1))
      print "cp %.2f %.2f %.2f" % (p.x, p.y, h1)

    # compute height of cp3 and cp4 interpolating for surrounding points
    if len(cps) > 0:
      cp = cps[0]
    else:
      cp = self.cp5
    ds = dist_2d(self.cp2, cp)
    ds1 = dist_2d(self.cp2, self.cp3)
    dz = cp.z - self.cp2.z
    dz1 = ds1/ds * dz
    self.cp3.z = self.cp2.z + dz1

    if len(cps) > 0:
      cp = cps[-1]
    else:
      cp = self.cp2
    ds = dist_2d(cp, self.cp5)
    ds1 = dist_2d(self.cp4, self.cp5)
    dz = cp.z - self.cp5.z
    dz1 = ds1/ds * dz
    self.cp4.z = self.cp5.z + dz1

    # bspline for segment between point i and i+1
    v = [self.cp1, self.cp2, self.cp3] + cps + [self.cp4, self.cp5, self.cp6]
    self.cps_3d = v
    samples = []
    sample_bspline(v, min_dist, samples)
    self.samples_3d = samples

    # update heights and cum distances
    d = 0
    self.dists_3d = []
    self.heights_3d = []

    for j in range(len(samples)):
      if j > 0:
        d += dist(samples[j], samples[j-1])
      self.dists_3d.append(d)
      self.heights_3d.append(
          hmap.get_height(gts.Point(samples[j].x, samples[j].y)))




class Tour:
  def __init__(self, hmap, min_dist = 50):
    self.hmap = hmap  # HeightLookup instance
    self.tps = []     # Tour points
    self.tss = []     # Tour segments
    self.min_dist = min_dist

  def set_tour_points(self, points):
    self.tps = [TourPoint(self.hmap, v) for v in points]
    self.tss = [TourSegment(self.hmap, self.tps[i], self.tps[i+1])
                for i in range(len(self.tps) - 1)]


  def compute_tour(self):
    for tp in self.tps:
      tp.setup_point()

    # compute direction after best positions for point are computed above
    for i in range(len(self.tps)):
      if i == 0:
        prev = None
      else:
        prev = self.tps[i-1].v
      if i == len(self.tps) - 1:
        next = None
      else:
        next = self.tps[i+1].v
      self.tps[i].setup_direction(prev, next)

    for i in range(len(self.tss)):
      self.tss[i].compute(self.min_dist)

    self.original_tps = [o for o in self.tps]
    self.original_tss = [o for o in self.tss]

  def get_tour_segments():
    pass


  def add_tour_point(self, index, point):
    if index < 1 or index > len(self.tps) - 1:
      # insert only between interior points
      print "Index out of range"
      return False

    print dir(point)
    print "Adding new point", index, point.x, point.y, point.z

    # create new point
    tp = TourPoint(self.hmap, point)
    tp.setup_point()

    prev = self.tps[index-1].v
    next = self.tps[index].v
    tp.setup_direction(prev, next)

    # create new segments
    tss1 = TourSegment(self.hmap, self.tps[index-1], tp)
    tss2 = TourSegment(self.hmap, tp, self.tps[index])

    # insert new point and segments
    self.tps.insert(index, tp)
    self.tss[index-1:index] = [tss1, tss2]

    tss1.compute(self.min_dist)
    tss2.compute(self.min_dist)

    return True


  def delete_tour_point(self, index):
    if index < 1 or index > len(self.tps) - 1:
      # insert only between interior points
      print "Index out of range"
      return False

    # compute new segment
    tss = TourSegment(self.hmap, self.tps[index-1], self.tps[index+1])
    tss.compute(self.min_dist)

    # remove original points and segments and insert new one
    self.tps[index:index+1] = []
    self.tss[index-1:index+1] = [tss]

    return True

  def reset_original_tour(self):
    self.tps = [o for o in self.original_tps]
    self.tss = [o for o in self.original_tss]


def plot_tour(tour, pat = 'r-', samples_2d = True):
  x = []
  y = []
  for ts in tour.tss:
    if samples_2d:
      samples = ts.samples_2d
    else:
      samples = ts.samples_3d
    x.extend([o.x for o in samples])
    y.extend([-o.y for o in samples])
  plt.plot(x, y, pat)
  plt.axis([-20000, 20000, -20000, 20000])


def plot_height_map(tour):
  cum_dist = 0
  for ts in tour.tss:
    plt.plot([cum_dist], [ts.tp1.v.z], 'rx')
    x = [cum_dist + d for d in ts.dists_3d]
    cum_dist += ts.dists_3d[-1]
    plt.plot(x, ts.heights_3d, 'b-')
    y = [o.z for o in ts.samples_3d]
    plt.plot(x, y, 'g-')
  plt.plot([cum_dist], [tour.tss[-1].tp2.v.z], 'rx')


def write_bez_cp(tour, fname="hw4.bez"):
  fout = open(fname, "w")
  for i in range(len(tour.tss)):
    ts = tour.tss[i]
    for j in range(len(ts.cps_3d) - 2):
      p0 = mid_point(ts.cps_3d[j], ts.cps_3d[j+1])
      p1 = ts.cps_3d[j+1]
      p2 = mid_point(ts.cps_3d[j+1], ts.cps_3d[j+2])
      print >> fout, i, p0.x, p0.y, p0.z
      print >> fout, i, p1.x, p1.y, p1.z
      print >> fout, i, p2.x, p2.y, p2.z
  fout.close()


def write_bspline_cp(tour, fname="hw4.spline"):
  fout = open(fname, "w")
  for i in range(len(tour.tss)):
    ts = tour.tss[i]
    for j in range(len(ts.cps_3d) - 2):
      p1 = ts.cps_3d[j]
      print >> fout, p1.x, p1.y, p1.z
    p1 = tour.tss[-1].cps_3d[-2]
    print >> fout, p1.x, p1.y, p1.z
    p1 = tour.tss[-1].cps_3d[-1]
    print >> fout, p1.x, p1.y, p1.z
  fout.close()

def write_stats(tour, fname="hw4.stats"):
  tour_len = 0
  num_cps = 0
  min_height = None
  for ts in tour.tss:
    tour_len += ts.dists_3d[-1]
    num_cps += len(ts.cps_3d) - 2
    for s in ts.samples_3d:
      h1 = s.z
      h2 = hmap.get_height(gts.Point(s.x, s.y))
      if min_height is None or h1-h2 < min_height:
        min_height = h1-h2
  num_cps += 2

  f = open(fname, "w")
  print >> f, "tour-length:", tour_len
  print >> f, "min-height:", min_height
  print >> f, "num-control-points", num_cps
  print "tour-length:", tour_len
  print "min-height:", min_height
  print "num-control-points", num_cps


def write_samples(tour, fname="hw4.sample"):
  f = open(fname, "w")
  for ts in tour.tss:
    for s in ts.samples_3d:
      print >> f, s.x, s.y, s.z
  f.close()


# global stuff

vertices = []
triangle_indices = []

# read terrain info
read_vertices(vertices, "hw4.vertices")
read_triangles(triangle_indices, "hw4.faces")
(vmin, vmax, vspan) = get_extents(vertices)
print "Terrain extent: ", vmin, vmax
print "Terrain span: ", vspan


# create HeightLookup map
hmap = HeightLookup(vertices, triangle_indices)

# read tour points
tour_points = []
read_tour(tour_points, "hw4.tour")
#read_tour(tour_points, "tour.tsp")

# create and initialize tour instance
tour = Tour(hmap, 10)
tour.set_tour_points(tour_points)
tour.compute_tour()

write_bez_cp(tour)
write_bspline_cp(tour)
write_stats(tour)
write_samples(tour)


if __name__ == "__main__":

  plt.figure(1)
  plot_tour(tour, 'r-', True)
  plot_tour(tour, 'b-', False)

  plt.figure(2)
  plot_height_map(tour)
  plt.show()
