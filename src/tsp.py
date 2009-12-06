from util import *
import random
import matplotlib.pyplot as plt

T = []
read_tour(T, "hw4.tour")

tour = range(len(T))

dist_cache = {}

def mydist(i, j):
  if i < j:
    a, b = i, j
  else:
    a, b = j, i
  return dist_cache.setdefault((a,b), dist_2d(T[a], T[b]))

def tour_length(tour):
  d = 0.0
  for i in range(len(tour) - 1):
    d += mydist(tour[i], tour[i+1])
  return d


def hill_climb(tour):
  t = tour[:]   # initial tour
  d = tour_length(t)

  while True:
    change = False
    for i in range(len(t)-1):
      # try swapping index i
      t1 = t[:]
      t1[i], t1[i+1] = t1[i+1], t1[i]
      d1 = tour_length(t1)
      if d1 < d:
        d = d1
        t = t1
        change = True
    if change == False:
      return t, d

def hill_climb_with_random_restart(tour):
  t, d = hill_climb(tour)
  print d

  for i in range(2000):
    t1 = t[:]
    a, b = random.sample(xrange(len(tour)), 2)
    t1[a], t1[b] = t1[b], t1[a]
    a, b = random.sample(xrange(len(tour)), 2)
    t1[a], t1[b] = t1[b], t1[a]
    a, b = random.sample(xrange(len(tour)), 2)
    t1[a], t1[b] = t1[b], t1[a]
    t1, d1 = hill_climb(t1)
    if d1 < d:
      t, d = t1, d1
      print i, d
  return t,d


t,d = hill_climb_with_random_restart(tour)

x = [T[i].x for i in t]
y = [T[i].y for i in t]
plt.plot(x, y, 'b-')
plt.axis([-20000, 20000, -20000, 20000])
plt.show()

ftsp = open("tour.tsp", "w")
for i in t:
  print >> ftsp, T[i].x, T[i].y, T[i].z
ftsp.close()
