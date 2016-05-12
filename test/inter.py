from lagranto import pyLagranto
from lagranto.caltra import caltra

filename = 'datadir/xjjhq/xjjhqa_036.pp'
spt1, uut1, vvt1, wwt1, p3t1 = caltra.load_winds(filename)

(nx, ny, nz, xmin, ymin,
 dx, dy, hem, per) = caltra.grid_parameters(filename)

print 'nx = ' + str(nx)
print 'ny = ' + str(ny)
print 'nz = ' + str(nz)

i, j, k = pyLagranto.inter.get_index4(
    360, 0, 5000, 0, p3t1, p3t1, spt1, spt1, 1, nx, ny, nz, xmin, ymin,
    dx, dy, -1000)

print i
print j
print k

i, j, k = pyLagranto.inter.get_index3(
    360, 0, 5000, 3, p3t1, spt1, nx, ny, nz, xmin, ymin, dx, dy)

print i
print j
print k
