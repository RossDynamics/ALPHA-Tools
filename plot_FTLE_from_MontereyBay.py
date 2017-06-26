# coding: utf-8
import numpy
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.tri as tri
nc = netCDF4.Dataset('MontereyBay.nc')
init_lon = nc.variables['initial_lon']
init_lat = nc.variables['initial_lat']
f = nc.variables['FTLE']
x = init_lon[::30]
y = init_lat[:30]
xx,yy = numpy.meshgrid(x,y)
f0 = f[:,-1,0]
f00 = numpy.reshape(f0,(30,30),'F')
a=plt.contourf(xx,yy,f00,300)
plt.colorbar(a)
plt.show()
