# coding: utf-8
import numpy
import netCDF4
import matplotlib.pyplot as plt
import matplotlib.tri as tri
nc = netCDF4.Dataset('montbay.nc')
init_lon = nc.variables['initial_lon']
init_lat = nc.variables['initial_lat']
f = nc.variables['FTLE']
print init_lon.shape
print init_lat.shape
#x = init_lon[::200]
x = init_lon[::30]
print x[:]
y = init_lat[:30]
print y[:]
xx,yy = numpy.meshgrid(x,y)
print f.shape
f0 = f[:,-1,0]
f00 = numpy.reshape(f0,(30,30),'F')
a=plt.contourf(xx,yy,f00,300)
plt.colorbar(a)
plt.show()
a=plt.contourf(xx,yy,f00,300)
plt.colorbar(a)
plt.show()
#get_ipython().magic('save current_session ~0/')
