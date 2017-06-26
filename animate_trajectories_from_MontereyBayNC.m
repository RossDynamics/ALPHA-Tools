clear all
close all
clc
ncfile='MontereyBay.nc';
ncid=netcdf.open(ncfile,'NC_NOWRITE');
%Grab NC data
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
%Loop through NC data infomation
for i = [1:1:nvars]
    i-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,i-1);
    varname
end
time = netcdf.getVar(ncid,0);
size(time)
lon = netcdf.getVar(ncid,3);
size(lon)
lat = netcdf.getVar(ncid,4);
size(lat)

maxlat = 38;
maxlon = -121;
minlat = 36;
minlon = -124;

figure
for i =1:length(time)
    clf
    scatter(lon(i,:),lat(i,:),'k.')
    axis([minlon maxlon minlat maxlat])
    %drawnow
    pause(0.05)
end
