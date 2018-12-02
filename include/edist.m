function [ adist, azi, bazi ] = edist(slat,slon,xlat,xlon)

%  [ adist, az, baz ] = edist(slat,slon,xlat,xlon)
% returns angular distance, azimuth, and back azimuth
%
%  Angular distance calculated using vector dot product.
%
%  Azimuth and back azimuth calculated using code converted
%  from the fortran routine edist which was adapted by Rene 
%  van Hissenhoven, 1984, using Andoyer Lambert distance formula.
%  Converted to matlab function by Dave Wilson, Jun, 2002
%
%      Reference: Thomas, Paul D., geodesic arc length on the 
%      reference ellopsiod to sencod order terms in the flattening;
%      J.G.R.,70,14,1965. Equation 1B, page 3331.  

lat1 = (90-slat)*pi/180;
lon1 = slon*pi/180;
x1 = sin(lat1)*cos(lon1);
y1 = sin(lat1)*sin(lon1);
z1 = cos(lat1);

elat=(90-xlat)*pi/180;
elon=xlon*pi/180;
x2 = sin(elat)*cos(elon);
y2 = sin(elat)*sin(elon);
z2 = cos(elat);

adist = real(acos(x1*x2 + y1*y2 + z1*z2)*180/pi);
clear lat1 lon1 elat elon

% the code below is converted from the fortran program edist
PI = 3.141592653589793238462643;
       f=1.0/298.26;
        a = 6378.135;
       rcd=0.0174532925199432957692 ;

%      f,a,rcd are respectively flattening,major radious, deg to rad.
%      for the international ellipsoid.

       flat=f;
       r=57.29577951308232;
       dtor=1.0/r;
       lat1=slat*dtor;
       lon1=slon*dtor;
       lat2=xlat*dtor;
       lon2=xlon*dtor;

       if lon1 < 0.0,
               lon1=-lon1;
       else
               lon1=2.0*PI-lon1;
       end

       if lon2<0.0,
               lon2=-lon2;
       else
               lon2=2.0*PI-lon2;
       end

       ncon=0;
%      interchange points
 
       if lon1>lon2,
               savlon=lon1;
               savlat=lat1;
               lon1=lon2;
               lat1=lat2;
               lon2=savlon;
               lat2=savlat;
               ncon=1;
       end

       if (lon2-lon1) > PI,
               savlon=lon1;
               savlat=lat1;
               lon1=lon2;
               lat1=lat2;
               lon2=savlon+2.0*PI;
               lat2=savlat;
               if ncon < 1,
                       ncon=1;
               else
                       ncon=0;
               end
       end

%      above checks for east+, west- longs.
%      and north+, south-, latitudes.

       phim=0.5*(lat1+lat2);      
       dphim=0.5*(lat2-lat1);
       dl=(lon2-lon1);
       dlm=0.5*dl;
       cdphi=cos(dphim);
       sdphi=sin(dphim);
       cdlm=cos(dlm);

       sk=sin(phim)*cos(dphim);
       bk=sin(dphim)*cos(phim);
       bh=cos(dphim)*cos(dphim)-sin(phim)*sin(phim);
       bl=sin(dphim)*sin(dphim)+bh*sin(dlm)*sin(dlm);

       cdel=1.0-2.0*bl;
       del=acos(cdel);

       t=del/sin(del);
       u=2.0*sk*sk/(1.0-bl);
       v=2.0*bk*bk/bl;
       be=60.0*cos(del);
       x=u+v;
       y=u-v;

       bd=8.0*(6.0+t*t);
       ba=4.0*t+(16.0+be*t/15.0);
       bb=2.0*bd;
       bc=2.0*t-0.5*(ba+be);
       delta=-(flat/4.0)*(t*x-3.0*y);

%      distance in kilometers

       range=a*sin(del)*(t+delta+(flat*flat/128.0)*(ba*x+bb*y+bc*x*x+bd*x*y+be*y*y));

%      azimuth given as east of north.
%      note that Thomas gives azi as west of south.

       sindl=sin(dl);
       da2pd=-0.5*flat*bh*(t+1.0)*(bk*sindl/bl);
       da2md=-0.5*flat*bh*(t-1.0)*(bk*sindl/(1.0-bl));
       cd2=cos(del/2.0);
       arg1=sdphi*cdlm/sin(del/2.0);
       arg2=cos(PI/2.0-phim)*sin(dlm)/cd2;
       if arg1 > 1.0,
               %disp(arg1)
               arg1=2.0-arg1;
       elseif arg1 < -1.0,
               %disp(arg1)
               arg1=-2.0-arg1;
       end
       a2pa1=asin(arg1)+PI;
       if (arg2>1.0),
               %disp(arg2)
               arg2=2.0-arg2;
       elseif arg2 < -1.0,
               %disp(arg2)
               arg2=-2.0-arg2;
       end  
       a2ma1=PI-acos(arg2);
       a12=a2pa1+a2ma1+da2pd-da2md;
       a21=a2pa1-a2ma1+da2pd+da2md;

       if (ncon>=1)
               azi=a21;
               bazi=a12;
       else
               azi=a12;
               bazi=a21;
       end
               azi=azi*180/PI;
               bazi=bazi*180/PI;

% this is the angular distance which the origonal edist gives
% adist=(range/a)*180/pi;
% the vector dot product is better

return;

