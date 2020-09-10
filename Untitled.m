clear all
clc
mu=3.986005e14;
J2=1.082626683*0.001;
re=6378140;
i0=pi/6;
j0=0;
e0=0;
w0=0;
f0=0;
a0=6.971*10^6;
p0=a0*(1-e0^2);
h0=sqrt(p0*mu);
r1=h0^2*cos(f0)/(mu*(1+e0*cos(f0)));
r2=h0^2*sin(f0)/(mu*(1+e0*cos(f0)));
r3=h0^2*sin(f0)*0;
v1=-mu/h0*sin(f0);
v2=mu/h0*(e0+cos(f0));
v3=mu/h0*0;
x0=(cos(j0)*cos(w0)-sin(j0)*sin(w0)*cos(i0))*r1+(-cos(j0)*sin(w0)-sin(j0)*cos(i0)*cos(w0))*r2+sin(j0)*sin(i0)*r3;
y0=(sin(j0)*cos(w0)+cos(j0)*cos(i0)*sin(w0))*r1+(-sin(j0)*sin(w0)+cos(j0)*cos(i0)*cos(w0))*r2-cos(j0)*sin(i0)*r3;
z0=sin(i0)*sin(w0)*r1+sin(i0)*cos(w0)*r2+cos(i0)*r3;
vx0=(cos(j0)*cos(w0)-sin(j0)*sin(w0)*cos(i0))*v1+(-cos(j0)*sin(w0)-sin(j0)*cos(i0)*cos(w0))*v2+sin(j0)*sin(i0)*v3;
vy0=(sin(j0)*cos(w0)+cos(j0)*cos(i0)*sin(w0))*v1+(-sin(j0)*sin(w0)+cos(j0)*cos(i0)*cos(w0))*v2-cos(j0)*sin(i0)*v3;
vz0=sin(i0)*sin(w0)*v1+sin(i0)*cos(w0)*v2+cos(i0)*v3;
v0=sqrt(vx0^2+vy0^2+vz0^2);
% [X,Y,Z]=sphere(30);
% RE=6.371*10^6;
% X=RE*X;
% Y=RE*Y;
% Z=RE*Z;
% grid on;
% mesh(X,Y,Z)
% surf(X,Y,Z)
% hold on;
T=2*pi*sqrt(a0^3/mu);
% t1=3600*24;
options=odeset('RelTol',1e-9,'AbsTol',1e-9);
m0=[x0,y0,z0,vx0,vy0,vz0];
% m0=[6.971e+06,0,0,0,6.5487e+03,3.7809e+03];
[t,f]=ode45('ready',0:0.01:20*T,m0,options);
x=f(:,1);
y=f(:,2);
z=f(:,3);
vx=f(:,4);
vy=f(:,5);
vz=f(:,6);
r=[x,y,z];
r1=sqrt(x.^2+y.^2+z.^2);
v=[vx,vy,vz];
v1=sqrt(vx.*vx+vy.*vy+vz.*vz);
vr=(x.*vx+y.*vy+z.*vz)./r1;
h2=[y.*vz-z.*vy;
    z.*vx-x.*vz;
    x.*vy-y.*vx];
hx=y.*vz-z.*vy;
hy=z.*vx-x.*vz;
hz=x.*vy-y.*vx;
h1=sqrt(hx.^2+hy.^2+hz.^2);
% i=acos(hz./h1);
% nx=hy;
% ny=-hx;
% nz=0;
% n1=sqrt(hx.^2+hy.^2);
% if(hx<=0)
% j=acos(hy./n1);
% else
%     j=2*pi-acos(hy./n1);
% end
ex=1/mu*((v1.*v1-mu./r1).*x-r1.*vr.*vx);
ey=1/mu*((v1.*v1-mu./r1).*y-r1.*vr.*vy);
ez=1/mu*((v1.*v1-mu./r1).*z-r1.*vr.*vz);
e=sqrt(ex.*ex+ey.*ey+ez.*ez);
p=h1.^2/mu;
ra=p./(1+e);
plot(t,ra);
% % % e=sqrt((2*mu-r1.*v1.*v1).*r1.*vr.*vr+(mu-r1.*v1.*v1).^2)/mu; 
% if(ez>=0)
% w=acos((nx.*ex+ny.*ey+nz.*ez)./n1./e)*180/pi; 
% else
%     w=360-acos((nx.*ex+ny.*ey+nz.*ez)./n1./e)*180/pi;
% end
% if(vr>=0)
% f=acos((h1.^2/mu./r1-1)./e);
% else
%     f=2*pi-acos((h1.^2./mu./r1-1)./e);
% end
% f1=sqrt(f.*f);
% plot(t,j);
% figure;
% plot(t,r1);
% xlabel('时间t');ylabel('轨道半径');
% figure;
% plot(t,j);
% xlabel('时间t');ylabel('升交点赤经j');
% figure;
% plot(t,w);
% xlabel('时间t');ylabel('近地点幅角w');
% figure;
% plot(t,e1);
% xlabel('时间t');ylabel('偏心率e');
% figure;
% plot(t,f1);
% xlabel('时间t');ylabel('真近点角f');
% plot3(x,y,z,'r')
