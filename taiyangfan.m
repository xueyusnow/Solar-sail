clear all
clc
mu=3.986005e14;
J2=1.082626683*0.001;
re=6378140;
i0=pi/4;
j0=0;
e0=eps;
w0=0;
a0=13.371*10^6;
p0=a0*(1-e0^2);
f0=acos((p0-a0)/e0/a0);
h0=sqrt(p0*mu);
r1=h0^2*cos(f0)/(mu*(1+e0*cos(f0)));
r2=h0^2*sin(f0)/(mu*(1+e0*cos(f0)));
r3=h0^2/mu*0;
v1=-mu/h0*sin(f0);
v2=mu/h0*(e0+cos(f0));
v3=mu/h0/(1+e0*cos(f0))*0;
x0=(cos(j0)*cos(w0)-sin(j0)*sin(w0)*cos(i0))*r1+(-cos(j0)*sin(w0)-sin(j0)*cos(i0)*cos(w0))*r2+sin(j0)*sin(i0)*r3;
y0=(sin(j0)*cos(w0)+cos(j0)*cos(i0)*sin(w0))*r1+(-sin(j0)*sin(w0)+cos(j0)*cos(i0)*cos(w0))*r2-cos(j0)*sin(i0)*r3;
z0=sin(i0)*sin(w0)*r1+sin(i0)*cos(w0)*r2+cos(i0)*r3;
vx0=(cos(j0)*cos(w0)-sin(j0)*sin(w0)*cos(i0))*v1+(-cos(j0)*sin(w0)-sin(j0)*cos(i0)*cos(w0))*v2+sin(j0)*sin(i0)*v3;
vy0=(sin(j0)*cos(w0)+cos(j0)*cos(i0)*sin(w0))*v1+(-sin(j0)*sin(w0)+cos(j0)*cos(i0)*cos(w0))*v2-cos(j0)*sin(i0)*v3;
vz0=sin(i0)*sin(w0)*v1+sin(i0)*cos(w0)*v2+cos(i0)*v3;
v0=sqrt(vx0^2+vy0^2+vz0^2);
T=2*pi*sqrt(a0^3/mu);
% t1=3600*24;
options=odeset('RelTol',1e-9,'AbsTol',1e-9);
m0=[x0,y0,z0,vx0,vy0,vz0];
[t,result1]=ode45(@YunDongFangCheng_sun,[0:1:100*T],m0,options); 
x=result1(:,1);
y=result1(:,2);
z=result1(:,3);
vx=result1(:,4);
vy=result1(:,5);
vz=result1(:,6);

% [X,Y,Z]=sphere(20);
% RE=0.6371e7;
% X=RE*X;
% Y=RE*Y;
% Z=RE*Z;
% mesh(X,Y,Z)%绘制地球
% surf(X,Y,Z); 
% xlabel('x');
% ylabel('y');
% hold on;
%将两个合并

r=[x,y,z];
rdemo=sum(abs(r).^2,2).^(1/2);%计算距离
v=[vx,vy,vz];
vdemo=sum(abs(v).^2,2).^(1/2);%速度大小
vr=(x.*vx+y.*vy+z.*vz)./rdemo;%径向速度
h=cross(r,v);%比角动量
hdemo=sum(abs(h).^2,2).^(1/2);%比角动量的模
hz=h(:,3);
i=acos(hz./hdemo);%求轨道倾角

k=[0,0,1];
[hrow,hcol]=size(h);
K=k(ones(1,hrow),:);%创建一个单位的矩阵

N=cross(K,h);
Ndemo=sum(abs(N).^2,2).^(1/2);%交线N的模
RAAN=acos(N(:,1)./Ndemo);%升交点赤经

p=hdemo.*hdemo/mu;
edemo=ones(hrow,1);
for i=1:hrow
    edemo(i)=1/mu*sqrt((2*mu-rdemo(i)*vdemo(i)^2)*rdemo(i)*vr(i)*vr(i)+(mu-rdemo(i)*vdemo(i)^2)^2);%偏心率
%     ra(i)=p(i)./(1-edemo(i));%半长轴
end
ra=p./(1-edemo);
% plot(t,ra);
% figure;
rp=p./(1+edemo);
% plot(t,rp);
% figure;
a=p./(1-edemo.*edemo);
% plot(t,a);
% maxy=max(a)

  
% plot3(x,y,z,'r')%绘制卫星轨道
% figure;
% plot(t,rdemo);
% plot(t,edemo);
% plot(a0,max(ra));
