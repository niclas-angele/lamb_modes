clear; close all; 

%%%parameters and frequency w 
w=3; 
h=1; 
ct=0.5;
cl=0.9;
mu=ct^2; 
lambda=cl^2-2*ct^2;

%%%lamb modes 
p=@(k) sqrt(w.^2./(cl^2)-k.^2);
q=@(k) sqrt(w.^2./(ct^2)-k.^2);
jns=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*sin(q(k)*h).^2+h*4*k.^2.*p(k).^2.*sin(p(k)*h).^2+sin(p(k)*h).*sin(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*sin(q(k)*h).*cos(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*sin(p(k)*h).*cos(q(k)*h)));
us=@(k,y) (1i*k*(q(k)^2-k^2)*sin(q(k)*h)*cos(p(k)*y)-2*1i*k*p(k)*q(k)*sin(p(k)*h)*cos(q(k)*y)); 
vs=@(k,y) (-p(k)*(q(k)^2-k^2)*sin(q(k)*h)*sin(p(k)*y)-2*k^2*p(k)*sin(p(k)*h)*sin(q(k)*y));
ss=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*sin(q(k)*h)*cos(p(k)*y)+4*mu*p(k)*q(k)*k^2*sin(p(k)*h)*cos(q(k)*y);
ts=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(-sin(q(k)*h)*sin(p(k)*y)+sin(p(k)*h)*sin(q(k)*y));
ua=@(k,y) 1i*k*(q(k)^2-k^2)*cos(q(k)*h)*sin(p(k)*y)-2*1i*k*p(k)*q(k)*cos(p(k)*h)*sin(q(k)*y);
va=@(k,y) p(k)*(q(k)^2-k^2)*cos(q(k)*h)*cos(p(k)*y)+2*k^2*p(k)*cos(p(k)*h)*cos(q(k)*y);
sa=@(k,y) -(q(k)^2-k^2)*(cl^2*k^2+lambda*p(k)^2)*cos(q(k)*h)*sin(p(k)*y)+4*mu*p(k)*q(k)*k^2*cos(p(k)*h)*sin(q(k)*y);
ta=@(k,y) 2*1i*k*mu*(q(k)^2-k^2)*p(k)*(cos(q(k)*h)*cos(p(k)*y)-cos(p(k)*h)*cos(q(k)*y));
jna=@(k) 1i*mu*k.*(q(k).^2+k.^2).*(h*(q(k).^2-k.^2).^2.*cos(q(k)*h).^2+h*4*k.^2.*p(k).^2.*cos(p(k)*h).^2-cos(p(k)*h).*cos(q(k)*h).*(1./p(k).*(q(k).^2-k.^2).*(q(k).^2-k.^2-8*p(k).^2).*cos(q(k)*h).*sin(p(k)*h)+4./q(k).*p(k).^2.*(2*q(k).^2-k.^2).*cos(p(k)*h).*sin(q(k)*h)));

%%%solution of the dispersion relation 
[S,A]=dispersionh(w,5*w,h);

%%%symmetric mode
k=S(2); 
x=linspace(-2.2,2.2,400); 
y=linspace(-1,1,200); 
hold on 
%position of the moving sections
pos=[-2,-1.5,-1,-0.5,0,0.5,1,1.5,2]; 
for i=1:length(pos)
    plot(pos(i)+0.0008*exp(1i*k*pos(i))*us(k,y), y+0.0008*exp(1i*k*pos(i))*vs(k,y),'k')
end
plot(x+0.0008*exp(1i*k*x)*us(S(2),-1),-1+0.0008*exp(1i*k*x)*vs(k,-1),'k')
plot(x+0.0008*exp(1i*k*x)*us(S(2),1),1+0.0008*exp(1i*k*x)*vs(k,1),'k')


%%%antisymmetric mode
figure
k=A(2); 
x=linspace(-1.1,1.1,400); 
y=linspace(-1,1,200); 
hold on 
%position of the moving sections
pos=[-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1]; 
for i=1:length(pos)
    plot(pos(i)+0.00004*exp(1i*k*pos(i))*ua(k,y), y+0.00004*exp(1i*k*pos(i))*va(k,y),'k')
end
plot(x+0.00004*exp(1i*k*x)*ua(k,-1),-1+0.00004*exp(1i*k*x)*va(k,-1),'k')
plot(x+0.00004*exp(1i*k*x)*ua(k,1),1+0.00004*exp(1i*k*x)*va(k,1),'k')
axis([-1.1,1.1,-1.3,1.3])
