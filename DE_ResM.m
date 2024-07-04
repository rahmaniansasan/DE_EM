function F=DE_ResM(x)
global Vdc alpha_0 Jm eps mu H VDC r
P=0;
F=zeros(2,1);
V=Vdc/H*sqrt(eps/mu);

F(1)= Jm*(x(1)-x(1)^(-5))/(Jm-2*x(1)^2-x(1)^(-4)+3) - V^2*x(2)^2*x(1)^3/(x(1)^4+r)^2-P;
F(2)= alpha_0*x(2)/x(1)^4-VDC;
