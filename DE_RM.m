function F=DE_RM(t,x,Par)
%global Vdc alpha_0 alpha_1 Jm eps mu H rho L Le C0 r Omega_n
F=zeros(4,1);

Omeg=1.2; %Par.Omega_n(2);

V=(Par.Vdc+Par.Vac)/Par.H*sqrt(Par.eps/Par.mu);

VDC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vdc;
VAC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vac;

F(1)=x(3);
F(2)=x(4);
F(3)=-Par.c*x(3)-Par.Jm*(x(1)-x(1)^(-5))/(Par.Jm-2*x(1)^2-x(1)^(-4)+3) + V^2*x(2)^2*x(1)^3/(x(1)^4+Par.r)^2 - Par.P;
F(4)=- Par.alpha_1*x(4)-Par.alpha_0*x(2)/x(1)^4+VDC+VAC*cos(Omeg*t);
