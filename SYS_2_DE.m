function F=SYS_2_DE(alph,x,t)
% Function used for round_F/round_alpha
global Vdc Vac alpha_0 alpha_1 Jm rho L H mu Le eps C0 c P r N

F=zeros(2*N,1);

Omg=alph;

V=(Vdc+Vac)/H*sqrt(eps/mu);

VDC=rho*L^2/(3*mu*Le*(Vdc+Vac)*C0)*Vdc;
VAC=rho*L^2/(3*mu*Le*(Vdc+Vac)*C0)*Vac;

F(1:N)=x(N+1:2*N);
F(3)=-c*x(3)-Jm*(x(1)-x(1)^(-5))/(Jm-2*x(1)^2-x(1)^(-4)+3) + V^2*x(2)^2*x(1)^3/(x(1)^4+r)^2 -P;
F(4)=-alpha_1*x(4)-alpha_0*x(2)/x(1)^4+VDC+VAC*cos(Omg*t);
