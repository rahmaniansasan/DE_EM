function dx=F2_DE(t,x)
% Augmented system with Etta and Alpha
global Vdc Vac alpha_0 alpha_1 Jm rho L H mu Le eps C0 c P r N Omega

V=(Vdc+Vac)/H*sqrt(eps/mu);

VDC=rho*L^2/(3*mu*Le*(Vdc+Vac)*C0)*Vdc;
VAC=rho*L^2/(3*mu*Le*(Vdc+Vac)*C0)*Vac;

dx=zeros(2*N*(2*N+2),1);

for i=1:2*N
    X(i,:)=x(2*N*i+1:2*N*(i+1));
end

DxF=NumJacob(@ SYS_1_DE,x(1:2*N),t);

dx(1:N)=x(N+1:2*N);
dx(3)=-c*x(3)-Jm*(x(1)-x(1)^(-5))/(Jm-2*x(1)^2-x(1)^(-4)+3) + V^2*x(2)^2*x(1)^3/(x(1)^4+r)^2 -P;
dx(4)=-alpha_1*x(4)-alpha_0*x(2)/x(1)^4+VDC+VAC*cos(Omega*t);

for i=1:2*N
    dx(2*N*i+1:2*N*(i+1))=DxF(i,:)*X;
end

dx(2*N*(2*N+1)+1:2*N*(2*N+2))=DxF*x(2*N*(2*N+1)+1:2*N*(2*N+2))+NumJacob(@ SYS_2_DE,Omega,x(1:2*N),t);
