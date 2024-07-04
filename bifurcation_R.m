clc
%global R Omega Beta5

% Omega=2;
% T=2*pi/Omega;

tend=2500;
X0=[1 0.01 0 0];
opt=odeset('RelTol',1e-7,'AbsTol',1e-7);
i=1;
X1=zeros(100,100);
X2=zeros(100,100);
tic
for Omega= 1.1782 : 0.0005 : 1.2284
    params.Omega=Omega;
    T=2*pi/Omega;
    point=tend:-T:tend-99*T;
    sol=ode45(@ DE_RM,[0 tend],X0,opt,params);
    lam=deval(sol,point,1);
    curnt=deval(sol,point,4);

    X1(:,i)=max(lam)-min(lam);
    X2(:,i)=max(curnt)-min(curnt);
    r(i)=Omega;
    X0=sol.y(:,end).';
    i=i+1
end
toc

%%% Bifurcation plot procedure

for i=1:length(r)  %% Suppose that r is your continuation parameter
    for j=1:size(X1,1)  
        plot(r(i),X1(j,i),'r.','MarkerSize',5);
        hold on
    end
end


%%% Poincare Map
% Please Note, that here, X1 and X2 are assumed to be in Matrix form containing data for different values of clearance size
figure
for i=1:size(X1,1)
    plot(X1(i,1),X2(i,2),'k.','MarkerSize',5);
    hold on
end




% FUNCTION USED

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

end


