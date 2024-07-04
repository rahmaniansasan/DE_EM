clc
clear
%global Vdc alpha_0 alpha_1 Jm eps mu H VDC rho L Le Re r C0 Omega_n N
Omega_n=[];

rho=1.2e3;
mu=6.71e4;
eps=6.198e-11;
L=50e-3;
H=1e-3;
Jm=100;
C0=eps*L^2/H;
t_str=L*sqrt(rho/(3*mu));
lamb_st=1+1e-6;
Omg_m=sqrt((Jm*(lamb_st-lamb_st^(-5))/(Jm-2*lamb_st^2-lamb_st^(-4)+3) )/1e-6);   % Untuned mechanical NF (Vdc=0)
Omg_e_dim=1/2 * Omg_m/t_str;                                                             % Electrical NF is syntonized with the Untuned Mechanical NF
Le=1/(Omg_e_dim^2*C0);  % 16e3;
Re=90e9;  % Ohms
r=0;

alpha_0=rho*H/(3*mu*Le*eps);
alpha_1=Re*L/(3*mu*Le)*sqrt(rho/(3*mu));
VDC=rho*L^2/(3*mu*Le*C0);


params.Vac=400;
params.alpha_0=alpha_0;
params.alpha_1=alpha_1;
params.Jm=Jm;
params.rho=rho;
params.L=L;
params.Le=Le;
params.H=H;
params.mu=mu;
params.eps=eps;
params.C0=C0;
params.c=0.15;
params.P=0;
params.r=r;
params.N=2;
params.alpha_1=alpha_1;
params.VDC=VDC;


% Vdc=0.5;   Vac=0.5;
i=1; N=2; j=1;
V=0; x0=[1 0.7];
option=optimset('TolFun',1e-7,'TolX',1e-7);

for Vdc=0 %1e-10 : 2 : 402
    
    V=Vdc/H*sqrt(eps/mu);
    params.Vdc=Vdc;
    [sol,err]=fsolve(@ DE_ResM,x0,option,params);
    
    if norm(err)>1e-3
        break
    end
    
    lamb_st=sol(1);
    for j=1:2
        f(j)=Jm*(lamb_st-lamb_st^(-5))/(Jm-2*lamb_st^2-lamb_st^(-4)+3);
        lamb_st=lamb_st+1e-6;
    end
    round_M=(f(2)-f(1))/1e-6;
    
    K(1,1)=round_M+5*V^2*sol(2)^2/sol(1)^6; K(1,2)=-2*V^2*sol(2)/sol(1)^5;
    K(2,1)=-4*alpha_0*sol(2)/sol(1)^5;      K(2,2)=alpha_0/sol(1)^4;
    %Mass=eye(2);
    [vv,II]=eig(K);
    Omega_n(i,:)=sqrt(diag(II));
    
    X(i,:)=sol;
    v(i)=V;
    %Omeg_0(i)=sqrt((FGDE(Lam(i)+1e-6)-FGDE(Lam(i)))/1e-6);

    
    x0=sol;
    i=i+1;
end

params.Omega_n=Omega_n(end,:);


%%%% FUNCTIONS USED

function F=DE_ResM(x,Par)
%global Vdc alpha_0 Jm eps mu H VDC r
%P=0;
F=zeros(2,1);
V=Par.Vdc/Par.H*sqrt(Par.eps/Par.mu);

F(1)= Par.Jm*(x(1)-x(1)^(-5))/(Par.Jm-2*x(1)^2-x(1)^(-4)+3) - V^2*x(2)^2*x(1)^3/(x(1)^4+Par.r)^2-Par.P;
F(2)= Par.alpha_0*x(2)/x(1)^4-Par.VDC;
end







