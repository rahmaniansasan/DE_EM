% Continuation of a periodic solution (Predictor-Corrector scheme)
% Forward sweep
clc
clear Wmax omeg
%global Omega % Vac

%eps=1e-2;                              %% Perturbation parameter
b1=-20;   b2=20; Om_st=0.78;    Om_stp=0.5*Omega_n(end,2)+0.26;    %0.5*Omega_n(end,2)-0.22;
Omega=Om_st; %OmeGa=3.49; %Omeg_0+eps*b1;                   % Frequency interval for sweeping
j=1;
etta_0=[1 0.01 0 0]; %1e-6*ones(1,2*N);
x0=zeros(1,2*N*(2*N+1)); x0(1:2*N)=etta_0; S0=eye(2*N);
for i=1:2*N
    x0(2*N*i+1:2*N*(i+1))=S0(i,:);
end

%%% Computing the starting point solution using Shooting scheme
d_etta=ones(1,2*N);
option=odeset('RelTol',1e-7,'AbsTol',1e-7);
Wmax=zeros(1,50); omeg=zeros(1,50); monodromy_Mat=zeros(2*N,2*N); rondx_rondeta=monodromy_Mat;

T=1* 2*pi/Omega;
params.Omega=Omega;
while norm(d_etta)>1e-4      %% !!!!!!!!!!!!!!! tol
    [t,y]=ode15s(@ F1_DE,[0 T],x0,option,params);
    for i=1:2*N
        monodromy_Mat(i,:)=y(end,2*N*i+1:2*N*(i+1));
    end
    %monodromy_Mat=[y(end,3) y(end,4);y(end,5) y(end,6)];
    d_etta=(monodromy_Mat-eye(2*N))\(etta_0-y(end,1:2*N)).';
    etta_0=etta_0+0.5*d_etta.';               %% 0.1 Coefficient
    x0(1:2*N)=etta_0;                         %% x0=[etta_0 1 0 0 1];
end
[t,y]=ode15s(@ F1_DE,[0 T],x0,option,params);
W_max=0;
for i=1
    W_max=y(:,i); W_max2=y(:,i+3);
end
Wmax(j)=max(W_max)-min(W_max);  % Response Amplitude
Wmax2(j)=max(W_max2)-min(W_max2);
omeg(j)=params.Omega; %OmeGa;
for i=1:2*N
    monodromy_Mat(i,:)=y(end,2*N*i+1:2*N*(i+1));
end
Eig_val=eig(monodromy_Mat);
eigen(:,j)=Eig_val; % max(abs(Eig_val));     %%%%% NEWLY added

if abs(Eig_val)<1.001
    plot(omeg(j),Wmax(j),'b.','MarkerSize',7);
    %plot(omeg(j),eigen(j),'b.','MarkerSize',4);
    hold on
else
    plot(omeg(j),Wmax(j),'r.','MarkerSize',7);
    %plot(omeg(j),eigen(j),'r.','MarkerSize',4);
    hold on  
end
X_init(j,:)=x0(1:2*N);
j=j+1

%%% Continuation scheme

X0=zeros(1,2*N*(2*N+2));
X0(1:2*N*(2*N+1))=x0;                % Starting point (obtained from shooting method, a solution)
ds=1e-2; d_ETTA=ones(1,2*N);
r=0.8;                               % Relaxation parameter

tic
while Omega < Om_stp %OmeGa < 3.525 %Omeg_0+eps*b2
    [t,y]=ode45(@ F2_DE,[0 T],X0,option,params);
    for i=1:2*N
        rondx_rondeta(i,:)=y(end,2*N*i+1:2*N*(i+1));
    end
    rondG_rondeta=rondx_rondeta-eye(2*N);
    rondG_rondalpha=y(end,2*N*(2*N+1)+1:2*N*(2*N+2)).';
    Z=-inv(rondG_rondeta)*rondG_rondalpha;
    alpha_prim=(1+norm(Z)^2)^(-0.5);       % +,- specify the path direction
    eta_prim=alpha_prim*Z;
    eta=X0(1:2*N)+ds*eta_prim.';
    Omega=Omega+ds*alpha_prim;    %OmeGa=OmeGa+ds*alpha_prim;
    T=1* 2*pi/Omega;
    params.Omega=Omega;
    X0(1:2*N)=eta;        %% Predicted values for the initial conditions leading to a periodic solution
    num=1;
    %%% Correction scheme (Newton-Raphson method)
    while norm(d_ETTA) > 1e-4   %% This Tolerance is IMPORTANT
        [t,y]=ode45(@ F1_DE,[0 T],X0(1:2*N*(2*N+1)),option,params);
        
        for i=1:2*N
            rondx_rondeta(i,:)=y(end,2*N*i+1:2*N*(i+1));
        end
        rondG_rondeta=rondx_rondeta-eye(2*N);
        
        d_ETTA=-inv(rondG_rondeta)*(y(end,1:2*N)-X0(1:2*N)).';
        X0(1:2*N)=X0(1:2*N)+r*d_ETTA.';
        num=num+1;
    end
    %%%
    [t,y]=ode45(@ F1_DE,[0 T],X0(1:2*N*(2*N+1)),option,params);
    W_max=0;
    for i=1
        W_max=y(:,i); W_max2=y(:,i+3);
    end
    Wmax(j)=max(W_max)-min(W_max);
    Wmax2(j)=max(W_max2)-min(W_max2);
    omeg(j)=params.Omega; %OmeGa;
    
    for i=1:2*N
        monodromy_Mat(i,:)=y(end,2*N*i+1:2*N*(i+1));
    end
    Eig_val=eig(monodromy_Mat);
    eigen(:,j)=Eig_val; %max(abs(Eig_val));     %%%%% NEWLY added
    
    if abs(Eig_val)<1.01
        plot(omeg(j),Wmax(j),'b.','MarkerSize',7);
        %plot(omeg(j),eigen(j),'b.','MarkerSize',4);
        hold on
    else
        plot(omeg(j),Wmax(j),'r.','MarkerSize',7);
        %plot(omeg(j),eigen(j),'r.','MarkerSize',4);
        hold on  
    end
    X_init(j,:)=X0(1:2*N);
    d_ETTA=ones(1,2*N);
    j=j+1
end
toc



%%% FUNCTIONS USED

function dx=F1_DE(t,x,Par)
% Augmented system with Etta
%global Vdc Vac alpha_0 alpha_1 Jm rho L H mu Le eps C0 c P r N Omega

% Par.c=0.15;  P=0;
% Vac=200; %Vdc;

V=(Par.Vdc+Par.Vac)/Par.H*sqrt(Par.eps/Par.mu);

VDC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vdc;
VAC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vac;

dx=zeros(2*Par.N*(2*Par.N+1),1);

for i=1:2*Par.N
    X(i,:)=x(2*Par.N*i+1:2*Par.N*(i+1));
end

DxF=NumJacob(@ SYS_1_DE,x(1:2*Par.N),t,Par);

dx(1:Par.N)=x(Par.N+1:2*Par.N);
dx(3)=-Par.c*x(3)-Par.Jm*(x(1)-x(1)^(-5))/(Par.Jm-2*x(1)^2-x(1)^(-4)+3) + V^2*x(2)^2*x(1)^3/(x(1)^4+Par.r)^2 - Par.P;
dx(4)=- Par.alpha_1*x(4)-Par.alpha_0*x(2)/x(1)^4+VDC+VAC*cos(Par.Omega*t);

for i=1:2*Par.N
    dx(2*Par.N*i+1:2*Par.N*(i+1))=DxF(i,:)*X;
end

end



function dx=F2_DE(t,x,Par)
% Augmented system with Etta and Alpha
%global Vdc Vac alpha_0 alpha_1 Jm rho L H mu Le eps C0 c P r N Omega

V=(Par.Vdc+Par.Vac)/Par.H*sqrt(Par.eps/Par.mu);

VDC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vdc;
VAC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vac;

dx=zeros(2*Par.N*(2*Par.N+2),1);

for i=1:2*Par.N
    X(i,:)=x(2*Par.N*i+1:2*Par.N*(i+1));
end

DxF=NumJacob(@ SYS_1_DE,x(1:2*Par.N),t,Par);

dx(1:Par.N)=x(Par.N+1:2*Par.N);
dx(3)=-Par.c*x(3)-Par.Jm*(x(1)-x(1)^(-5))/(Par.Jm-2*x(1)^2-x(1)^(-4)+3) + V^2*x(2)^2*x(1)^3/(x(1)^4+Par.r)^2 - Par.P;
dx(4)=- Par.alpha_1*x(4)-Par.alpha_0*x(2)/x(1)^4+VDC+VAC*cos(Par.Omega*t);

for i=1:2*Par.N
    dx(2*Par.N*i+1:2*Par.N*(i+1))=DxF(i,:)*X;
end

dx(2*Par.N*(2*Par.N+1)+1:2*Par.N*(2*Par.N+2))=DxF*x(2*Par.N*(2*Par.N+1)+1:2*Par.N*(2*Par.N+2))+NumJacob(@ SYS_2_DE,Par.Omega,x(1:2*Par.N),t,Par);

end




function F=SYS_1_DE(x,t,Par)
% Function used for round_F/round_x
%global Vdc Vac alpha_0 alpha_1 Jm rho L H mu Le eps C0 c P r N Omega

F=zeros(2*Par.N,1);

V=(Par.Vdc+Par.Vac)/Par.H*sqrt(Par.eps/Par.mu);

VDC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vdc;
VAC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vac;

F(1:Par.N)=x(Par.N+1:2*Par.N);
F(3)=-Par.c*x(3)-Par.Jm*(x(1)-x(1)^(-5))/(Par.Jm-2*x(1)^2-x(1)^(-4)+3) + V^2*x(2)^2*x(1)^3/(x(1)^4+Par.r)^2 - Par.P;
F(4)=- Par.alpha_1*x(4)-Par.alpha_0*x(2)/x(1)^4+VDC+VAC*cos(Par.Omega*t);

end



function F=SYS_2_DE(alph,x,t,Par)
% Function used for round_F/round_alpha
%global Vdc Vac alpha_0 alpha_1 Jm rho L H mu Le eps C0 c P r N

F=zeros(2*Par.N,1);

Omg=alph;

V=(Par.Vdc+Par.Vac)/Par.H*sqrt(Par.eps/Par.mu);

VDC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vdc;
VAC=Par.rho*Par.L^2/(3*Par.mu*Par.Le*(Par.Vdc+Par.Vac)*Par.C0)*Par.Vac;

F(1:Par.N)=x(Par.N+1:2*Par.N);
F(3)=-Par.c*x(3)-Par.Jm*(x(1)-x(1)^(-5))/(Par.Jm-2*x(1)^2-x(1)^(-4)+3) + V^2*x(2)^2*x(1)^3/(x(1)^4+Par.r)^2 - Par.P;
F(4)=- Par.alpha_1*x(4)-Par.alpha_0*x(2)/x(1)^4+VDC+VAC*cos(Omg*t);

end


