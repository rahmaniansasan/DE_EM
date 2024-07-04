% Pseudo-Arclength Continuation scheme 
clc
%global Omega % Vac
jj=input('insert the number of star point: ');
alpha_st=omeg(jj); %% alpha_str
eta_st=X_init(jj,:);
rondx_rondeta=zeros(2*N,2*N);   monodromy_Mat=rondx_rondeta;
X0(1:2*N)=eta_st;
Omega=alpha_st; %OmeGa=alpha_st;
T=1* 2*pi/Omega;
params.Omega=Omega;

dd=zeros(1:length(eta_st)); dd=0.008*ones(1,2*N); %% !! Tuning parameter (VERY IMPORTANT)

[t,y]=ode45(@ F2_DE,[0 T],X0,option,params);
for i=1:2*N
    rondx_rondeta(i,:)=y(end,2*N*i+1:2*N*(i+1));
end
rondG_rondeta=rondx_rondeta-eye(2*N);
rondG_rondalpha=y(end,2*N*(2*N+1)+1:2*N*(2*N+2)).';
Z=-inv(rondG_rondeta)*rondG_rondalpha;
alpha_str_prim=(1+norm(Z)^2)^(-0.5);       % +,- specify the path direction
eta_str_prim=alpha_prim*Z;

% ds=0.005;
% eta_1=eta_st+ds*eta_str_prim.';
% alpha_1=Omega+ds*alpha_str_prim;   Not needed here
%%%
% OmeGa=alpha_1;
X0(1:2*N)=eta_st+dd;
[t,y]=ode45(@ F1_DE,[0 T],X0(1:2*N*(2*N+1)),option,params);
W_max=0;
for i=1
    W_max=y(:,i); W_max2=y(:,i+3);
end
WMAX=max(W_max)-min(W_max);
plot(Omega,WMAX,'bo','MarkerSize',5)

T=1* 2*pi/Omega;   d_etta=ones(1,2*N);  x0=X0(1:2*N*(2*N+1));  etta_0=x0(1:2*N);
params.Omega=Omega;
while norm(d_etta)>1e-4      %% !!!!!!!!!!!!!!! tol
    [t,y]=ode15s(@ F1_DE,[0 T],x0,option,params);
    for i=1:2*N
        monodromy_Mat(i,:)=y(end,2*N*i+1:2*N*(i+1));
    end
    %monodromy_Mat=[y(end,3) y(end,4);y(end,5) y(end,6)];
    d_etta=(monodromy_Mat-eye(2*N))\(etta_0-y(end,1:2*N)).';
    etta_0=etta_0+0.4*d_etta.';               %% 0.1 Coefficient
    x0(1:2*N)=etta_0;                         %% x0=[etta_0 1 0 0 1];
end
[t,y]=ode15s(@ F1_DE,[0 T],x0,option,params);
W_max=0;
for i=1
    W_max=y(:,i); W_max2=y(:,i+3);
end
Wmax(j)=max(W_max)-min(W_max);
Wmax2(j)=max(W_max2)-min(W_max2);
omeg(j)=params.Omega;
for i=1:2*N
    monodromy_Mat(i,:)=y(end,2*N*i+1:2*N*(i+1));
end
Eig_val=eig(monodromy_Mat);
if abs(Eig_val)<1 %1.001
    plot(omeg(j),Wmax(j),'g*','MarkerSize',4);
    hold on
else
    plot(omeg(j),Wmax(j),'r.','MarkerSize',4);
    hold on  
end
X_init(j,:)=X0(1:2*N);
j=j+1;





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






