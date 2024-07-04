% Frequency response curve near primary resonance (using N modes)
% Long time simulation method
clc
clear lam_1 lam_2 omeg
global OMEG

% eps=1e-2;                   %% Perturbation parameter
% b1=-900; b2=5; hs=3;      % !!!!!   IMPORTANT   Freq-resolution (h) =0.01
% sgma=b1: hs : b2;
Omeg_int=[1: 1e-3 : 1.46]; %Omeg_0+eps*sgma;   % Frequency interval for sweeping
j=1;
x0=[1 0];  X_init=x0; % Initial guess

option=odeset('RelTol',1e-9,'AbsTol',1e-9);
lam_1=zeros(1,50); lam_2=zeros(1,50); omeg=zeros(1,50); w1=0; w2=0;
for OMEG=Omeg_int
    [t,y]=ode45(@ DE_RM_2D,[0 300],x0,option);
    w1=y(:,1);
    w2=y(:,2);
    
    lam_1(j)=max(abs(w1(floor(0.8*length(t)):end))) - min(abs(w1(floor(0.8*length(t)):end)));
    lam_2(j)=max(w2(floor(0.8*length(t)):end)) - min(w2(floor(0.8*length(t)):end));
    omeg(j)=OMEG;
    x0=y(end,:);
    plot(OMEG,lam_1(j),'b.','MarkerSize',5); hold on
    w1=0; w2=0;
    j=j+1
end