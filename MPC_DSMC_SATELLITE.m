clear
clc

R2D = 180/pi;
ToRPM = 60/2/pi;

%% ========================================================================
%% Satellite Dynamic Parameters:
J=[20 1.5 1;1.5 18 1.8;1 1.8 12];
Jw=diag([10;6;4])*10^-3; % Wheel momentums
WheelMax = 2000; % [rpm] Wheel max speed
WheelMax = WheelMax/ToRPM; % [rad/s] Wheel max speed
Umax = 0.8; % N.m. maximum torque
Umin = -Umax;

mu = 398600; % km^3/s^2
Rc = 6371; % km
rc = 700+Rc; % km
wo = sqrt(mu/(rc^3));

%% ========================================================================
%% Simulation/General Parameters:

Omega_0 = 0*[0.01;-0.03;0.02];
Angle_0 = 0*[10;0;-5]*pi/180;

Ts = 0.1;   % Sample Time
Tsim = 60; % (sec) Simulation Time
Nsim = floor(Tsim/Ts);

Flag_PreKnown = 0;  % (=1): Pre-Known
                    % (=0): Pre-UnKnown
%% ========================================================================
%% Uncertainty Parameters:
Flag_Uncert = 1;  % (1: On) (0: off)

d_f = Flag_Uncert*0.5;  % 50% uncertainty
d_g = (1+Flag_Uncert*0.3); 
% Note: d_g =1 means no uncertainty
%       d_g =1.5 means 50% uncertainty
I_hat = d_g*J;
I_hat_inv = inv(I_hat);

%% ========================================================================
%% Noise & Disturbance Parameters:
Flag_Dist = 0; % (=0): No Disturbance.
                     % (=1): Apply Disturbance.
Flag_Noise = 0; % (=0): No Noise
                % (=1): Apply Noise
Flag_Filter = 0; % (=0): No Filter
                 % (=1): Apply Filter

Tau_dist = [-0.4; 0.4; -0.4];
NoisePower = Flag_Noise*7;  % Sensor noise power in percent (7%)

%% ========================================================================
%% ------- Apply Saturation on RW (reaction wheels) -----------------
Flag_Sat = 0; % (=1): Saturated
                    % (=0): NotSaturated

%% ========================================================================
%% -------- Filter parameters:
%    H(z^-1) = (1-a)/2*(1+z^-1)/(1-az^-1)
%    H(z^-1) = abar^-1*(1+z^-1)/(1-az^-1)
a = 0.6;
abar = 2/(1-a);
% -------------------------------------------------------------------------

%% ========================================================================
%% ===== MPC Controller parameters ========================================
n=3; % number of states
m=3; % number of inputs (RW)
A = eye(n);
B = eye(n,m);
C = eye(n);

Np = 20;  % Prediction Horizon for MPC
% We assume: Nc = Np
Q = 0.01*eye(n);
R = 100*eye(n);

P = dare(A,B,Q,R);  % NOTE: Function "dare(.)" solves Riccati Eq. (49) in Paper.
% P is equal to "X" in this function.
% See: "Help dare" in MATLAB.
P = (P+P')/2; % Just to ensure symmetry
Qbar = blkdiag(ndiag(Q,Np-1),P); % Eq. (50)
Rbar = ndiag(R,Np); % Eq. (50)

%% Constraints of RWs: -----------
Wmax = 0.2*[1;1;1];
Wmin = -Wmax;
Vmax  = Ts*Wmax/sqrt(6);  % See Lemma 3 in Paper
Vmin = -Vmax;

dVmax = 10000*[1;1;1];
dVmin = -dVmax;

%% ----------------------------------------------------------------------


%% ========================================================================
% -------- Prediction Matrices  --------------------
% -------------------------------------------------------------------------
% See: Eq. (51) in Paper:
%   y_fut = Hdv*dV_fut + F_fut
% where:
%   F_fut = H1*LI*v(t-1) + LI*d(t) + H2*Theta(t)
% -------------------------------------------------------------------------
CI = [eye(m*Np); -eye(m*Np)];
EI0 = [eye(m), zeros(m,m*(Np-1))];
EI = EI0;
for i=2:Np
    EI0 = [eye(m), EI0(:,1:end-m)];
    EI = [EI;EI0];
end
LI = ncat(eye(m),Np);
Wdv = [ncat(dVmax,Np);ncat(-dVmin,Np)];
Wv  = [ncat(Vmax,Np);ncat(-Vmin,Np)];

WW = [Wdv; Wv];
Gdv = [CI; CI*EI];
Gv = [ncat(zeros(m),2*Np); CI*LI];

H10 = [C*B, zeros(m,m*(Np-1))];
H1 = H10;
H2 = C*A;
for i=1:Np-1
    H10 = [C*A^i*B, H10(:,1:end-m)];
    H1 = [H1;H10];
    H2 = [H2;C*A^(i+1)];
end
Hdv = H1*EI;
% -------------------------------------------------------------------------

%% ========================================================================
%% -------- Prediction with Filter:
%   H(z^-1) = (1-a)/2*(1+z^-1)/(1-az^-1)
%   H(z^-1) = abar^-1*(1+z^-1)/(1-az^-1)
% -------------------------------------------------------------------------
% See: Eqs. (52)-(57) in Paper:
%   y_fut = Hdv_bar*dV_fut + Ff_fut_bar
% where:
%   Ff_fut_bar = Ca*Ff_fut -a*abar*La*yf(t)-La*y(t)+...
%                Hdv_bar*La*dv(t-1) + a*abar*Hdv_bar*La*dvf(t-1)
% -------------------------------------------------------------------------
C1 = [eye(3)      zeros(n,3*(Np-1))];
C2 = [abar*eye(3)      zeros(n,3*(Np-1))];
D = [eye(3); zeros(3*(Np-1),3)];
if Np>=2
    C10 = [eye(3)  eye(3)    zeros(3,3*(Np-2))];
    C20 = [-a*abar*eye(3)  abar*eye(3)    zeros(3,3*(Np-2))];
    C1 = [C1; C10];
    C2 = [C2; C20];
end
for i=3:Np
    C10 = [zeros(3,3) C10(:,1:end-3)];
    C20 = [zeros(3,3) C20(:,1:end-3)];
    C1 = [C1; C10];
    C2 = [C2; C20];
end
Ca = inv(C1)*C2;
La = inv(C1)*D;

Hdv_bar = Ca*Hdv*inv(Ca);
S1 = (Hdv_bar'*Qbar*Hdv_bar+Rbar);  % NOTE: {H=2*S1} IN QP: min {0.5 X'HX + f'X}
S1 = 0.5*(S1+S1'); % Just to ensure symmetry
%% -------------------------------------------------------------------------

%% =========================================================================
% ===== DSMC Controller (Inner Loop) =======
K1 = 1000*eye(3);  
K2 = 10;
K = 5;


%% ========================================================================


%% ************************************************************************
%        Reference Generation:
% ************************************************************************
Ref0 = [10;-10;-5]*pi/180;
Ref1 = [0;-10;-5]*pi/180;
Ref2 = [0;10;-5]*pi/180;
Ref3 = [0;10;-15]*pi/180;

T0 = 15;  % (sec)
T1 = 3;  % (sec)
T2 = 2;  % (sec)
T3 = Tsim-(T0+T1+T2);  % (sec)

N0 = floor(T0/Ts);
N1 = floor(T1/Ts);
N2 = floor(T2/Ts);
N3 = floor(T3/Ts);

Yd = [ncat(Ref0,N0);ncat(Ref1,N1);ncat(Ref2,N2);ncat(Ref3,N3)];
Yd = [Yd;ncat(Ref3,Np)];
%% ************************************************************************


%% ************************************************************************
%        START Simulation:
% ************************************************************************

%% Initial Values at t=0:
Omega(:,1)=Omega_0;      % Satellite's angular rate
Angle = Angle_0;         % Satellite's Euler angles
Omega_w(:,1)=zeros(3,1); % Wheels' angular rate
W = zeros(3,Nsim);
dt = zeros(3,Nsim);
V = zeros(3,Nsim);
Wd = zeros(3,Nsim);
Omega_d = zeros(3,Nsim);
S_Omega = zeros(3,Nsim);
u = zeros(3,Nsim);

Vt_1 = zeros(3,1);   % V(t-1)
Fft_1 = zeros(3,1);  % Ff(t-1)
Ft_1 = zeros(3,1);   % F(t-1)
yt_1 = zeros(3,1);   % y(t-1)
yft_1 = zeros(3,1);  % yf(t-1)
dVft_1 = zeros(3,1); % dVf(t-1)
dVt_1 = zeros(3,1);  % dV(t-1)
Omega_d_1 = zeros(3,1); % Omega_d(t-1)
Omega_e(:,1) = Omega_d_1-Omega(:,1);
options = optimset('Display','off');
for i=1:Nsim
    Phi = Angle(1,i);
    Teta = Angle(2,i);
    Sai = Angle(3,i);
    
    % See: first paragraph after Eq. (4)
    jo = [cos(Teta)*sin(Sai);
        cos(Phi)*cos(Sai)+sin(Phi)*sin(Teta)*sin(Sai);
        -sin(Phi)*cos(Sai)+cos(Phi)*sin(Teta)*sin(Sai)];
    W(:,i) = Omega(:,i)- wo*jo;
    
    % See: Eq. (5)
    R = [1  tan(Teta)*sin(Phi)  tan(Teta)*cos(Phi);
        0      cos(Phi)        -sin(Phi);
        0  sin(Phi)/cos(Teta)  cos(Phi)/cos(Teta)];
    
    % Rinv = R^-1

    Rinv = [1       0       -sin(Teta);
        0   cos(Phi)  sin(Phi)*cos(Teta);
        0   -sin(Phi) cos(Phi)*cos(Teta)];
    
    y_model = [Phi; Teta; Sai];
    
    %----------- Add Measurment Noise -------------------------------------
    Noise = y_model.*(rand(3,1)-0.5)*2*NoisePower/100;
    y_plant = y_model+Noise;
    %----------------------------------------------------------------------
    
    dt(:,i)= y_plant - y_model; % Disturbance
    
    %----------------------------------------------------------------------
    if Flag_PreKnown==1
        Theta_d_fut = Yd((i-1)*m+1:(i-1+Np)*m);
    else
        Theta_d_fut = ncat(Yd((i-1)*m+1:(i-1)*m+3),Np);
    end
    %----------------------------------------------------------------------
    
    yt = y_plant;  % measured output => yt = Theta(t)
    if Flag_Filter==1
        % Apply Filter H(z) to y(t):
        %   H(z^-1) = abar^-1*(1+z^-1)/(1-az^-1)
        %   yft = H(z^-1)*yt >>>
        %   --> (1-az^-1)*yft = abar^-1*(1+z^-1)*yt >>>
        %   --> yft = abar^-1*(yt+yt_1)+a*yft_1;
        yft = abar^-1*(yt+yt_1)+a*yft_1; % Filtered Output
    else
        yft = yt; % No-Filtered Output
    end
    %----------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % -------------- Solve MPC: (See Eq. (58) in paper) -------------------
    %     min { dV_fut'*S1*dV_fut + S2*dV_fut }
    %      s.t. Gdv*dV_fut <= W - Gv*v(t-1)
    % ---------------------------------------------------------------------
    
    F_fut = H1*LI*Vt_1 + LI*dt(:,i) + H2*yft; % See Eq. (51)
    Ff_fut = inv(Ca)*F_fut + a*abar*inv(Ca)*La*Fft_1 + inv(Ca)*La*Ft_1; % See Eq. (56)
    
    Ff_fut_bar = Ca*Ff_fut - a*abar*La*yft - La*yt + Hdv_bar*La*dVt_1 + Hdv_bar*a*abar*La*dVft_1; % See Eq. (57)
    
    S2 = 2*(Ff_fut_bar-Theta_d_fut)'*Qbar*Hdv_bar;
    
    % ------------ Quadratic Programming -------------------
    %     min { 0.5*x'*H*x + f'*x }
    %      s.t. AA*x <= bb
    % ---------------------------------------------------------------------
    H = 2*S1;
    f = S2';
    AA = Gdv;
    bb = WW - Gv*Vt_1;
    dV_fut = quadprog(H,f,AA,bb,[],[],[],[],[],options);
    
    dVt = dV_fut(1:m);
    dVft = abar^-1*(dVt+dVt_1)+a*dVft_1; % Filtered input
    
    V(:,i) = dVt + Vt_1;
    Vt_1 = V(:,i);
    Wd(:,i) = 1/Ts*Rinv*V(:,i);
    Omega_d(:,i) = Wd(:,i) +  wo*jo; % it must be >>>  Wd(:,i) + wo*jo;
    
    %---- Updating "y(t-1) & yf(t-1)" for next sample time
    yt_1 = yt;
    yft_1 = yft;
    %-------------------------------------
    %---- Updating "F(t-1) & Ff(t-1)" for next sample time
    Ft_1 = F_fut(1:3);
    Fft_1 = Ff_fut(1:3);
    %-------------------------------------
    %---- Updating "dVf_1_past" for next sample time
    dVt_1 = dVt;
    dVft_1 = dVft;
    %-------------------------------------
    
    
    % ------------------------------------------------------------------
    % -------------- DSMC Design ---------------------------------------
    %      u_dot = inv(K2)*(K1*(Omega_d_dot - fhat - I_hat*u) + K*Sign(S_Omega))
    %      u = Ts*u_dot + u
    % ------------------------------------------------------------------
    
    % Adding Disturbance to the system
    t0_d = 35;   % Disturbance start time
    t1_d = 37.5; % Disturbance stop time
    t2_d = 40;   % Disturbance stop time
    t=i*Ts;
    if Flag_Dist && t>=t0_d && t<=t1_d
        Tau_d = Tau_dist; % External disturbance Torque
    elseif Flag_Dist && t>=t1_d && t<=t2_d
        Tau_d = -Tau_dist; % External disturbance Torque
    else
        Tau_d = zeros(3,1); % External disturbance Torque
    end
    
    Tau_gg = zeros(3,1);  % Tau_gg = 3wo^2*cross(ko,J*ko): it is neglected, very small.
    
    f_Omega = -inv(J)*cross(Omega(:,i),J*Omega(:,i))+inv(J)*Tau_d; % See: (9)
    if i==1
        Omega_d_1 = Omega_d(:,i);
    end
    Omega_d_dot = (Omega_d(:,i) - Omega_d_1)/Ts;
    Omega_d_1 = Omega_d(:,i);
    fhat = -inv(I_hat)*cross(Omega(:,i),I_hat*Omega(:,i))+0*inv(I_hat)*Tau_d;% inv(I_hat)*Taud >> Measurable disturbance
    Omega_e(:,i) = Omega_d(:,i)-Omega(:,i);

    S_Omega(: , i) = K1*Omega_e(: , i) - K2*u(: , i);
    u_dot(: , i) = inv(K2)*(K1*(Omega_d_dot - fhat - I_hat_inv*u(: , i)) + K*(S_Omega(: , i)));
    u(: , i+1) = Ts*u_dot(: , i) + u(: , i);
    
    % ------- Convert Control input 'u' to the Reaction Wheels speed
    %   u = -Jw*Omega_w_dot - cross(Omega, Jw*Omega_w) + Tau_gg
    
    % Apply External Saturation:
    if Flag_Sat == 1
        u(:,i)=Usaturate(u(:,i),Umax,Umin);
    end
    
    Omega_w_dot= inv(Jw)*( -u(:,i)- cross(Omega(:,i), Jw*Omega_w(:,i))+Tau_gg );
    Omega_w(:,i+1) = Ts*Omega_w_dot + Omega_w(:,i);
    
    % ------------- Apply RW Saturation ----------------------------------
    if Flag_Sat == 1
        [u(:,i), Omega_w(:,i+1)] = Saturate(Omega_w(:,i+1),Omega_w(:,i),Omega(:,i),WheelMax,Ts,Jw,Tau_gg);
    end
    
    % ---------------------------------------------------------------------
    
    %% ------------- Apply Control Law-------------------------------------
    Omega(:,i+1) = Ts*(f_Omega + inv(J)*u(:,i)) + Omega(:,i);
    Angle(:,i+1) = Ts*R*W(:,i)+y_model;
    % ---------------------------------------------------------------------
    
    Ym(:,i)=y_model;
    Yp(:,i)=y_plant;
    Yf(:,i)=yft;
end
%% ************************************************************************


%% ************************************************************************
%% PLOT Results:
%% ************************************************************************
Fig=[1 2 3];

figure(Fig(1))
clf
subplot(3,1,1);  hold on;
Time = 0:Ts:Tsim-Ts;
dN =1;% round(Nsim/500);
plot(Time(1:dN:Nsim),R2D*Angle(1,1:dN:Nsim),'--k','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*Angle(2,1:dN:Nsim),':b','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*Angle(3,1:dN:Nsim),'-r','LineWidth',1.5);

plot(Time,R2D*Yd(1:m:m*Nsim)','r','LineWidth',1);
plot(Time,R2D*Yd(2:m:m*Nsim)','m','LineWidth',1);
plot(Time,R2D*Yd(3:m:m*Nsim)','b','LineWidth',1);
ylabel({'$(\phi, \theta, \psi)$ (deg)'},'Interpreter','latex','FontSize',12);
legend('\phi','\theta','\psi')

subplot(3,1,2); hold on;
plot(Time(1:dN:Nsim),R2D*W(1,1:dN:Nsim),'--k','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*W(2,1:dN:Nsim),':b','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*W(3,1:dN:Nsim),'-r','LineWidth',1.5);
ylabel('$\omega(t)$ (deg/s)','Interpreter','latex','FontSize',12); legend('\omega_x','\omega_y','\omega_z')


subplot(3,1,3); hold on;
plot(Time(1:dN:Nsim),R2D*Omega_e(1,1:dN:Nsim),'--k','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*Omega_e(2,1:dN:Nsim),':b','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*Omega_e(3,1:dN:Nsim),'-r','LineWidth',1.5);
xlabel('Time (sec)'); ylabel('$\Omega_e(t)$ (deg/s)','Interpreter','latex','FontSize',12);
legend('\Omega_{ex}','\Omega_{ey}','\Omega_{ez}')

figure(Fig(2))
clf
subplot(3,1,1); hold on;
plot(Time(1:dN:Nsim),R2D*V(1,1:dN:Nsim),'--k','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*V(2,1:dN:Nsim),':b','LineWidth',1.5);
plot(Time(1:dN:Nsim),R2D*V(3,1:dN:Nsim),'-r','LineWidth',1.5);
ylabel('$v(t)$ (deg/s)','Interpreter','latex','FontSize',12);
legend('v_{x}','v_{y}','v_{z}')

subplot(3,1,2);  hold on;
plot(Time(1:dN:Nsim),ToRPM*Omega_w(1,1:dN:Nsim),'--k','LineWidth',1.5);
plot(Time(1:dN:Nsim),ToRPM*Omega_w(2,1:dN:Nsim),':b','LineWidth',1.5);
plot(Time(1:dN:Nsim),ToRPM*Omega_w(3,1:dN:Nsim),'-r','LineWidth',1.5);
ylabel('$\Omega_w(t)$ (rpm)','Interpreter','latex','FontSize',12);
legend('\Omega_{wx}','\Omega_{wy}','\Omega_{wz}')

subplot(3,1,3);  hold on;
plot(Time(1:dN:Nsim),u(1,1:dN:Nsim),'--k','LineWidth',1.5);
plot(Time(1:dN:Nsim),u(2,1:dN:Nsim),':b','LineWidth',1.5);
plot(Time(1:dN:Nsim),u(3,1:dN:Nsim),'-r','LineWidth',1.5);
xlabel('Time (sec)'); ylabel('$\Omega_e(t)$ (deg/s)','Interpreter','latex','FontSize',12);
ylabel('$u(t)$ (Nm)','Interpreter','latex','FontSize',12);
legend('u_{x}','u_{y}','u_{z}')

%% ----------------------------------
figure(Fig(3))
clf
hold on;
plot(Time(1:dN:Nsim),S_Omega(1,1:dN:Nsim),'--k','LineWidth',1.5);
plot(Time(1:dN:Nsim),S_Omega(2,1:dN:Nsim),':b','LineWidth',1.5);
plot(Time(1:dN:Nsim),S_Omega(3,1:dN:Nsim),'-r','LineWidth',1.5);
ylabel('S_{\Omega}(t)','FontSize',12);
legend('S_{{\Omega}x}','S_{{\Omega}y}','S_{{\Omega}z}')

xlabel('Time (sec)','FontSize',12);



%% ========================================================================
%% Extra Required Functions:
%% ========================================================================
function Pn = ndiag(P,n)
Pn = P;
for i=2:n
    Pn = blkdiag(Pn,P);
end
end
%% ========================================================================
function Pn = ncat(P,n)
Pn = P;
for i=2:n
    Pn = [Pn;P];
end
end
%% ========================================================================
function [Usat,  OmegaWsat]= Saturate(OmegaW,OmegaW_1,Omega,WheelMax,Ts,Jw,Tau_gg)
for i=1:3
    if OmegaW(i) > WheelMax
        OmegaWsat(i,1) = WheelMax;
    elseif OmegaW(i) < -WheelMax
        OmegaWsat(i,1) = -WheelMax;
    else
        OmegaWsat(i,1) = OmegaW(i);
    end
end
OmegaW_dot = (OmegaWsat-OmegaW_1)/Ts;
Usat = -Jw*OmegaW_dot - cross(Omega, Jw*OmegaWsat) + Tau_gg;
end
%% ========================================================================
function [Usat]= Usaturate(U,Umax,Umin)
for i=1:3
    if U(i) > Umax
        Usat(i,1) = Umax;
    elseif U(i) < Umin
        Usat(i,1) = Umin;
    else
        Usat(i,1)=U(i);
    end
end
end