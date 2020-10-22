clc 
clear
close all

%% Data

%Apollo 11 re-entry initial conditions:
vo = 39461/3.6;
ho = 121000;
gamma0 = deg2rad(6.49); 

%Apollo geometric data
Izz = 5273*1.355818;
m = 5560;    %[kg]

dm = 154*0.0254; %[m]
S = 0.25*pi*dm^2;
Nr = 184.8*0.0254;

Re = 6371e3;  %[m]
dt = 1/50; 
n = round(1600/dt);

%% Coefficient data

% PITCHING MOMENT DATA FROM ZHANG ET. AL. 

alpha_cmdata2 = [-180 -170 -160 -150 -140 -130 -120 -110 -100 -90 -70 -50 -45 -40 -35 -30 -25 -20 -15 -10 0];
for i = 1:size(alpha_cmdata2,2)
    alpha_cmdata(i) = 180 + alpha_cmdata2(i)
end

Ma_cmdata = [0.6 1.5 2 5 10];
Cmz_data = [0 0.043 0.08 0.085 0.07 0.06 0.059 0.033 0.01 0.018 0.047 0.067 0.071 0.075 0.077 0.07 0.05 0.025 0.01 -0.008 -0.03;
            0.012 0.028 0.048 0.06 0.063 0.065 0.059 0.048 0.04 0.032 0.034 0.032 0.035 0.039 0.037 0.031 0.02 0.01 -0.01 -0.025 -0.045;
            0.01 0.018 0.03 0.045 0.056 0.059 0.054 0.045 0.03 0.028 0.033 0.036 0.037 0.033 0.028 0.02 0.013 0 -0.012 -0.022 -0.045;
            0.01 0.013 0.02 0.03 0.045 0.048 0.041 0.03 0.024 0.022 0.033 0.042 0.039 0.03 0.021 0.013 0.005 -0.006 -0.016 -0.028 -0.045;
            0.01 0.013 0.02 0.03 0.045 0.048 0.041 0.03 0.024 0.022 0.033 0.042 0.039 0.03 0.021 0.013 0.005 -0.006 -0.016 -0.028 -0.045];


alpha_cddata = [0 30 60 90 120 150 180];
Ma_cddata = [0.4 0.9 1.1 2.49 5 9];
CD_data = [0.6 0.55 0.4 0.3 0.4 0.75 0.97;
           0.75 0.65 0.65 0.45 0.5 0.95 1.1;
           1 0.8 0.8 0.65 0.75 1.1 1.3;
           0.9 0.95 0.9 0.6 0.6 1.2 1.5;
           0.2 0.4 0.38 0.01 0.01 0.65 1;
           0.2 0.4 0.38 0.05 0.01 0.6 1];
       
alpha_cldata = [0 20 60 80 140 180];
Ma_cldata = [0.4 0.9 1.1 2.49 5 9];
CL_data = [0 0.2 0.5 -0.2 1.3 0;
           0 0.1 0.2 -1 1.2  0;
           0 -0.1 0 -0.9 1.7 0;
           0 0.1 -0.15 -0.25 0.5 0;
           0 0.2 -0.15 -0.25 0.5 0;
           0 0.22 -0.15 -0.25 0.5 0];
 


alphas = 10:20:170


%% Calculations



for j=1:size(alphas,2)
    
    %% Vector Pre-initialization
linS = {'-','--','-.',':'}; 
t = linspace(0,1600,n);
v = zeros(1,n);
gamma = zeros(1,n);
h = zeros(1,n);
r = zeros(1,n);
Q = zeros(1,n);
M = zeros(1,n);
Reynolds = zeros(1,n);
Pt = zeros(1,n);
ht = zeros(1,n);
q = zeros(1,n); 
alpha = zeros(1,n);
alpha_dot = zeros(1,n);


%% Initial Conditions
%alpha(1)=-30;
h(1) = ho;
[T,P,rho,g_h,mu] = isa_calc(h(1));
v(1)=vo;
M(1) = v(1)/sqrt(1.4*287*T);
Reynolds(1) = rho*v(1)*dm/mu;
cp = 1034.09-0.2849*T+7.817e-4*T^2-4.971e-7*T^3+1.077e-10*T^4;
ht(1) = cp*T +v(1)^2/2;
q(1) = 1.83e-4*v(1)^3*sqrt(rho/Nr);
Q(1) = 1/2*rho*v(1)^2;
p_shock = P*(1+(2.8/2.4)*(M(1)^2-1));
Pt(1) = p_shock*(1+0.2*((1+0.2*M(1)^2)/(1.4*M(1)^2-0.2)))^(1.4/0.4);

    
    gamma(1) =  gamma0;
    alpha(1) =  deg2rad(alphas(j));
    alpha_dot(1) = 0;
    state = 1;
    i_count = 1;
    
for i=2:n
    [T,P,rho,g_h,mu] = isa_calc(h(i-1));
    
    if(M(i-1)<=9)
        Cmz = interp2(alpha_cmdata,Ma_cmdata,Cmz_data,rad2deg(alpha(i-1)),M(i-1),'linear');
        CL = interp2(alpha_cldata,Ma_cldata,CL_data,rad2deg(alpha(i-1)),M(i-1),'linear');
        CD = interp2(alpha_cddata,Ma_cddata,CD_data,rad2deg(alpha(i-1)),M(i-1),'linear');
    end
    if(M(i-1)>9)
        Cmz = interp2(alpha_cmdata,Ma_cmdata,Cmz_data,rad2deg(alpha(i-1)),9,'linear');
        CL = interp2(alpha_cldata,Ma_cldata,CL_data,rad2deg(alpha(i-1)),9,'linear');
        CD = interp2(alpha_cddata,Ma_cddata,CD_data,rad2deg(alpha(i-1)),9,'linear');
    end
    if(M(i-1))<0.6
        Cmz = interp2(alpha_cmdata,Ma_cmdata,Cmz_data,rad2deg(alpha(i-1)),0.6,'linear');
        CL = interp2(alpha_cldata,Ma_cldata,CL_data,rad2deg(alpha(i-1)),0.6,'linear');
        CD = interp2(alpha_cddata,Ma_cddata,CD_data,rad2deg(alpha(i-1)),0.6,'linear');
    end
    
    
    if (state==2 && i_count<= 500)
        CL = CL - (2*CL)*(i-i_min)/500;
    end
    if (state==2 && i_count > 500 && i_count<3500)
        CL = -CL;
    end
    
    if (state==2 && i_count >= 3500 && i_count<4000)
        CL = -CL + (2*CL)*(i_count-3500)/500;
    end

    
   
     L_D = CL/CD;
     
     Moment = 1/2*rho*v(i-1)^2*S*dm*Cmz;
     alpha_dot(i) = alpha_dot(i-1) + (Moment/Izz)*dt;
     alpha(i) = alpha(i-1) + alpha_dot(i)*dt;
     if (alpha(i)<0 || alpha(i)>pi)
         break
     end
     if (t(i)>60 && prod(alpha(1:i)<pi/2)==1) %burn-up cond
         break
     end
          
    Qd = 1/2*rho*v(i-1)^2;
    beta = m*g_h/(CD*S);
        
    v(i) = v(i-1) +  g_h*(-Qd/beta + sin(gamma(i-1)))*dt;
    h(i) = h(i-1) - v(i-1) * sin(gamma(i-1))*dt;
    r(i) = r(i-1) +  ((Re*v(i-1)*cos(gamma(i-1)))/(Re+h(i-1)))*dt;
    gamma(i) = gamma(i-1) + (((-Qd*g_h)/(beta)*L_D + cos(gamma(i-1))*(g_h-((v(i-1)^2)/(Re+h(i-1)))))/v(i-1))*dt;
    
    if( h(i)<8000 || h(i)>200e3)
        break
    end
    
    if (h(i)>h(i-1))
        state = 2;
        i_min = i;
        i_count = i_count+1;
    end
    
    Q(i) = Qd;
    M(i) = v(i)/sqrt(1.4*287*T);
    
    Reynolds(i) = rho*v(i)*dm/mu;
    if M(i)>=1
        p_shock = P*(1+(2.8/2.4)*(M(i)^2-1));
        Pt(i) = p_shock*(1+0.2*((1+0.2*M(i)^2)/(1.4*M(i)^2-0.2)))^(1.4/0.4);
    end
    if (M(i)>=0.3 && M(i)<1)
        Pt(i) = P*(1+0.2*M(i)^2)^(1.4/0.4);
    end
    if (M(i)<0.3)
        Pt(i) = P+(1/2)*rho*v(i)^2; 
    end
    cp = 1034.09-0.2849*T+7.817e-4*T^2-4.971e-7*T^3+1.077e-10*T^4;
    ht(i) = cp*T +v(i)^2/2;
    q(i) = 1.83e-4*(v(i))^3*sqrt(rho/Nr);
end

t(t==0) = NaN;
r(r==0) = NaN;
h(h==0) = NaN;
v(v==0) = NaN;
alpha(alpha==0) = NaN;
gamma(gamma==0) = NaN;

%% Plots

 figure(4)

hold on
plot(h/1000,v,'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alpha(1)))) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Velocity $\;(m/s)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if j==size(alphas,2)
grid  
xlim([8 121])
end
% 

figure(1)
hold on
plot(h/1000,t,'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alpha(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Entry Time$\;(m/s)$','Interpreter','latex')
legend('location','best','Interpreter','latex')

if j==size(alphas,2)
grid  
xlim([8 121])
end
% 
figure(3)

hold on
plot(h/1000,r/1000,'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alpha(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Range$\;(km)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
%xlim([0 inf])
if j==size(alphas,2)
    xlim([8 121])
grid  
end
% 

a = diff(v)./diff(t);

figure(5)
hold on
plot(h(1:i-1)/1000,-a(1:i-1)/9.81,'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alphas(1)))) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Deceleration$\;(g)$','Interpreter','latex')
legend('location','best','Interpreter','latex')

if j==size(alphas,2)
grid  
xlim([8 121])
end
figure(2)
hold on
plot(h/1000,rad2deg(alpha),'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alpha(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Angle of attack $\;(^o)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if j==size(alphas,2)
grid 
xlim([8 121])
end

figure(7)
hold on
plot(t,rad2deg(alpha),'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alpha(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Entry time $\;(s)$','Interpreter','latex')
ylabel('Angle of attack $\;(^o)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if j==size(alphas,2)
grid 
xlim([0 t(i)])
end



figure(6)
hold on
plot(h/1000,rad2deg(gamma),'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alpha(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Flight path angle $\gamma$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if j==size(alphas,2)
grid
xlim([8 121])
end

end