clc 
clear
close all

%% Data

%Apollo 11 re-entry initial conditions:
%vo = 11032;       %[m/s]
vo = 39461/3.6;
%ho = 400000*0.3048;   %[m]
ho = 121000;
gamma0 = deg2rad(6.49); 

%Apollo geometric data
%Izz = 4736*14.59390/0.09290304; %[kg/m^2]
Izz = 5273*1.355818;
m = 5560;    %[kg]

dm = 154*0.0254; %[m]
S = 0.25*pi*dm^2;
Nr = 184.8*0.0254;

% S = 2*pi*6*0.2728219713;       %[m^2] Àrea del casquet semiesfèric = 2piRh on h es alçada del casquet (0.27)
% Long = 107.5*0.3048;        %[m]
% Nr = 1*0.3048;           %[m]
% dm = 2.5;
% c = 0.756*dm;
% CD = 0.84;
% CL = 0.84;
% beta = W/(CD*S);         %[]
 

Re = 6371e3;  %[m]
dt = 1/50; 
n = round(1600/dt);

%% Coefficient data
alpha_cmdata = [0 60 70 130 180];
Ma_cmdata = [0.4 0.9 1.2 2.41 5 9];
Cmz_data = [0.04 -0.02 -0.085 0.04 -0.075;
            0.045 -0.01 -0.085 0.035 -0.07;
            0.06 -0.01 -0.03 0.015 -0.085;
            0.05 0.005 0 0.022 -0.08;
            0.045 0.005 0 0.025 -0.1;
            0.040 0 -0.005 0.025 -0.085];

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
 


alphas = 130:10:180
%alphas = [70 140];

% tf =1;
% 
% fprintf ('The calculations are made for 3 different initial flight path angles: \n \n \t 0.1º \n \t 1.0º \n \t 2.5º \n\n')
% prompt1 = 'Do you wish to introduce an additional initial flight path angle? Y/N: ';
% str = input(prompt1, 's');
% tf1 = ischar(str);
% 
% if tf1==1
% 
%     if str == 'Y' || str == 'y'
%         
%         prompt2 = '\nIntroduce an additional initial flight path angle ([º]): ';
%         aux = input(prompt2);
%         tf2 = ischar(aux);
%         if tf2==1
%            fprintf '\n\nFollow the instructions please :)\n Exiting... \n'
%            pause(3)
%            clc 
%            return
%         else
%         %gamma_extra = aux;
%         gammas = [0.1 1.0 2.5 aux]; 
%         end
%     elseif str == 'N'|| str == 'n'
% 
%         gammas = [0.1 1.0 2.5];
% 
%     else 
%         fprintf '\n\nFollow the instructions please :)\nExiting... \n'
%         pause(3)
%         clc
%         return
%     end
%     
% else
%      fprintf '\n\nFollow the instructions please :)\nExiting... \n'
%      pause(2)
%      clc
%      return
% end

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
    if(M(i-1))<0.4
        Cmz = interp2(alpha_cmdata,Ma_cmdata,Cmz_data,rad2deg(alpha(i-1)),0.4,'linear');
        CL = interp2(alpha_cldata,Ma_cldata,CL_data,rad2deg(alpha(i-1)),0.4,'linear');
        CD = interp2(alpha_cddata,Ma_cddata,CD_data,rad2deg(alpha(i-1)),0.4,'linear');
    end
    
%     
%     if (state==2 && i_count<= 500)
%         CL = CL - (2*CL)*(i-i_min)/500;
%     end
%     if (state==2 && i_count > 500 && i_count<3500)
%         CL = -CL;
%     end
%     
%     if (state==2 && i_count >= 3500 && i_count<4000)
%         CL = -CL + (2*CL)*(i_count-3500)/500;
%     end

    
   
     L_D = CL/CD;
     
     Moment = 1/2*rho*v(i-1)^2*S*dm*Cmz;
     alpha_dot(i) = alpha_dot(i-1) + (Moment/Izz)*dt;
     alpha(i) = alpha(i-1) + alpha_dot(i)*dt;
     if (alpha(i)<0 || alpha(i)>pi)
         break
     end
%      if (t(i)>60 && prod(alpha(1:i)<pi/2)==1) %burn-up cond
%          break
%      end
          
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
% 
% %figure(3)
% figure(2)
% subplot(2,2,1)
% %subplot(1,2,1)
% hold on
% plot(h/1000, Q,'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',rad2deg(alphas(1))))
% set(gca,'TickLabelInterpreter','latex');
% xlabel('Height$\;(km)$','Interpreter','latex')
% ylabel('Dynamic Pressure $\;(Pa)$','Interpreter','latex')
% legend('location','best','Interpreter','latex')
% xlim([0 inf])
% if j==size(alphas,2)
% grid  
% end
% 
% subplot(2,2,2)
% %subplot(1,2,2)
% hold on
% plot(h/1000, M,'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',rad2deg(alphas(1)))) 
% set(gca,'TickLabelInterpreter','latex');
% xlabel('Height$\;(km)$','Interpreter','latex')
% ylabel('Mach Number','Interpreter','latex')
% legend('location','best','Interpreter','latex')
% xlim([0 inf])
% if j==size(alphas,2)
% grid  
% end
% 
% %figure(4)
% subplot(2,2,3)
% %subplot(1,2,1)
% hold on
% plot(h/1000, Reynolds,'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',rad2deg(alphas(1)))) 
% set(gca,'TickLabelInterpreter','latex');
% xlabel('Height$\;(km)$','Interpreter','latex')
% ylabel('Reynolds Number','Interpreter','latex')
% legend('location','best','Interpreter','latex')
% xlim([0 inf])
% if j==size(alphas,2)
% grid  
% end
% 
% subplot(2,2,4)
% %subplot(1,2,2)
% hold on
% plot(h/1000, Pt,'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',rad2deg(alphas(1)))) 
% set(gca,'TickLabelInterpreter','latex');
% xlabel('Height$\;(km)$','Interpreter','latex')
% ylabel('Stagnation Point Pressure (Pa)','Interpreter','latex')
% legend('location','best','Interpreter','latex')
% xlim([0 inf])
% if j==size(alphas,2)
% grid  
% end
% 
% %figure(5)
% figure(3)
% subplot(2,2,1)
% %subplot(1,2,1)
% hold on
% plot(h/1000, ht,'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',rad2deg(alphas(1))))
% set(gca,'TickLabelInterpreter','latex');
% xlabel('Height$\;(km)$','Interpreter','latex')
% ylabel('Stagnation Point Enthalpy $\;(J/kg)$','Interpreter','latex')
% legend('location','best','Interpreter','latex')
% xlim([0 inf])
% if j==size(alphas,2)
% grid  
% end
% 
% subplot(2,2,2)
% %subplot(1,2,2)
% hold on
% plot(h/1000, -q,'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',rad2deg(alphas(1))))
% set(gca,'TickLabelInterpreter','latex');
% xlabel('Height$\;(km)$','Interpreter','latex')
% ylabel('Stagnation Point Heat Transfer $\;(W/m^2)$','Interpreter','latex')
% legend('location','best','Interpreter','latex')
% xlim([0 inf])
% if j==size(alphas,2)
% grid  
% end
% 
% %figure(6)
% subplot(2,2,3)
% %subplot(1,2,1)
% hold on
% plot(h/1000, Q.*v,'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',rad2deg(alphas(1))))
% set(gca,'TickLabelInterpreter','latex');
% xlabel('Height$\;(km)$','Interpreter','latex')
% ylabel('Dynamic Energy $\;(W/m^2)$','Interpreter','latex')
% legend('location','best','Interpreter','latex')
% xlim([0 inf])
% if j==size(alphas,2)
% grid  
% end

%subplot(2,2,4)
%subplot(1,2,1)
figure(2)
hold on
%plot(h/1000, rad2deg(alpha),'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',alphas(1)))
plot(h/1000,rad2deg(alpha),'DisplayName',sprintf('$\\alpha_0=%g^o$',rad2deg(alpha(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Angle of attack $\;(^o)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if j==size(alphas,2)
grid 
%xlim([8 121])
end

figure(7)
hold on
%plot(h/1000, rad2deg(alpha),'linestyle', linS{j},'DisplayName',sprintf('$\\alpha=%g$',alphas(1)))
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