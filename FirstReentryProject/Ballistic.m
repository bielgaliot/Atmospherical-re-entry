clc
close all 
clear

%% Data

RE = 6371000;
h0 = 250000*0.3048; %m
delta_h = -1;
linS = {'-','--','-.', ':','-'}; 

v0 = 22500*0.3048; %m/s
gamma0 = 12*pi/180; %rad

n = round(h0/abs(delta_h));

h = linspace(h0,0,n);
v = zeros(1,n);
gamma = zeros(1,n);
t = zeros(1,n);
r = zeros(1,n);
Q = zeros(1,n);
M = zeros(1,n);
Pt = zeros(1,n);
ht = zeros(1,n);
Re = zeros(1,n);
q = zeros(1,n);
E = zeros(1,n);

v(1) = v0;
gamma(1) = gamma0;
t(1) = 0;
r(1) = 0;
[T,P,rho,g,mu] = isa_calc(h(1));
%[rho,~,T,P,~,~,~] = atmos(h(1));
Q(1) = (1/2)*v(1)^2*rho;
M(1) = v(1)/sqrt(1.4*287*T);
%Pt(1) = P+(1/2)*rho*v(1)^2;

        p_shock = P*(1+(2.8/2.4)*(M(1)^2-1));
        Pt(1) = p_shock*(1+0.2*((1+0.2*M(1)^2)/(1.4*M(1)^2-0.2)))^(1.4/0.4);

%Cp de l'aire en funció de T:
%1034.09-0.2849*T+7.817e-4*T^2-4.971e-7*T^3+1.077e-10*T^4
cp = 1034.09-0.2849*T+7.817e-4*T^2-4.971e-7*T^3+1.077e-10*T^4;
ht(1) = cp*(T+1/(2*cp)*v(1)^2);
% viscositat de l'aire (mu) en funció de T: (1.458e-6*T^(1.5))/(T+110.4)
Re(1) = rho*v(1)*0.3048/mu;
q(1) = 1.83e-4*(v(1))^3*sqrt(rho/0.3048);
E(1) = Q(1)*v(1);

fprintf ('The calculations are made for 4 different values of beta (W/(cd*A) [lbf/ft^2]): \n \n \t 100 \n \t 500 \n \t 1000 \n \t 5000 \n\n')
prompt1 = 'Do you wish to introduce an additional value for beta ([lbf/ft^2])? Y/N: ';
str = input(prompt1, 's');
tf1 = ischar(str);

if tf1==1

    if str == 'Y' || str == 'y'
        
        prompt2 = '\nIntroduce an additional value for beta (lbf/ft^2): ';
        aux = input(prompt2);
        tf2 = ischar(aux);
       
        if tf2==1
           fprintf '\n\nFollow the instructions please :)\n Exiting... \n'
           pause(3)
           clc 
           return
        else
        fprintf '\n\n*Note that units for the BC in the legend will be in the International System\n'
        betas = [100 500 1000 5000 aux]; 
        end
    elseif str == 'N'|| str == 'n'
        fprintf '\n\n*Note that units for the BC in the legend will be in the International System\n'
        betas = [100 500 1000 5000];

    else 
        fprintf '\n\nFollow the instructions please :)\nExiting... \n'
        pause(2.5)
        clc
        return
    end
    
else
     fprintf '\n\nFollow the instructions please :)\nExiting... \n'
     pause(2.5)
     clc
     return
end

%% Calculations
for i=1:size(betas,2)

beta = betas(i)*4.4482216/0.3048^2;

for j=2:n
    [T,P,rho,g,mu] = isa_calc(h(j));
    v(j) = v(j-1) + (g*rho*v(j-1)/(2*beta*sin(gamma(j-1)))-g/v(j-1))*delta_h;
    gamma(j) = gamma(j-1) + (1/(RE+h(j))*1/tan(gamma(j-1))-g/((v(j-1))^2*tan(gamma(j-1))))*delta_h;
    t(j) = t(j-1)-1/(v(j-1)*sin(gamma(j-1)))*delta_h;
    r(j) = r(j-1)-RE/(RE+h(j-1))*1/tan(gamma(j-1))*delta_h;
    Q(j) = (1/2)*v(j)^2*rho;
    cp = 1034.09-0.2849*T+7.817e-4*T^2-4.971e-7*T^3+1.077e-10*T^4;
    M(j) = v(j)/sqrt(1.4*287*T);
    if M(j)>=1
        p_shock = P*(1+(2.8/2.4)*(M(j)^2-1));
        Pt(j) = p_shock*(1+0.2*((1+0.2*M(j)^2)/(1.4*M(j)^2-0.2)))^(1.4/0.4);
    end
    if (M(j)>=0.3 && M(j)<1)
        Pt(j) = P*(1+0.2*M(j)^2)^(1.4/0.4);
    end
    if (M(j)<0.3)
        Pt(j) = P+(1/2)*rho*v(j)^2;
    end
    ht(j) = cp*(T+1/(2*cp)*v(j)^2);
    Re(j) = rho*v(j)*0.3048/mu;
    q(j) = 1.83e-4*(v(j))^3*sqrt(rho/0.3048);
    E(j) = Q(j)*v(j);
end

%% Plots

figure(1)
subplot(2,2,1)
hold on
plot(h/1000,v,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle',linS{i}) %per ordre: beta actual a la llegenda, gruix de línia i tipus de línia
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Velocity $\;(m/s)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor  %només s'ha d'activar a l'última iteració pq si no es va activant i desactivant cada vegada
end


a = diff(v)./diff(t);
subplot(2,2,2)
hold on
plot(h(1:n-1)/1000,-a/9.81,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Deceleration $\;(g)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

subplot(2,2,3)
hold on
plot(h/1000,Q,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Dynamic pressure $\;(Pa)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

subplot(2,2,4)
hold on
plot(h/1000,M,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Mach number','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end


figure(2)
subplot(2,2,1)
hold on
plot(h/1000,Re,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Reynolds number $\;(1/m)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

subplot(2,2,2)
hold on
plot(h/1000,Pt,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Pressure $\;(Pa)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

subplot(2,2,3)
hold on
plot(h/1000,ht,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Enthalpy $\;(J/kg)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

subplot(2,2,4)
hold on
plot(h/1000,-q,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Stagnation point heat transfer $\;(W/m^2)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

figure(3)
subplot(2,2,1)
hold on
plot(h/1000,t,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Entry time $\;(s)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

subplot(2,2,2)
hold on
plot(h/1000,r/1000,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Range $\;(km)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

subplot(2,2,3)
hold on
plot(h/1000,E,'DisplayName',sprintf('$\\beta=%g$',beta),'linew',1,'linestyle', linS{i})
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Dynamic energy $\;(W/m^2)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
if i==size(betas,2)
    grid on
    grid minor
end

end