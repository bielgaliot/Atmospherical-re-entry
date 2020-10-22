clc 
clear
close all

%% Data
vo = 23000*0.3048;       %[m/s]
ho = 250000*0.3048;      %[m/s]
W = 200000*4.4482216;    %[N]
S = 2690*0.092903;       %[m^2]
Long = 107.5*0.3048;        %[m]
Nr = 1*0.3048;           %[m]
Re = 6371e3;             %[m]
CD = 0.84;
CL = 0.84;
beta = W/(CD*S);         %[]

dt = 1/100; 
n = round(1600/dt);

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

%% Initial Conditions
h(1) = ho;
[T,P,rho,g_h,mu] = isa_calc(h(1));
v(1)=vo;
M(1) = v(1)/sqrt(1.4*287*T);
Reynolds(1) = rho*v(1)*Long/mu;
cp = 1034.09-0.2849*T+7.817e-4*T^2-4.971e-7*T^3+1.077e-10*T^4;
ht(1) = cp*T +v(1)^2/2;
q(1) = 1.83e-4*v(1)^3*sqrt(rho/Nr);
Q(1) = 1/2*rho*v(1)^2;
p_shock = P*(1+(2.8/2.4)*(M(1)^2-1));
Pt(1) = p_shock*(1+0.2*((1+0.2*M(1)^2)/(1.4*M(1)^2-0.2)))^(1.4/0.4);

tf =1;

fprintf ('The calculations are made for 3 different initial flight path angles: \n \n \t 0.1º \n \t 1.0º \n \t 2.5º \n\n')
prompt1 = 'Do you wish to introduce an additional initial flight path angle? Y/N: ';
str = input(prompt1, 's');
tf1 = ischar(str);

if tf1==1

    if str == 'Y' || str == 'y'
        
        prompt2 = '\nIntroduce an additional initial flight path angle ([º]): ';
        aux = input(prompt2);
        tf2 = ischar(aux);
        if tf2==1
           fprintf '\n\nFollow the instructions please :)\n Exiting... \n'
           pause(3)
           clc 
           return
        else
        %gamma_extra = aux;
        gammas = [0.1 1.0 2.5 aux]; 
        end
    elseif str == 'N'|| str == 'n'

        gammas = [0.1 1.0 2.5];

    else 
        fprintf '\n\nFollow the instructions please :)\nExiting... \n'
        pause(3)
        clc
        return
    end
    
else
     fprintf '\n\nFollow the instructions please :)\nExiting... \n'
     pause(2)
     clc
     return
end

%% Calculations

for j=1:size(gammas,2)

    gamma(1) =  deg2rad(gammas(j));
    
for i=2:n
    
    [T,P,rho,g_h,mu] = isa_calc(h(i-1));
    Qd = 1/2*rho*v(i-1)^2;
    L = 1/2*rho*v(i-1)^2*S*CL;
    D = 1/2*rho*v(i-1)^2*S*CD;
    v(i) = v(i-1) +  g_h*(-Qd/beta + sin(gamma(i-1)))*dt;
    h(i) = h(i-1) - v(i-1) * sin(gamma(i-1))*dt;
    r(i) = r(i-1) +  ((Re*v(i-1)*cos(gamma(i-1)))/(Re+h(i-1)))*dt;
    gamma(i) = gamma(i-1) + (((-Qd*g_h*L)/(beta*D) + cos(gamma(i-1))*(g_h-((v(i-1)^2)/(Re+h(i-1)))))/v(i-1))*dt;
    
    Q(i) = Qd;
    M(i) = v(i)/sqrt(1.4*287*T);
    Reynolds(i) = rho*v(i)*Long/mu;
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

%% Plots

figure(1)
subplot(2,2,1)
%subplot(1,2,1)
hold on
plot(h/1000,v,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1)))) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Velocity $\;(m/s)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 80])
if j==size(gammas,2)
grid  
end

subplot(2,2,2)
%subplot(1,2,2)
hold on
plot(h/1000,t,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Entry Time$\;(m/s)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

%figure(2)
subplot(2,2,3)
%subplot(1,2,1)
hold on
plot(h/1000,r/1000,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Range$\;(km)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

%subplot(1,2,2)
a = diff(v)./diff(t);
subplot(2,2,4)

hold on
plot(h(1:n-1)/1000,-a/9.81,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1)))) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Deceleration$\;(g)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

%figure(3)
figure(2)
subplot(2,2,1)
%subplot(1,2,1)
hold on
plot(h/1000, Q,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Dynamic Pressure $\;(Pa)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

subplot(2,2,2)
%subplot(1,2,2)
hold on
plot(h/1000, M,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1)))) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Mach Number','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

%figure(4)
subplot(2,2,3)
%subplot(1,2,1)
hold on
plot(h/1000, Reynolds,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1)))) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Reynolds Number','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

subplot(2,2,4)
%subplot(1,2,2)
hold on
plot(h/1000, Pt,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1)))) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Pressure (Pa)','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

%figure(5)
figure(3)
subplot(2,2,1)
%subplot(1,2,1)
hold on
plot(h/1000, ht,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Enthalpy $\;(J/kg)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

subplot(2,2,2)
%subplot(1,2,2)
hold on
plot(h/1000, -q,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Heat Transfer $\;(W/m^2)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

%figure(6)
subplot(2,2,3)
%subplot(1,2,1)
hold on
plot(h/1000, Q.*v,'linestyle', linS{j},'DisplayName',sprintf('$\\gamma=%g$',rad2deg(gamma(1))))
set(gca,'TickLabelInterpreter','latex');
xlabel('Height$\;(km)$','Interpreter','latex')
ylabel('Dynamic Energy $\;(W/m^2)$','Interpreter','latex')
legend('location','best','Interpreter','latex')
xlim([0 inf])
if j==size(gammas,2)
grid  
end

end