close all 

RE = 6371000;
h0 = 250000*0.3048; %m
delta_h = -1;

v0 = 23000*0.3048; %m/s
gamma0 = 1.5*pi/180; %rad

W = 2662.8*4.44822; %N
A = 30.27*0.3048^2; %m^2
Cd = 1.6;
beta = W/(A*Cd);
Long = 6.2*0.3048; %m

n = round(h0/abs(delta_h));

h = linspace(h0,0,n);
v = zeros(1,n);
gamma = zeros(1,n);
t = zeros(1,n);
r = zeros(1,n);
M = zeros(1,n);
Q = zeros(1,n);
q = zeros(1,n);
E = zeros(1,n);
Re = zeros(1,n);
Pt = zeros(1,n);
ht = zeros(1,n);

v(1) = v0;
gamma(1) = gamma0;
t(1) = 0;
r(1) = 0;
[T,P,rho,g,mu] = isa_calc(h(1));
Q(1) = (1/2)*v(1)^2*rho;
M(1) = v(1)/sqrt(1.4*287*T);
p_shock = P*(1+(2.8/2.4)*(M(1)^2-1));
Pt(1) = p_shock*(1+0.2*((1+0.2*M(1)^2)/(1.4*M(1)^2-0.2)))^(1.4/0.4);
cp = 1034.09-0.2849*T+7.817e-4*T^2-4.971e-7*T^3+1.077e-10*T^4;
ht(1) = cp*(T+1/(2*cp)*v(1)^2);
Re(1) = rho*v(1)*Long/mu;
q(1) = 1.83e-4*(v(1))^3*sqrt(rho/0.3048);
E(1) = Q(1)*v(1);

for j=2:n
    [T,P,rho,g,mu] = isa_calc(h(j));
    v(j) = v(j-1) + (g*rho*v(j-1)/(2*beta*sin(gamma(j-1)))-g/v(j-1))*delta_h;
    gamma(j) = gamma(j-1) + (1/(RE+h(j))*1/tan(gamma(j-1))-g/((v(j-1))^2*tan(gamma(j-1))))*delta_h;
    t(j) = t(j-1)-1/(v(j-1)*sin(gamma(j-1)))*delta_h;
    r(j) = r(j-1)-RE/(RE+h(j-1))*1/tan(gamma(j-1))*delta_h;
    
    M(j) = v(j)/sqrt(1.4*287*T);
    Q(j) = (1/2)*v(j)^2*rho;
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
    Re(j) = rho*v(j)*Long/mu;
    q(j) = 1.83e-4*(v(j))^3*sqrt(rho/0.3048);
    E(j) = Q(j)*v(j);
    ht(j) = cp*(T+1/(2*cp)*v(j)^2);
end

figure (1)
subplot(2,2,1)
hold on
plot(h/1000,v,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Velocity $\;(m/s)$','Interpreter','latex')
grid on
grid minor

subplot(2,2,2)
hold on
plot(h/1000,t,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Entry time $\;(s)$','Interpreter','latex')
grid on
grid minor

subplot(2,2,3)
hold on
plot(h/1000,r/1000,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Range $\;(km)$','Interpreter','latex')
grid on
grid minor

a = diff(v)./diff(t);
subplot(2,2,4)
hold on
plot(h(2:size(h,2))/1000,-a/9.81,'linew',1)
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Deceleration $\;(g)$','Interpreter','latex')
grid on
grid minor

figure (2)
subplot(2,2,1)
hold on
plot(h/1000,M,'linew',1)
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Mach number','Interpreter','latex')
grid on
grid minor

subplot(2,2,2)
hold on
plot(h/1000,Q,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Dynamic Pressure $\;(Pa)$','Interpreter','latex')
grid on
grid minor

subplot(2,2,3)
hold on
plot(h/1000,Pt,'linew',1)
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Pressure $\;(Pa)$','Interpreter','latex')
grid on
grid minor

subplot(2,2,4)
hold on
plot(h/1000,Re,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Reynolds number','Interpreter','latex')
grid on
grid minor

figure (3)
subplot(2,2,1)
hold on
plot(h/1000,-q,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Heat Transfer $\;(W/m^2)$','Interpreter','latex')
grid on
grid minor

subplot(2,2,2)
hold on
plot(h/1000,E,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Dynamic Energy $\;(W/m^2)$','Interpreter','latex')
grid on
grid minor

subplot(2,2,3)
hold on
plot(h/1000,ht,'linew',1) 
set(gca,'TickLabelInterpreter','latex');
xlabel('Altitude $\;(km)$','Interpreter','latex')
ylabel('Stagnation Point Enthalpy $\;(J/kg)$','Interpreter','latex')
grid on
grid minor