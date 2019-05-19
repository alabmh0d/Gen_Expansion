close all;
clc;
clear all;
%load Saudi Demand data for 2015 (GW)
load_data = xlsread('demand.xlsx');
wind = xlsread('wind.xlsx');
solar = xlsread('solar.xlsx');

% matrix elements indices are arranged as follows:
% Tech/Fuel	CO      D       G       HFO
% CC        1       4       7       10
% GT        2       5       8       11
% ST        3       6       9       12

%Saudi Thermal Capacity 2015 (GW)   (don't forget IPP)   
capacity = [ 2.232	0	6.462	0;
    13.463	6.46	5.663	0;
    0	0	6.876	9.665];

%Heat Rate (TBTU/GWh)
heat_rate = [0.010306	nan	0.006527	nan;
    0.013281	0.013219	0.012663 nan;
    nan	nan	0.01	0.009];
%Fuel cost [New Saudi Cost (MSR/10^12*BTU)]
fuel_cost = [9.27;20.2;8;7.32];%[Crude Oil;Diesel;Gas;Heavy Fuel Oil]

%Fuel cost (MSR/GWh)		
gen_fuel_cost = zeros(size(heat_rate,1),size(heat_rate,2));
for i = 1:size(heat_rate,1)
    for j = 1:size(heat_rate,2)
        gen_fuel_cost (i,j) = heat_rate(i,j)*fuel_cost(j);
    end
end

%Conventional Variable O&M (MSR/GWh)				
O_M=[0.006	nan	0.007	nan;
    0.01	0.01	0.007	nan;
    nan	nan	0.002	0.003];

% Capex MSR/GW[CC;GT;ST;FPV;SPV;DPV;CSP;W;N;GEO];
capex = [3563;2719;6290;4375;4900;5425;21550;6113;24375;18438];

%Ramping Limits (p.u. of Capacity per Hour) 				
ramp =[0.6	nan	4.2	nan;
    4.2	7.2	6	nan;
    nan	nan	1.2	1.02];
t = 8760;

%Renewable Variable O&M  (MSR/GWh)	
% PV
% CSP
% W
% N
% GEO

RE_Var = [0
0.00041
0.01954
0.00195
0];

% Renweable Fixed O&M MSR / Gwy
% FPV
% SPV
% DPV
% W
% CSP
% N
% GEO

Re_fixed = [37.2
44.7
55.95
44.925
182.9625
351.6375
654.3375];



%% General Observation of demand data
% Peak demand info
[peak,I] = max(load_data);
peak_hour = mod(I,24);
peak_day = ceil(I/24);

%average Load , and load factor
average_load = mean(load_data);
load_factor = average_load./peak;
load_energy = sum(load_data);
%% LDC 
LDC=sortrows(load_data,-1);
LDC=LDC';
n=1:t;
[p,~,mu]=polyfit(n,LDC,4);
m=polyval(p,n,[],mu);
figure (1)
plot(n,LDC,'k-',n,m,'-r')
legend('LDC Curve','Fitted')
grid on
figure(2)
plot(n,LDC)
grid on

%% Operating costs of exisitng capacity
% need to get fuel costs plus O&M costs that satisfies G>L
% priority of dispatching to units with cheap fuel and O&M
capacity=capacity(:);

% conv_var_cost = gen_fuel_cost+O_M;
% [cheap_operation, index]=sortrows(conv_var_cost(:),1);
% gen = [];
% i = 0;
% inc = 0.25; %0.25 GW
% cost = 0;
% while gen_hr<load_energy
%     i = i +1;
%     gen(i) = gen+capacity(index(i));
%     imp = LDC<gen;
%     [max, ind] = max(imp);
%     (8760-ind)*gen(i)+sum(imp.*LDC);
%     
%         cost = cost + capacity(index(i))*cheap_operation(i)*t;
% 
% end

%% reducing load data
load_data_norm = load_data;%/max(load_data);
B1 = reshape(load_data_norm(1:8736),[168,52]);
S1 = mean([B1(:,1:6),B1(:,49:52)] ,2);
S2 = mean(B1(:,7:17),2);
S3 = mean(B1(:,18:32),2);
S4 = mean(B1(:,33:52),2);

L = [S1;S2;S3;S4];
t = 1:672;
figure,
plot(t',L)
title('Load Profile with 4 Representative Weeks')
ylabel('Normalized Load')
xlabel('Time (hrs)')
grid on 

%% 8 GW wind net load
wind_gen = wind(:,7);
net_load1=load_data-wind_gen;
net_load_norm1 = net_load1;%/max(net_load1);
B2 = reshape(net_load_norm1(1:8736),[168,52]);
W1 = mean([B2(:,1:6),B2(:,49:52)] ,2);
W2 = mean(B2(:,7:17),2);
W3 = mean(B2(:,18:32),2);
W4 = mean(B2(:,33:52),2);

W = [W1;W2;W3;W4];
t = 1:672;
figure,
plot(t',W)
title('Net Load with 8GW Wind')
ylabel('Normalized Load')
xlabel('Time (hrs)')
grid on

%% 8 GW solar net load
solar_gen = solar(:,7);
net_load2=load_data-solar_gen;
net_load_norm2 = net_load2;%/max(net_load2);
B3 = reshape(net_load_norm2(1:8736),[168,52]);
P1 = mean([B3(:,1:6),B3(:,49:52)] ,2);
P2 = mean(B3(:,7:17),2);
P3 = mean(B3(:,18:32),2);
P4 = mean(B3(:,33:52),2);

P = [P1;P2;P3;P4];
t = 1:672;
figure,
plot(t',P)
title('Net Load with 8GW Solar PV')
ylabel('Normalized Load')
xlabel('Time (hrs)')
grid on

%% 4 GW solar + 4 GW wind net load
hybrid = solar(:,6)+solar(:,5)+wind(:,6)+wind(:,5);
net_load3=load_data-hybrid;
net_load_norm3 = net_load3;%/max(net_load3);
B4 = reshape(net_load_norm3(1:8736),[168,52]);
H1 = mean([B4(:,1:6),B4(:,49:52)] ,2); % 10 weeks winter
H2 = mean(B4(:,7:17),2);% 11 weeks spring
H3 = mean(B4(:,18:37),2);% 20 weeks summer
H4 = mean(B4(:,38:48),2);% 11 weeks Automn

H = [H1;H2;H3;H4];
t = 1:672;
figure,
plot(t',H)
title('Net Load with 4 GW solar + 4 GW wind')
ylabel('Normalized Load')
xlabel('Time (hrs)')
grid on



