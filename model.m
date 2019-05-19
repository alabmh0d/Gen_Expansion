% Done by Mohammed Alabdullah for Bulk Power Planning Term Paper
% May 2018
% 
% Some parts of the belwo code are based on the work done by Szilard Nement gamax.hu
%%
clear
clc
close all
% load plant properties table and load forecast
%load unitCommData
 %load load_profile %load profile
% load wind %net load profile with 8 GW wind
% load solar %net load profile with 8 GW solar pv
% load hybrid %net load profile with 4 GW wind and 4 GW solar pv
wind = xlsread('wind.xlsx');
solar = xlsread('solar.xlsx');
load_data = xlsread('demand.xlsx');

%% wind net load
wind_gen = 2*wind(:,1);
solar_gen = 7.5*solar(:,1);

net_load1=load_data-solar_gen-wind_gen;
net_load_norm1 = net_load1;%/max(net_load1);
B2 = reshape(net_load_norm1(1:8736),[168,52]);
W1 = mean([B2(:,1:6),B2(:,49:52)] ,2); % 10 weeks winter
W2 = mean(B2(:,7:17),2);% 11 weeks spring
W3 = mean(B2(:,18:37),2);% 20 weeks summer
W4 = mean(B2(:,38:48),2);% 11 weeks Automn

W = [W1;W2;W3;W4];
weight = [10 11 20 11];
h = 24*7;
ww = [repmat(weight(1),h,1);repmat(weight(2),h,1);repmat(weight(3),h,1);repmat(weight(4),h,1)];

load info2 %plants information
disp(info2)

%Cases: L (normal load), H (with hybrid renewables), P (PV), W (wind)
load_profile = W;%(400:672,:);

% we'll also create a few variables for time throughout the week
hours = numel(load_profile);
time = (1:hours)';


% add 5% reserver
GenerationTarget = load_profile * 1.05;

% pull apart plant properties
Fuel_Cost            = info2{1,:}+info2{11,:};% includes variable O&M Cost
Startup_Cost         = info2{2,:}/1000;
Operating_Cost       = info2{3,:}/1000;
MinGeneration_Level  = info2{4,:}.*info2{13,:}/1000;
MaxGeneration_Level  = info2{5,:}.*info2{13,:}/1000;
Ramp_Up              = info2{6,:}.*info2{13,:}/1000;
Ramp_Down            = info2{7,:}.*info2{13,:}/1000;
MinimumUp_Time       = info2{8,:};
MinimumDown_Time     = info2{9,:};

%candidate plants
plants2 = {'low_flex','med_flex','high_flex'};
Pmin = [100 50 10]/1000;%GW
Pmax = [200 200 200]/1000;%GW
Ramp=[70 100 150]/1000;%GW/hr
Inv=[2117 536 409]*3.75; % MSR/GW
FuelCost2=[18 250 3700]*3.75/1000; %MSR/GWh

nPlants = width(info2);
plants = info2.Properties.VariableNames;

InstalledCapacity   = sum(MaxGeneration_Level);

figure
hold on
plot(time,load_profile)
plot(time,GenerationTarget)
plot(xlim,[InstalledCapacity InstalledCapacity],'k:')
hold off
title('Load Profile')
ylabel('Load (GW)')
xlabel('Time (hrs)')
legend('Load Profile','Generation Target','Maximum Generation Level',...
    'Location','best')


%% here you need to modify info to match same type of plants

nSlots = hours*nPlants;
idxHr2ToEnd=2:hours;
%increase dimensions of any thing with a gnerator. 
maxGen_Const = repmat(MaxGeneration_Level,hours,1);
minGen_Const = repmat(MinGeneration_Level,hours,1);

maxGen_Const2 = repmat(Pmax,hours,1);
minGen_Const2 = repmat(Pmin,hours,1);

%Define the optimization problem and the optimization variables
powerprob = optimproblem;

% amount of power generated in an hour by a plant
power1 = optimvar('power',hours,plants,'LowerBound',0,'UpperBound',maxGen_Const);
% indicator if plant is operating during an hour 
isOn1 = optimvar('isOn',hours,plants,'Type','integer','LowerBound',0,'UpperBound',1);
% indicator if plant is starting up during an hour
startup1 = optimvar('startup',hours,plants,'Type','integer','LowerBound',0,'UpperBound',1);

% indicator if candidate plant to be invested 
ej = optimvar('ej',hours,plants2,'Type','integer','LowerBound',0,'UpperBound',1);
% amount of power generated in an hour by a added plant
power2 = optimvar('power2',hours,plants2,'LowerBound',0,'UpperBound',maxGen_Const2);


%% Define the objective function
% costs
powerCost1 = sum(power1*Fuel_Cost',1);
isOnCost1 = sum(isOn1*Operating_Cost',1);
startupCost1 = sum(startup1*Startup_Cost',1);
OC = powerCost1+isOnCost1+startupCost1;

AOC = sum(power2*FuelCost2',1);
ROI = 0.05;%rate of investment return
n = 1:30; % assume 30 years lifecycle of generators
investment = Pmax.*Inv*sum(1./(1+ROI).^n)^(-1);
AIC = sum(ej*investment');

% set objective
powerprob.Objective = OC;
%powerprob.Objective = OC+AOC+AIC;

%% constraints

powerprob.Constraints.isDemandMet = sum(power1,2) >= GenerationTarget;
%powerprob.Constraints.isDemandMet = sum(power1,2)+sum(power2,2) >= GenerationTarget;

%plant operating status to power generation
% only gen power when plant is on
% if isOn=0 power must be zero, if isOn=1 power must be less than maxGenConst
powerprob.Constraints.powerOnlyWhenOn = power1 <= maxGen_Const.*isOn1; 
% powerprob.Constraints.powerOnlyWhenOn2 = power2 <= maxGenConst2.*ej; 

% if on, meet MinGenerationLevel
% if isOn=0 power >= 0, if isOn=1 power must be more than minGenConst
powerprob.Constraints.meetMinGenLevel = power1 >= minGen_Const.*isOn1; 
% powerprob.Constraints.meetMinGenLevel2 = power2 >= minGenConst2.*ej; 

%enforce startup=1 when moving from off to on
% no need to enforce startup=0 at other times since minimizing cost forces it
powerprob.Constraints.startupConst = -isOn1(idxHr2ToEnd-1,:) + isOn1(idxHr2ToEnd,:) - startup1(idxHr2ToEnd,:) <= 0;
showconstr(powerprob.Constraints.startupConst(1:3))

%Ramprate limit constraints
% rampup limit
RampUpConst = repmat(Ramp_Up,hours-1,1);
% RampUpConst2 = repmat(Ramp,nHours-1,1);
powerprob.Constraints.rampupConst = -power1(idxHr2ToEnd-1,:) + power1(idxHr2ToEnd,:) <= RampUpConst(idxHr2ToEnd-1,:) + ...
     max(minGen_Const(idxHr2ToEnd,:)-RampUpConst(idxHr2ToEnd-1,:),0).*startup1(idxHr2ToEnd,:);
 
%  
%  powerprob.Constraints.rampupConst = -power1(idxHr2ToEnd-1,:)-power2(idxHr2ToEnd-1,:) + power1(idxHr2ToEnd,:)+power2(idxHr2ToEnd,:) <= RampUpConst(idxHr2ToEnd-1,:) + ...
%      max(minGenConst(idxHr2ToEnd,:)-RampUpConst(idxHr2ToEnd-1,:),0).*startup1(idxHr2ToEnd,:);
%  
 
% rampdown limit
RampDownConst = repmat(Ramp_Down,hours-1,1);
powerprob.Constraints.rampdownConst = power1(idxHr2ToEnd-1,:) - power1(idxHr2ToEnd,:) <= max(minGen_Const(idxHr2ToEnd,:),RampDownConst(idxHr2ToEnd-1,:)) - ...
    max(minGen_Const(idxHr2ToEnd,:)-RampDownConst(idxHr2ToEnd-1,:),0).*isOn1(idxHr2ToEnd,:);
showconstr(powerprob.Constraints.rampdownConst(1:3))

%Minimum uptime and downtime constraints
% min uptime
powerprob.Constraints.minUptimeConst = optimconstr(hours,plants);
for jj = 1:nPlants
    for kk = 1:hours % based on possible startups; no penalty at end for running over
        if kk > hours-MinimumUp_Time(jj)
            sumidx = kk:hours;
        else
            sumidx = kk:kk+MinimumUp_Time(jj)-1;
        end
        powerprob.Constraints.minUptimeConst(kk,jj) = ...
            startup1(kk,jj) - sum(isOn1(sumidx,jj)/length(sumidx)) <= 0;
    end
end

% min downtime
powerprob.Constraints.minDowntimeConst = optimconstr(hours,plants);
for jj = 1:nPlants
    for kk = 2:hours % based on possible startups; no penalty at beginning
        if kk <= MinimumDown_Time(jj)
            sumidx = 1:kk-1;
        else
            sumidx = kk-MinimumDown_Time(jj):kk-1;
        end
        powerprob.Constraints.minDowntimeConst(kk,jj) = ...
            startup1(kk,jj) + sum(isOn1(sumidx,jj)/length(sumidx)) <= 1;
    end
end
showconstr(powerprob.Constraints.minDowntimeConst(1:5,:))

%% Find the optimal generation 
% options for the optimization algorithm, here we set the max time it can run for
options = optimoptions('intlinprog','MaxTime',10);
% call the optimization solver to find the best solution
[sol,TotalCost,exitflag,output] = solve(powerprob,options);

%Reshape and Visualize Optimized Schedule Data
Scheduled_LoadProfile = sol.power;
loadPercentage = Scheduled_LoadProfile./repmat(MaxGeneration_Level,hours,1);
figure
%schedulePlot(loadPercentage',[],[],fliplr(plants));

weighted_cost = sum((ww'*Scheduled_LoadProfile).*Fuel_Cost)+sum((ww'*sol.isOn).*Operating_Cost)+sum((ww'*sol.startup).*Startup_Cost);

%Compare Scheduled Generation with Generation Target
figure;
hold on
plot(time,GenerationTarget)
plot(time,sum(Scheduled_LoadProfile,2))
plot(xlim,[InstalledCapacity InstalledCapacity],'k:')
hold off
ylabel('Load (MW)');
xlabel('Time (hrs)');
legend('Generation_Target','Scheduled_Generation','Maximum_Generation Level',...
    'Location','best');