function [item_final, vendor, act_torque, T_us, RPM_us, min_dist,act_vol,act_mass,act_NLRPM,act_price,act_voltage,act_type] =...
    euclid_dist_fnctn2(min_rpm, min_torque, ideal_vol, max_vol, Input_Voltage,Price_Lim,Mass_Lim)
% [item_final, price_final, vol_final] = euclid_dist_fnctn(Input_RPM, Input_Torque, Input_Volume)
%
%  Input_RPM = 900 ;
% %  Input_St_Current = 1.3;
% Input_Torque = 4.8 ;
% Input_Volume = 300;
% Input_Mass = 400 ;
% Input_Price = 100 ;
%
% a = [Input_Torque,Input_RPM,Input_Volume] ;

Input_RPM = min_rpm;
Input_Torque = min_torque;
Input_Volume = ideal_vol;
%Price_Lim = 100; % [USD]
%Price_Lim = 100000;
%Mass_Lim = 200; % [g]
%Mass_Lim = 2000;
% min_rpm = 0;
% min_torque = 8;
% max_vol = 4000;

b_lim = [min_rpm, min_torque, max_vol];

%%

[num,txt,raw] = xlsread('TSADB', 'Sheet1', 'A2:N326');
RPM = num(:,8);
% St_Current = num(:,10);
Torque = num(:,11);
Volume = num(:,6);
Voltage = num(:,12);
Mass = num(:,7);
Price = num(:,13);
T_us = Torque;
RPM_us= RPM;

% figure
% scatter3(Torque, RPM, Volume);
% xlabel('$\tau_s$ [kg$\cdot$cm]');
% xlim([0, 1000]);
% ylabel('$\omega_{NL}$ [RPM]');
% ylim([0, 1000])
% zlabel('Volume, $\Gamma$ [mm$^3$]');
% zlim([0, 2e5]);
% 
% figure
% plot(Torque, RPM,'o')
% xlim([0,5])
% ylim([0,3000])
%% mean
mean_Torque = mean(Torque','omitnan');
mean_RPM = mean(RPM','omitnan');
% mean_St_Current = mean(St_Current','omitnan')
mean_Volume = mean(Volume','omitnan');
% mean_Mass = mean(Mass','omitnan')
% mean_Price = mean(Price','omitnan')
%%
%%standard deviation
sd_Torque = std(Torque','omitnan');
sd_RPM = std(RPM','omitnan');
% sd_St_Current = std(St_Current','omitnan')
sd_Volume = std(Volume','omitnan');
% sd_Mass = std(Mass','omitnan')
% sd_Price = std(Price','omitnan')

%%make the NAn to zeros
Torque(isnan(Torque))=0;
RPM(isnan(RPM))=0;
% St_Current(isnan(St_Current))=0;
Volume(isnan(Volume))=0;
% Mass(isnan(Mass))=0;
% Price(isnan(Price))=0;

%%
b_us =  [RPM, Torque, Volume];
%scaling the data
Torque = (Torque-mean_Torque)/ sd_Torque ;
RPM = (RPM-mean_RPM)/ sd_RPM ;
% St_Current = (St_Current-mean_St_Current)/ sd_St_Current ;
Volume = (Volume-mean_Volume)/ sd_Volume ;
indices_nan_Torque = find(isnan(Torque));
indices_nan_RPM = find(isnan(RPM));
%indices_nan_St_Current = find(isnan(St_Current) == 1);
indices_nan_Volume = find(isnan(Volume));
% indices_nan_Mass = find(isnan(Mass) == 1);
% indices_nan_Price = find(isnan(Price) == 1);
a_norm= [(Input_RPM-mean_RPM)/sd_RPM, (Input_Torque-mean_Torque)/ sd_Torque,...
    (Input_Volume-mean_Volume)/sd_Volume] ;

%%
k = 1;
num_rec = 5; % number of recommendations
iteration = size(RPM,1);
distance = zeros(iteration,1)+Inf;
% distance = 0;
for i = 1: iteration
    b = [RPM(i), Torque(i), Volume(i)];
    % distance(i) = norm(a_norm-b);
    
    % if the actual values are within the limits
    % b_lim(1) = lower limit
    % b_lim(2) = lower limit
    % b_lim(3) = upper limit
    if (b_lim(1) <= b_us(i,1) && b_lim(2) <= b_us(i,2) && b_lim(3) >= b_us(i,3))
        if Price(i) <= Price_Lim && Mass(i) <= Mass_Lim
            if logical(Input_Voltage)
                if Voltage(i) == Input_Voltage% also check the voltage
                % that's good, save that value for the recommendation
                    distance(i) = norm(a_norm-b);
                % disp('In Bounds')
                % k = k + 1;
                end
            else
                distance(i) = norm(a_norm-b);
            end
        end
    else
        % disp('Out of bounds')
    end
end
    item_final = cell(1,num_rec);
    vendor = cell(1,num_rec);
    act_torque = cell(1,num_rec);
    act_vol = cell(1,num_rec);
    act_mass = cell(1,num_rec);
    act_NLRPM = cell(1,num_rec);
    act_price = cell(1,num_rec);
    act_voltage = cell(1,num_rec);
    act_type = cell(1,num_rec);
if (sum(~(isinf(distance)))) % If there's at least one non-infinite entry.
    if (sum(~(isinf(distance))) < num_rec)
        % Adjust the number of recommendations based on the number of
        % non-infinite entries.
        num_rec = sum(~(isinf(distance)));
    end
    [min_dist,min_ind] = mink(distance, num_rec);
    for i=1:num_rec
        item_final{i} = raw{min_ind(i),2}; % 2 = Product Number
        vendor{i} = raw{min_ind(i),1}; % 1 = Vendor
        act_vol{i} = raw{min_ind(i),7}; % 7 = Boxed Volume
        act_mass{i} = raw{min_ind(i),8}; % 8 = mass
        act_NLRPM{i} = raw{min_ind(i),9}; % 9 = no-load RPM
        act_torque{i} = raw{min_ind(i), 12}; % 12 = Stall Torque
        act_voltage{i} = raw{min_ind(i),13}; % 13 = actual voltage
        act_price{i} = raw{min_ind(i),14}; % 14 = price
        act_type{i} = raw{min_ind(i),3}; % type of motor
    end
else
    min_dist = zeros(1,num_rec);
    min_ind = zeros(1,num_rec);
    for i = 1:num_rec
        min_dist(i) = Inf;
        item_final{i} = '0';
        vendor{i} = '0';
        act_vol{i} = '0';
        act_mass{i} = '0';
        act_NLRPM{i} = '0';
        act_torque{i} = '0';
        act_voltage{i} = 0;
        act_price{i} = '0';
        act_type{i} = '0';
    end
end
end
