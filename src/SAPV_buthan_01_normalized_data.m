%% This script is a mathemathical model of the simulink-version of the microgrid.
% It is very simplified, and is to give a fast simulation of the situation
% over the year, based on simple input-data.
% This script is a combination of the script SAPV_buthan_01[...] from Stefano Mandelli
% and 'fullYear_script' from Hkon Duus. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTRODUCTION
% Input Notes
% - Input time series (Solar, Load) must be already yearlized (i.e. 365x24=8760
% values needed)
% - Solar: Incident global solar radiation on PV array required
% 
% Plant Notes
% - Inverter sized on peak power
% - Ideal battery:
%         simulation: water tank with charge and discharge efficiency, with
%         limit on power / energy ratio 
%         lifespan: discharge energy given by battery cycle compared with 
%         discharged energy during simulation
% - Ideal PV array: temp effect can be considered
%
% Optimization Notes:
% - Find the optimum plant (minimum NPC) given a max accepted value of LLP


%% Note on the data
% the data from 2009-2013 are data collected from Statnett, and manipulated
% in the script resize. The data called LoadCurve_scaled_2000 is the
% standard data picked from the model of Stefano Mandelli and is not
% manipulated in any way. It is only renamed for running-purposes.

% LoadCurve_scaled_1 is a constant array with the same energy over the year
% as 2000, but with same hourly values every day year around.

% LoadCurve_scaled_2 is a fictive load, with different loads in weeks and
% weekends

% LoadCurve_scaled_3 is the constructed curve from Bhutan

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1
% INITIALIZATION

clear all
close all
beep off
tic                                                           % Start timer for the script

x_llp = linspace(1,20,20);                                    % range of LLP_target (Loss of Load Probability) in [%]
loadCurve_titles = [100];                                     % array with names (i.e. name them by year) of all load curves to be imported. In the case '[100]' there is just one called '100'. 
makePlot = 1;                                                 % set to 1 if plots are desired
columns = length(loadCurve_titles) * 6;                       % since we will be interested in 6 variables at the end
MA_opt_norm_bhut_jun15_20_10 = zeros(length(x_llp), columns); % initialization of the optimal-solution matrix

% Simulation input data
min_PV = 280;               % Min PV power simulated [kW]
max_PV = 300;               % Max PV power simulated [kW]
step_PV = 2;                % PV power simulation step [kW]
min_batt = 50;              % Min Battery capacity simulated [kWh]
max_batt = 800;             % Max Battery capacity simulated [kWh]
step_batt = 5;              % Battery capacity simulation step [kWh]
   
% Computing Number of simulations
n_PV = ((max_PV - min_PV) / step_PV) + 1;                     % N. of simulated PV power sizes (i.e. N. of iteration on PV)
n_batt = ((max_batt - min_batt) / step_batt) + 1;             % N. of simulated Battery capacity (i.e. N. of iteration on Batt)

load_curves_counter = 0;                                      % counter for the number of load curves
    
for year = loadCurve_titles                                   % outer loop going through all the different data sets
    
    clearvars -except x_llp a_x makePlot MA_opt_norm_bhut_jun15_20_10 year loadCurve_titles load_curves_counter min_PV max_PV step_PV n_PV min_batt max_batt step_batt n_batt

    load_curves_counter = load_curves_counter + 1;
    
    % importing 3 data files that describe one year with hourly resolution i.e. 24 x 365 = (8760)-row vectors.                                                
    path_to_dataBase = 'C:\Users\MicrogridProject\Documents\MATLAB\Jeemijn\matlab-microgrid-components\dataBase\';
    irr = importdata([path_to_dataBase, 'solar_data_Phuntsholing_baseline.mat']);                      % Use \ for Windows and / for Mac and Linux
    filename = ([path_to_dataBase, 'LoadCurve_normalized_single_3percent_',num2str(year),'.mat']);   % Average hourly global radiation (beam + diffuse) incident on the PV array [kW/m2]. Due to the simulation step [1h], this is also [kWh/m2]
    Load = importdata(filename);                                                            % Import Load curve 
    T_amb = importdata([path_to_dataBase, 'surface_temp_phuent_2004_hour.mat']);                       % Import ambient temperature data
         
    % Declaration of simulation variables
    EPV = zeros(n_PV, n_batt);              % Energy PV (EPV): yearly energy produced by the PV array [kWh]
    ELPV = zeros(length(irr), n_PV, n_batt);% Energy Loss PV (ELPV): energy produced by the PV array not exploited (i.e. dissipated energy) per time period for each combination of PV and battery [kWh] (Does not include charging losses) 
                                            % N.B. time is the first dimension since later on plot() cannot plot values in 3rd dimension (even if 1st and 2nd dim are scalar) but can plot only values in 1st and 2nd dimension
    LL = zeros(length(irr), n_PV, n_batt);  % Energy not provided to the load: Loss of Load (LL) per time period for each combination of PV and battery [kWh]
    batt_balance = zeros(1,length(irr));    % Powerflow in battery. Positive flow out from battery, negative flow is charging
    num_batt = zeros(n_PV, n_batt);         % number of batteries employed due to lifetime limit
    SoC = zeros(1,size(Load,2));            % to save step-by-step SoC (State of Charge) of the battery
    IC = zeros(n_PV, n_batt);               % Investment Cost (IC) []
    YC = zeros(n_PV, n_batt);               % Operations & Maintenance & replacement; present cost []

    %% System components 
    % System details and input variables are as follows

    % PV panels
    eff_BoS = 0.85;             % Balance Of System: account for such factors as soiling of the panels, wiring losses, shading, snow cover, aging, and so on
    T_ref = 20;                 % Nominal ambient test-temperature of the panels [C] % todo in fullYear script this was 25 and in SAPV it was 20. which one?
    T_nom = 47;                 % Nominal Operating Cell Temperature [C]
    coeff_T_pow = 0.004;        % Derating of panel's power due to temperature [/C]
    irr_nom = 0.8;              % Irradiation at nominal operation [kW / m^2]

    % Battery
    SoC_min = 0.4;              % minimum allowed State Of Charge. This depends on the battery type. The choice of SoC_min influences the lifetime of the batteries.
    SoC_start = 1;              % setting initial State Of Charge
    eff_char = 0.85;            % charge efficiency
    eff_disch = 0.9;            % discharge efficiency
    max_y_repl = 5;             % maximum year before battery replacement
    batt_ratio = 0.5;           % ratio power / energy of the battery (a measure for how fast battery can be (dis)charged)

    % Inverter
    eff_inv = 0.9;              % inverter efficiency

    % Economics
    costPV = 1000;              % PV panel cost [/kW] (source: Uganda data)
    costINV = 500;              % Inverter cost [/kW] (source: MCM_Energy Lab + prof. Silva exercise, POLIMI)
    costOeM_spec = 50;          % Operations & Maintenance cost for the overall plant [/kW*year] (source: MCM_Energy Lab)
    coeff_cost_BoSeI = 0.2;     % Installation (I) and BoS cost as % of cost of PV+battery+Inv [% of Investment cost] (source: Masters, Renewable and Efficient Electric Power Systems,)

    % Battery cost defined as: costBatt_tot = costBatt_coef_a * battery_capacity [kWh] + costBatt_coef_b (source: Uganda data)
    costBatt_coef_a = 140;      % variable cost [per kWh]  %132.78;
    costBatt_coef_b = 0;        % fixed cost
    LT = 20;                    % plant LifeTime [year] 
    r_int = 0.06;               % rate of interest defined as (HOMER) = nominal rate - inflation

    % Simulation input data
    min_PV = 280;               % Min PV power simulated [kW]
    max_PV = 300;               % Max PV power simulated [kW]
    step_PV = 2;                % PV power simulation step [kW]
    min_B = 50;                 % Min Battery capacity simulated [kWh]
    max_B = 800;                % Max Battery capacity simulated [kWh]
    step_B = 5;                 % Battery capacity simulation step [kWh]

	% info not being used:
    % P_mod = 250;                                % Module power in [W]
    % n_mod = ceil(P_syst_des * 1e3 / P_mod);     % Number of modules required for the system of given size
    % a_module = 1.65;                            % Module area in [m^2]
    % V_oc = 34.4;                                % Open circuit voltage [V]
    % P_syst = n_mod * P_mod;                     % Actually installed capacity [W]
    % cycl_B_SoC_min = 2000;                      % number of charge/discharge battery cycle

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 2
    % SYSTEM SIMULATION AND PERFORMANCE INDICATORs COMPUTATION

    %% Plant simulation
    % iterate over all PV power sizes from min_PV to max_PV
    for PV_i = 1 : n_PV                                                 
        PVpower_i = min_PV + (PV_i - 1) * step_PV;                      % iteration on PV power
        T_cell = T_amb + irr .* (T_nom - T_ref) / irr_nom;              % Cell temperature as function of ambient temperature [C]
        eff_cell = 1 - coeff_T_pow .* (T_cell - T_ref);                 % cell efficiency as function of temperature
        P_pv = irr .* PVpower_i .* eff_cell .* eff_BoS;                 % array with Energy from the PV (EPV) for each time step throughout the year. see p.191 of thesis Stefano Mandelli
        
        batt_balance = Load / eff_inv - P_pv;                           % array containing the power balance of the battery for each time step throughout the year (negative value is charging battery) [kWh]
        
        % iterate over all battery capacities from min_batt to max_batt
        for batt_i = 1 : n_batt                                               
            batt_cap_i = min_batt + (batt_i - 1) * step_batt;           % iteration on battery capacity
            EPV(PV_i, batt_i) = sum(P_pv, 2);                           % computing EPV value
            SoC(1) = SoC_start;                                         % setting initial state of charge
            Pow_max = batt_ratio * batt_cap_i;                          % maximum power acceptable by the battery
            Den_rainflow = 0;                                           % counter for number of cycles battery goes through. Needed for CyclesToFailure()

            % iterate through the timesteps of one year
            for t = 1 : size(Load,2)                                    
                if t > 8 
                    if batt_balance(t-1) > 0 && batt_balance(t-2) > 0 && batt_balance(t-3) > 0 && batt_balance(t-4) > 0 && batt_balance(t-5) > 0 && batt_balance(t-6) > 0 && batt_balance(t-7) > 0 && batt_balance(t-8) > 0 && batt_balance(t) < 0   % battery has been charged the previous 8 hours but not this hour.
                       DoD = 1 - SoC(1,t);                              % Depth of Discharge (DoD) is the opposite of State of Charge (SoC)
                       cycles_failure = CyclesToFailure(DoD);
                       Den_rainflow = Den_rainflow + 1/(cycles_failure);
                    end
                end
                
                % charging the battery
                if batt_balance(t) < 0                                   % PV-production is larger than Load. Battery will be charged
                    flow_from_batt = batt_balance(t) * eff_char;         % energy flow that will be stored in the battery i.e. including losses in charging (negative number since charging) [kWh]    % todo this is now negative -> important for plots?
                    if (abs(batt_balance(t))) > Pow_max && SoC(t) < 1    % in-flow exceeds the battery power limit
                        flow_from_batt = Pow_max * eff_char;
                        ELPV(t,PV_i, batt_i) = ELPV(t,PV_i, batt_i) + (abs(batt_balance(t))- Pow_max);
                    end
                    SoC(t+1) = SoC(t) + abs(flow_from_batt) / batt_cap_i;
                    if SoC(t+1) > 1
                        ELPV(t,PV_i, batt_i) = ELPV(t,PV_i, batt_i) + (SoC(t+1) - 1) * batt_cap_i / eff_char;
                        SoC(t+1) = 1;
                    end
                else
                    % discharging the battery
                    flow_from_batt = batt_balance(t) / eff_disch;                                           % total energy flow from the battery i.e. including losses in charging (positive number since discharging) [kWh]    %todo this is now positive -> important for plots?
                    if batt_balance(t) > Pow_max && SoC(t) > SoC_min                                        % checking the battery power limit
                        flow_from_batt = Pow_max / eff_disch;
                        LL(t,PV_i, batt_i) = LL(t,PV_i, batt_i) + (batt_balance(t) - Pow_max) * eff_inv;    % adding the part to LL (Loss of Load) due to exceeding the battery discharging speed
                    end
                    SoC(t+1) = SoC(t) - flow_from_batt / batt_cap_i;
                    if SoC(t+1) < SoC_min
                        LL(t,PV_i, batt_i) = LL(t,PV_i, batt_i) + (SoC_min - SoC(t+1)) * batt_cap_i * eff_disch * eff_inv; % adding the part to LL (Loss of Load) due to not enough energy in battery (using that battery must stay at SoC_min)
                        SoC(t+1) = SoC_min;
                    end
                end
            end

            if batt_i == 2 && PV_i == 2                                 % temporary bad solution
                batt_balance_pos = subplus(batt_balance);               % batt_balance_pos becomes a vector only containing positive values in batt_balance i.e. only interested in when discharging. Negative values = 0
                LL_this = LL(:,PV_i,batt_i);                            % Loss of Load matrix as function of time for these fixed values of PV and battery.
                abs(sum(LL_this) / sum(Load));                          % Finds percentage of Load not served (w.r.t. kWh)
                length(LL_this(find(LL_this<0))) / length(LL_this);     % System Average Interruption Frequency Index (SAIFI), how many hours are without power  (w.r.t. hours)

                if makePlot == 1
                    figure(1)
                    plot(Load,'Color',[72 122 255] / 255)
                    hold on
                    plot(P_pv,'Color',[255 192 33] / 255)
                    hold on
                    plot(batt_balance_pos,'Color',[178 147 68] / 255)
                    hold off
                    xlabel('Time over the year [hour]')
                    ylabel('Energy [kWh]')
                    title('Energy produced and estimated load profile over the year')
                    legend('Load profile','Energy from PV', 'Energy flow in battery')

                    figure(2)
                    plot(Load(1:24),'Color',[72 122 255] / 255)
                    hold on
                    plot(P_pv(1:24),'Color',[255 192 33] / 255)
                    hold on
                    plot(batt_balance_pos(1:24),'Color',[178 147 68] / 255)
                    hold off
                    xlabel('Time over the day [hour]')
                    ylabel('Energy [kWh]')
                    title('Energy produced and estimated load profile of 1 January')
                    legend('Load profile','Energy from PV', 'Energy flow in battery')
                    
                    figure(3)    
                    plot(ELPV(:,PV_i, batt_i) ./ batt_cap_i + 1,'Color',[142 178 68] / 255)
                    hold on
                    plot(LL(:,PV_i, batt_i) ./ batt_cap_i + SoC_min,'Color',[255 91 60] / 255)
                    hold on
                    plot(SoC,'Color',[64 127 255] / 255)
                    hold off
                    xlabel('Time over the year [hour]')
                    ylabel('Power refered to State of Charge of the battery')
                    legend('Overproduction, not utilized', 'Loss of power', 'State of charge')                    
                end
            end

            %% Economic Analysis
            % Investment cost
            costBatt_tot = costBatt_coef_a * batt_cap_i + costBatt_coef_b;              % battery cost
            peak = max(Load);                                                           % peak Load
            costINV_tot = (peak/eff_inv) * costINV;                                     % inverter cost, inverter is designed on the peak power value
            costPV_tot = costPV * PVpower_i;
            costBoSeI = coeff_cost_BoSeI * (costBatt_tot + costINV_tot + costPV_tot);   % cost of Balance of System (BoS) and Installation
            IC(PV_i,batt_i) = costPV_tot + costBatt_tot + costINV_tot + costBoSeI;      % Investment Cost (IC)
            costOeM = costOeM_spec * PVpower_i;                                         % Operations & Maintenance & replacement present cost during plant lifespan
            years_to_go_batt = 1/Den_rainflow;                                          % batteries should be replaced after this number of years
            if years_to_go_batt > max_y_repl
               years_to_go_batt =  max_y_repl;
            end
            num_batt(PV_i,batt_i) = ceil(LT / years_to_go_batt);

            for k = 1 : LT
                if k > years_to_go_batt
                    YC(PV_i,batt_i) = YC(PV_i,batt_i) + costBatt_tot / ((1 + r_int)^years_to_go_batt);                % computing present values of battery
                    years_to_go_batt = years_to_go_batt + years_to_go_batt;
                end
                YC(PV_i,batt_i) = YC(PV_i,batt_i) + costOeM / ((1 + r_int)^k);                                        % computing present values of Operations & Maintenance
            end
            YC(PV_i,batt_i) = YC(PV_i,batt_i) - costBatt_tot * ( (years_to_go_batt - LT) / years_to_go_batt ) / (1 + r_int)^(LT); % salvage due to battery life i.e. estimating how much the batteries are worth after the lifetime of the system
            YC(PV_i,batt_i) = YC(PV_i,batt_i) + costINV_tot / ((1 + r_int)^(LT / 2));                                 % cost of replacing inverter. Assumption: lifetime inverter is half of lifetime system LT
        end
    end

    % Computing Indicators
    NPC = IC + YC;                                                          % Net Present Cost 
    CRF = (r_int * ((1 + r_int)^LT)) / (((1 + r_int)^LT) - 1);              % Capital Recovery Factor
    total_loss_load = squeeze(sum(LL,1));                                   % squeeze() throws away all matrix dimensions with size 1 (in this case the time that has been summed over)
    LLP = total_loss_load / sum(Load, 2);                                   % Loss of Load Probability w.r.t. total load
    LCoE = (NPC * CRF)./(sum(Load, 2) - total_loss_load);                   % Levelized Cost of Energy i.e. cost per kWh (here in ) of building and operating the plant over an assumed life cycle. This is important as we want it to be competitive with the grid LCoE. See eqn. (7.6) in thesis Stefano Mandelli.
        
    save('results.mat')

    % iterate over Loss of Load Probabilities (LLP)
    for a_x = 1 : length(x_llp) 
        LLP_target = x_llp(a_x)/100;                                        % gives LLP in [%] 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 3
        % LOOKING FOR THE OPTIMUM PLANT AS REGARDS THE TARGETED LLP

        LLP_var = 0.3;                                                                              % accepted error band near targeted LLP value
        [posPV, posBatt] = find( (LLP_target - LLP_var) < LLP & LLP < (LLP_target + LLP_var) );     % find possible systems with targeted LLP (within error band). Recall that LLP is a (n_PV x n_batt)-matrix. Example of this syntax: http://se.mathworks.com/help/matlab/ref/find.html#budq84b-1
        NPC_opt = min( diag(NPC(posPV, posBatt)) );                                                 % finds the system within the targeted set that has the minimal NPC
        
        for i = 1 : size(posPV, 1)
            if NPC(posPV(i), posBatt(i)) == NPC_opt
                PV_opt = posPV(i);
                Batt_opt = posBatt(i);
            end
        end

        kW_opt = (PV_opt - 1) * step_PV + min_PV;
        kWh_opt = (Batt_opt - 1) * step_batt + min_batt;
        LLP_opt = LLP(PV_opt, Batt_opt);
        LCoE_opt = LCoE(PV_opt, Batt_opt);
        IC_opt = IC(PV_opt, Batt_opt);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 4
        % PLOTTING

        if makePlot == 1
%         if false
%             figure(4);
%             mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, NPC);
%             title('Net Present Cost');
%             set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
%             xlabel('Battery Bank size [kWh]');
%             ylabel('PV array size [kW]');

            figure(5);
            mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, LLP);
            title('Loss of Load Probability');
            set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
            xlabel('Battery Bank size [kWh]');
            ylabel('PV array size [kW]');

            figure(6);
            mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, LCoE);
            title('Levelized Cost of Energy');
            set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
            xlabel('Battery Bank size [kWh]');
            ylabel('PV array size [kW]');

            figure(7);
            mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, num_batt);
            title('Num. of battery employed due to lifetime limit');
            set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
            xlabel('Battery Bank size [kWh]');
            ylabel('PV array size [kW]');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 5
        % Making optimal-solution matrix

        if isempty(NPC_opt) == 1
            NPC_opt = NaN;
        end

        opt_sol = [LLP_opt NPC_opt kW_opt kWh_opt LCoE_opt IC_opt];

        MA_opt_norm_bhut_jun15_20_10(a_x, ((6 * load_curves_counter - 5) : 6 * load_curves_counter)) = opt_sol;     % the : operator sets the range of y-coordinates that the array opt_sol will take in the matrix MA_opt_norm_bhut_jun15_20_10
    end
end

toc % End timer
save('MA_opt_norm_bhut_jun15_20_10.mat','MA_opt_norm_bhut_jun15_20_10')
