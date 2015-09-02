%% This script is a mathemathical model of the simulink-version of the microgrid.
% It is very simplified, and is to give a fast simulation of the situation
% over the year, based on simple input-data.
% This script is a combination of the script SAPV_buthan_01[...] from Stefano Mandelli
% and 'fullYear_script' from Hakon Duus. 

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
% manipulated in any way. It is only renamed for running purposes.

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
tic                                         % Start timer for the script

% choose a mode to run the script: 
% 1 for fixed LLP, 2 for fixed cost, 3 for computing for all LLP and all cost
mode = 3;

if mode == 1
    % specify if mode = 1 (fixed LLP):
    LLP_fixed = 50;                         % aimed LLP in [%]. The program will find the lowest budget for this LLP.
    
    disp(['The program runs in mode 1 (fixed LLP): for a fixed LLP of ',num2str(LLP_fixed), '% it looks for the lowest cost.']);
    range_LLP = LLP_fixed;                  % in this mode the range only consists of 1 fixed value
elseif mode == 2
    % specify if mode = 2 (fixed cost):
    budget_fixed = 800000;                  % aimed budget in [EUR]. The program will find the lowest LLP for this budget.
    
    disp(['The program runs in mode 2 (fixed cost): for a fixed cost of EUR ', num2str(budget_fixed), ' it looks for the lowest LLP.']);
    range_LLP = [];                         % in this mode we do not iterate over LLP (only over cost) so we skip the for loop over range_LLP.
elseif mode == 3
    % specify if mode = 3 (all LLP, all cost):
    range_LLP = linspace(0,100,101);        % range of LLP_target (Loss of Load Probability) in [%]
    
    disp(['The program runs in mode 3 (all LLP, all cost). It runs for all LLP between ',num2str(range_LLP(1)),'% and ',num2str(range_LLP(end)),'% in ',num2str(length(range_LLP)-1),' steps.']);
else
    error('ERROR: please specify a mode to run the program in. Choose mode=1 for fixed LLP, mode=2 for fixed cost or mode=3 for computing for all LLP and all cost.');
end

makePlot = 0;                                                       % set to 1 if plots are desired
loadCurve_titles = [100];                                           % array with names (i.e. name them by year) of all load curves to be imported. In the case '[100]' there is just one called '100'. 
columns = length(loadCurve_titles) * 6;                             % since we will be interested in 6 variables at the end
MA_opt_norm_bhut = zeros(length(range_LLP), columns);               % initialization of the optimal-solution matrix

% Simulation input data
min_PV = 0;                 % Min PV power simulated [kW]
max_PV = 300;               % Max PV power simulated [kW]
step_PV = 5;                % PV power simulation step [kW]
min_batt = 0;               % Min Battery capacity simulated [kWh]
max_batt = 800;             % Max Battery capacity simulated [kWh]
step_batt = 10;             % Battery capacity simulation step [kWh]

% Computing Number of simulations
n_PV = ((max_PV - min_PV) / step_PV) + 1;                     % N. of simulated PV power sizes (i.e. N. of iteration on PV)
n_batt = ((max_batt - min_batt) / step_batt) + 1;             % N. of simulated Battery capacity (i.e. N. of iteration on Batt)

if mod(max_PV - min_PV,step_PV) ~= 0 
    error('ERROR: The range of PV sizes is not divisible by the step size step_PV.')
end

if mod(max_batt - min_batt,step_batt) ~= 0 
    error('ERROR: The range of battery sizes is not divisible by the step size step_batt.')
end

load_curves_counter = 0;                                      % counter for the number of load curves
    
for year = loadCurve_titles                                   % outer loop going through all the different data sets
    
    clearvars -except x_llp a_x makePlot MA_opt_norm_bhut year loadCurve_titles load_curves_counter min_PV max_PV step_PV n_PV min_batt max_batt step_batt n_batt range_LLP LLP_fixed mode budget_fixed

    load_curves_counter = load_curves_counter + 1;
        
    % importing 3 data files that describe one year with hourly resolution i.e. 24 x 365 = (8760)-row vectors.                                                
    path_to_dataBase = '/Users/jeemijn/Desktop/NTNU_microgrids/matlab/matlab-microgrid-components/dataBase/';
    irr = importdata([path_to_dataBase, 'solar_data_Phuntsholing_baseline.mat']);                       % Use \ for Windows and / for Mac and Linux
    filename = ([path_to_dataBase, 'LoadCurve_normalized_single_3percent_',num2str(year),'.mat']);      % Average hourly global radiation (beam + diffuse) incident on the PV array [kW/m2]. Due to the simulation step [1h], this is also [kWh/m2]
    Load = importdata(filename);                                                                        % Import Load curve 
    T_amb = importdata([path_to_dataBase, 'surface_temp_phuent_2004_hour.mat']);                        % Import ambient temperature data
    
    % Declaration of simulation variables
    EPV = zeros(n_PV, n_batt);              % Energy PV (EPV): yearly energy produced by the PV array [kWh]
    ELPV = zeros(length(irr), n_PV, n_batt);% Energy Loss PV (ELPV): energy produced by the PV array not exploited (i.e. dissipated energy) per time period for each combination of PV and battery [kWh] (Does not include charging losses) 
                                            % N.B. time is the first dimension since later on plot() cannot plot values in 3rd dimension (even if 1st and 2nd dim are scalar) but can plot only values in 1st and 2nd dimension
    LL = zeros(length(irr), n_PV, n_batt);  % Energy not provided to the load: Loss of Load (LL) per time period for each combination of PV and battery [kWh]
    batt_balance = zeros(1,length(irr));    % Powerflow in battery. Positive flow out from battery, negative flow is charging
    num_batt = zeros(n_PV, n_batt);         % number of batteries employed due to lifetime limit
    SoC = zeros(1,size(Load,2));            % to save step-by-step SoC (State of Charge) of the battery
    IC = zeros(n_PV, n_batt);               % Investment Cost (IC) [EUR]
    YC = zeros(n_PV, n_batt);               % Operations & Maintenance & replacement; present cost [EUR]
    accept(1:length(range_LLP)) = -1;       % To check whether the found optimal solution is the best one. Declaring a value of -1 everywhere since we are filling it with 0 and 1.

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
                    if batt_cap_i == 0
                        SoC(t+1) = SoC(1);                              % to undo division by zero
                    end
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
                    if batt_cap_i == 0
                        SoC(t+1) = SoC(1);                              % to undo division by zero
                    end
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
                    plot(P_pv,'Color',[255 192 33] / 255)
                    hold on
                    plot(batt_balance_pos,'Color',[178 147 68] / 255)
                    hold on
                    plot(Load,'Color',[72 122 255] / 255)
                    hold off
                    set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
                    xlabel('Time over the year [hour]')
                    ylabel('Energy [kWh]')
                    title('Energy produced and estimated load profile over the year (2nd steps PV and Batt)')
                    legend('Energy from PV', 'Energy flow from battery','Load profile')

                    % integrate figure(1) to find rough LLP estimate
                    free = min(Load, P_pv);                                         % energy for free, i.e. directly from PV without battery intervenience, is the area under this graph. 
                                                                                    % N.B. Assumption: both Load and P_pv are positive functions
                    time = 1:length(irr);
                    free_area = trapz(time, free);                                  % discrete integration using trapeziums. Might not be the best solution for non-linear data. See http://se.mathworks.com/help/matlab/math/integration-of-numeric-data.html
                    area_to_batt = trapz(time, P_pv) - free_area;                   % by definition this should be positive
                    area_load_needed_from_batt = trapz(time, Load) - free_area;     % by definition this should be positive

                    unmet_load = area_load_needed_from_batt - area_to_batt;
                    unmet_load_perc = unmet_load / trapz(time, Load) * 100;          % equal to Loss of Load Probability. But rough estimate since only comparing totals of load and P_pv! And SoC at end of the day influences next day. (Negative means overproduction)

                    % plot functions for an average day in figure(2)
                    Load_av = DailyAverage(Load);
                    P_pv_av = DailyAverage(P_pv);
                    batt_balance_pos_av = DailyAverage(batt_balance_pos);           % This average is misleading/not so useful since it is influenced by state of charge of previous days.
                    
                    figure(2)
                    plot(Load_av,'Color',[72 122 255] / 255)
                    hold on
                    plot(P_pv_av,'Color',[255 192 33] / 255)
                    hold on
                    plot(batt_balance_pos_av,'Color',[178 147 68] / 255)
                    hold off
                    set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
                    xlabel('Time over the day [hour]')
                    ylabel('Energy [kWh]')
                    title('Energy produced and estimated load profile of an average day (2nd steps PV and Batt)')
                    legend('Load profile','Energy from PV', 'Energy flow in battery')
                                        
                    figure(3)    
                    plot(ELPV(:,PV_i, batt_i) ./ batt_cap_i + 1,'Color',[142 178 68] / 255)         % not nice for batt_cap_i == 0
                    hold on
                    plot(LL(:,PV_i, batt_i) ./ batt_cap_i + SoC_min,'Color',[255 91 60] / 255)
                    hold on
                    plot(SoC,'Color',[64 127 255] / 255)
                    hold off
                    set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
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
    for a_x = 1 : length(range_LLP) 
        LLP_target = range_LLP(a_x)/100;                                    % gives LLP in [%] 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 3
        % LOOKING FOR THE OPTIMUM PLANT AS REGARDS THE TARGETED LLP

        clear posPV posBatt NPC_opt kW_opt kWh_opt LLP_opt LCoE_opt IC_opt this_LLP this_LLP_costs
        
        LLP_var = 0.005;                                                                            % accepted error band near targeted LLP value. Note that values are not in % but in [0,1]. So LLP_var = 0.005 means that we look for an LLP of 14% within 0.135 to 0.145 (i.e. 13.5% to 14.5%)
        [posPV, posBatt] = find( (LLP_target - LLP_var) < LLP & LLP < (LLP_target + LLP_var) );     % find possible systems with targeted LLP (within error band). Recall that LLP is a (n_PV x n_batt)-matrix. Example of this syntax: http://se.mathworks.com/help/matlab/ref/find.html#budq84b-1
        
        if isempty(posPV)
            continue;                           % exit this loop if no values found for this LLP_target
        end
        
        NPC_opt = min( diag(NPC(posPV, posBatt)) );                                                 % finds the system within the targeted set that has the minimal NPC
        
        for i = 1:size(posPV, 1)
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
        
        

        % check that the optimal values are within the search range of PV/battery:
        
        LLP_flipped = flipud(LLP);                  % flipping LLP matrix up-down s.t. find() will search through the matrix in the order along the LLP isopleths
        NPC_flipped = flipud(NPC);                  % flipping NPC matrix up-down s.t. matrix coordinates correspond with LLP_flipped
        this_LLP = find((LLP_target - LLP_var) < LLP_flipped & LLP_flipped < (LLP_target + LLP_var));       % gives index nrs of LLP matrix that correspond to this LLP isopleth
        this_LLP_costs = NPC_flipped(this_LLP);                                                                     % gives cost (NPC) along this LLP isopleth
        
        [ymax,xmax,ymin,xmin] = extrema_no_boundaries(this_LLP_costs);
        accept(a_x) = ~isempty(xmin);                    % accept if there is at least one minimum in the cost curve along the LLP isopleth
        disp(['max: ', num2str(length(xmax)), '  min: ',num2str(length(xmin)), '  accept: ',num2str(accept(a_x)),'  LLP: ',num2str(LLP_target*100)])
        
        % Plot the cost along this LLP isopleth with local minima and maxima
        % The optimal system is within the search range if at least 1 minimum
        
        figure(9)
        x = 1:length(this_LLP_costs);
        plot(x,this_LLP_costs,'-o')        
        hold on
        plot(x(xmax),ymax,'r*',x(xmin),ymin,'g*')        
        hold off
        set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
        title(['Extrema of costs along the LLP = ',num2str(LLP_target*100),'% isopleth']);
        xlabel(['Nr of system along the LLP = ',num2str(LLP_target*100),'% isopleth']);
        ylabel('Costs (NPC) [EUR]');
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 4
        % Making optimal-solution matrix

        if isempty(NPC_opt) == 1
            NPC_opt = NaN;
        end

        opt_sol = [LLP_opt NPC_opt kW_opt kWh_opt LCoE_opt IC_opt];

        MA_opt_norm_bhut(a_x, ((6 * load_curves_counter - 5) : 6 * load_curves_counter)) = opt_sol;     % the : operator sets the range of y-coordinates that the array opt_sol will take in the matrix MA_opt_norm_bhut        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 5
    % PLOTTING (for last value of LLP in a_x)

    if makePlot == 1
%         if false
        figure(4);
        mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, NPC);
        set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
        title('Net Present Cost');
        xlabel('Battery Bank size [kWh]');
        ylabel('PV array size [kW]');

        figure(5);
        mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, LLP);
        set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
        title('Loss of Load Probability');
        xlabel('Battery Bank size [kWh]');
        ylabel('PV array size [kW]');

        figure(6);
        mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, LCoE);
        set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
        title('Levelized Cost of Energy');
        xlabel('Battery Bank size [kWh]');
        ylabel('PV array size [kW]');

        figure(7);
        mesh(min_batt : step_batt : max_batt, min_PV : step_PV : max_PV, num_batt);
        set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
        title('Num. of battery employed due to lifetime limit');
        xlabel('Battery Bank size [kWh]');
        ylabel('PV array size [kW]');
    end
end

%% PART 6 
% make a plot of resulting systems PV vs Batt size, where the color of
% dots gives LLP value and isopleths of budget are added in black

% plotting dots for each system (x,y)=(PV, Batt) with color corresponding to
% LLP
nr_of_systems = n_PV * n_batt;
x_values = zeros(1,nr_of_systems);
y_values = zeros(1,nr_of_systems);
colors = zeros(1,nr_of_systems);

dots_counter = 0;
for i = 1:n_PV
    for j = 1:n_batt
        dots_counter = dots_counter + 1;
        this_PV = min_PV + (i - 1) * step_PV;
        this_batt = min_batt + (j - 1) * step_batt;
        
        x_values(dots_counter) = this_batt;                     % we want to plot batt on x-axis and PV on y-axis (in NPC and LLP matrices it is the other way around)
        y_values(dots_counter) = this_PV;                   
%         colors(dots_counter) = LLP(i,j)*100;                    % choose colour of the dot according to value of Loss of Load Probability in [%]. 
        colors(dots_counter) = roundn(LLP(i,j)*100,1);          % choose colour of the dot according to value of Loss of Load Probability in [%]. Rounded to steps of 10% s.t. colour differences in the plot can be seen better.
    end
end

figure(8);
if dots_counter > 0
    scatter(x_values, y_values, [], colors, 'filled')            
    hold on
end

% rounding costs in order to compare up to x digits with certain values of isopleths
NPC_rounded_flipped = roundn(flipud(NPC),4);            % round the cost to 10^4 EUR
min = roundn(min(NPC(:)),4);                            % lowest cost that occurs rounded to 10^4 EUR
max = roundn(max(NPC(:)),4);                            % highest cost that occurs rounded to 10^4 EUR
budget_range = linspace(min, max, 10);

% plotting isopleths of equal cost (NPC) in black
for cost = budget_range
    [cost_x, cost_y_flipped] = find(NPC_rounded_flipped == cost);       % find all systems (x,y) = (PV_i, batt_i) for this cost
    cost_x = min_PV + (cost_x - 1) * step_PV;           % convert to correct values of PV i.e. in [kW] instead of their order 1, 2, 3, ...
    cost_y_flipped = max_batt - (cost_y_flipped - 1) * step_batt;
    plot(cost_y_flipped, cost_x,'k-o','linewidth',1.1)
end

hold off
bar = colorbar;
set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
ylabel(bar,'Loss of Load Probability [%]')
xlabel('Battery bank size [kWh]')
ylabel('PV array size [kW]')
title('The LLP of each (batt, PV) system in colour. Isopleths of equal cost (NPC) are in black.')

%%
toc % End timer
save('MA_opt_norm_bhut.mat','MA_opt_norm_bhut')

if mode == 1
    disp(['For the fixed LLP of ',num2str(LLP_fixed),'% the optimal system costs EUR ',num2str(round(NPC_opt)),' and is given in MA_opt_norm_bhut.mat.'])
    disp('The columns mean LLP_opt NPC_opt PV_opt[kW] batt_opt[kW] LCoE_opt IC_opt, respectively.')
    if accept(1) == 0
        warning('The true optimal system may lie outside the range of searched PV/batt sizes (accept=0).')
    end
elseif mode ==3
    disp('For each LLP in the range the optimal system (lowest cost) is given in MA_opt_norm_bhut.mat.')
    disp('Each row gives the optimal system for one value of LLP, in the order of LLP_range.')
    disp('The columns mean LLP_opt NPC_opt PV_opt[kW] batt_opt[kW] LCoE_opt IC_opt, respectively.')
    
    problematic_LLPs = range_LLP(find(accept == 0));
    if ~isempty(problematic_LLPs)
       warning(['For the ',num2str(length(problematic_LLPs)),' values of LLP given in the vector ''problematic_LLPs'' (in %), the true optimal system may lie outside the range of searched PV/batt sizes (accept=0).']) 
    end
end
