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
% GET INPUT DATA

clear all
close all
beep off
tic
% LLP_target = input('Required maximum accepted value of LLP [%]: ') / 100; % enter user defined targeted value of LLP (Loss of Load Probability)
x_llp=linspace(1,20,20);                                    % range of LLP_target (Loss of Load Probability) in [%]
loadCurve=[100];                                            % number of data sets
makePlot=0;                                                 % set to 1 if plots are desired
columns=length(loadCurve)*6;
MA_opt_norm_bhut_jun15_20_10=zeros(length(x_llp),columns);  % initialization of the optimal-solution matrix
counter=0;
    
% tags={'LLP_opt', 'NPC', 'PV opt', 'Bat_opt', 'LCoE'};
for year=loadCurve                                          % outer loop going through all the different data sets
    
    clearvars -except x_llp a_x makePlot MA_opt_norm_bhut_jun15_20_10 year loadCurve counter

    counter=counter+1;
            
    % load(filename, variable_name)    
    load('solar_data_Phuntsholing_baseline.mat', 'solar_data_Phuntsholing_baseline')    % Average hourly global radiation (beam + diffuse) incident on the PV array [kW/m2]. Due to the simulation step [1h], this is also [kWh/m2]
    % load('Tcell_Soroti_h_year.mat', 'Tcell_Soroti_h_year')                            % Cell Temperature [°C]
    irr = solar_data_Phuntsholing_baseline;
    filename=(['LoadCurve_normalized_single_3percent_',num2str(year),'.mat'])
    data=importdata(filename);                                                          % Import load data 
    T_amb=importdata('surface_temp_phuent_2004_hour.mat');                              % Import ambient temperature data
    T_nom=47;                                                                           % Nominal Operating Cell Temperature [°C]
    Load = data;                                                                        % Load Curve
    T_cell = T_amb+irr.*(T_nom-20)/0.8;                                                 % Cell temperature as function of ambient temperature
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 2
    % SYSTEM COMPONENTS INPUT DATA

    % PV array
    eta_BoS = 0.85;             % Balance Of System: account for such factors as soiling of the panels, wiring losses, shading, snow cover, aging, and so on
    T_ref = 25;                 % Nominal ambient test-temperature of the panels [C]
    coeff_T_pow = 0.004;        % Derating of panel's power due to temperature

    % Battery
    SoC_min = 0.4;              % minimum allowed State Of Charge
    SoC_start = 1;              % setting initial State Of Charge
    eta_char = 0.85;            % charge efficiency
    eta_disch = 0.9;            % discharge efficiency
    cycl_B_SoC_min = 2000;      % number of charge/discharge battery cycle
    max_y_repl = 5;             % maximum year before battery replacement
    batt_ratio = 0.5;           % ratio power / energy of the battery

    % Inverter
    eta_inv = 0.9;              % inverter efficiency

    % Economics
    costPV = 1000;              % PV panel cost [€/kW] (source: Uganda data)
    costINV = 500;              % Inverter cost [€/kW] (source: MCM_Energy Lab + prof. Silva exercise, POLIMI)
    costOeM_spec = 50;          % O&M cost for the overall plant [€/kW*year] (source: MCM_Energy Lab)
    coeff_cost_BoSeI = 0.2;     % Installation and BoS cost as % of cost of PV+B+Inv [% of Investment cost] (source: Masters, “Renewable and Efficient Electric Power Systems,”)

    % Battery cost defined as: costB = costB_coef_a * battery_capacity [kWh] + costB_coef_b (source: Uganda data)
    costB_coef_a = 140;         %132.78;
    costB_coef_b = 0;
    VU = 20;                    % plant lifetime [year] 
    r_int = 0.06;               % rate of interest defined as (HOMER) = nominal rate - inflation

    % Simulation input data
    min_PV = 285;               % Min PV power simulated [kW]
    max_PV = 300;               % Max PV power simulated [kW]
    step_PV = 5;                % PV power simulation step [kW]
    min_B = 50;                 % Min Battery capacity simulated [kWh]
    max_B = 800;                % Max Battery capacity simulated [kWh]
    step_B = 2;                 % Battery capacity simulation step [kWh]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 3
    % INPUT DATA ELABORATION

    % Computing Number (N) of simulations
    n_PV = ((max_PV - min_PV) / step_PV) + 1;   % N. of simulated PV power sizes (i.e. N. of iteration on PV)
    n_B = ((max_B - min_B) / step_B) + 1;       % N. of simulated Battery capacity (i.e. N. of iteration on Batt)
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PART 4
    % SYSTEM SIMULATION AND PERFORMANCE INDICATORs COMPUTATION

    % Declaration of simulation variables
    EPV = zeros(n_PV, n_B);         % Energy PV (EPV): yearly energy produced by the PV array [kWh]
    ELPV = zeros(n_PV, n_B);        % Energy Loss PV (ELPV): yearly energy produced by the PV array not exploited (i.e. dissipated energy) [kWh]
    LL = zeros(n_PV, n_B);          % Energy not provided to the load: Loss of Load (LL) [kWh]
    IC = zeros(n_PV, n_B);          % Investment Cost (IC)
    YC = zeros(n_PV, n_B);          % O&M & replacement present cost
    num_B = zeros(n_PV, n_B);

    % single battery simulation variable
    SoC = zeros(1,size(Load,2));    % to save step-by-step SoC (State of Charge)

    n = 0;
    %% Plant simulation
    % iterate over all PV power sizes from min_PV to max_PV
    for PV_i = 1 : n_PV                                                 
        n = n + 1;
        PVpower_i = min_PV + (PV_i - 1) * step_PV;                      % iteration on PV power
        eta_cell = 1 - coeff_T_pow .* (T_cell - T_ref);                 % cell efficiency as function of temperature
        P_pv = irr .* PVpower_i .* eta_cell .* eta_BoS;                 % array with Energy from the PV (EPV) for each time step throughout the year
        batt_balance = Load / eta_inv - P_pv;                           % array containing the power balance of the battery for each time step throughout the year (negative value is charging battery) [kWh]
        
        % iterate over all battery capacities from min_B to max_B
        for B_i = 1 : n_B                                               
            Bcap_i = min_B + (B_i - 1) * step_B;                        % iteration on battery capacity
            EPV(PV_i, B_i) = sum(P_pv, 2);                              % computing EPV value
            SoC(1, 1) = SoC_start;                                      % setting initial state of charge
            Pow_max = batt_ratio * Bcap_i;                              % maximum power acceptable by the battery
            Den_rainflow = 0;
            
            % iterate through the timesteps in one year
            for k = 1 : size(Load,2)                                    
                if k > 8
                    if batt_balance(1, k-1) > 0 && batt_balance(1, k-2) > 0 && batt_balance(1, k-3) > 0 && batt_balance(1, k-4) > 0 && batt_balance(1, k-5) > 0 && batt_balance(1, k-6) > 0 && batt_balance(1, k-7) > 0 && batt_balance(1, k-8) > 0 && batt_balance(1, k) < 0 
                       DOD = 1 - SoC(k);
                       cycles_failure = CyclesToFailure(DOD);
                       Den_rainflow = Den_rainflow + 1/(cycles_failure);
                    end
                end
                if batt_balance(1, k) < 0                                   % charging battery
                    EB_flow = batt_balance(1, k) * eta_char;                % [kWh] to the battery % todo this is now negative 
                    PB_flow = (EB_flow / eta_char);
                    if PB_flow > Pow_max && SoC(1, k) < 1                   % checking the battery power limit
                        EB_flow = Pow_max * eta_char;
                        ELPV(PV_i, B_i) = ELPV(PV_i, B_i) + (PB_flow - Pow_max);
                    end
                    SoC(1, k+1) = SoC(1, k) + EB_flow / Bcap_i;
                    if SoC(1, k+1) > 1
                        ELPV(PV_i, B_i) = ELPV(PV_i, B_i) + (SoC(1, k+1) - 1) * Bcap_i / eta_char;
                        SoC(1, k+1) = 1;
                    end
                else
                    % discharging battery
                    EB_flow = batt_balance(1, k) / eta_disch;                           % [kWh] from the battery %todo this is now positive
                    PB_flow = (abs(EB_flow) * eta_disch);
                    if PB_flow > Pow_max && SoC(1, k) > SoC_min                         % checking the battery power limit
                        EB_flow = Pow_max / eta_disch;
                        LL(PV_i, B_i) = LL(PV_i, B_i) + (PB_flow - Pow_max) * eta_inv;  % computing numerator of LLP indicator
                    end
                    SoC(1, k+1) = SoC(1, k) + EB_flow / Bcap_i;
                    if SoC(1, k+1) < SoC_min
                        LL(PV_i, B_i) = LL(PV_i, B_i) + (SoC_min - SoC(1, k+1)) * Bcap_i * eta_disch * eta_inv; % computing numerator of LLP indicator
                        SoC(1, k+1) = SoC_min;
                    end
                end
            end

            % Economic Analysis

            % Investment cost
            costB = costB_coef_a * Bcap_i + costB_coef_b;                           % battery cost
            Pmax = max(Load);                                                       % Peak Load
            costI = (Pmax/eta_inv) * costINV;                                       % Inverter cost, inverter is designed on the peak power value
            costBoSeI = coeff_cost_BoSeI * (costB + costI + costPV * PVpower_i);    % cost of Balance Of Systems and Installation
            IC(PV_i,B_i) = costPV * PVpower_i + costB + costI + costBoSeI;          % Investment Cost
            costOeM = costOeM_spec * PVpower_i;                                     % O&M & replacement present cost during plant lifespan
            y_rep_B = 1/Den_rainflow;

            if y_rep_B > max_y_repl
               y_rep_B =  max_y_repl;
            end

            num_B(PV_i,B_i) = ceil(VU / y_rep_B);
            cfr = y_rep_B;

            for k = 1 : VU
                if k > cfr
                    YC(PV_i,B_i) = YC(PV_i,B_i) + costB / ((1 + r_int)^cfr);        % computing present values of battery
                    cfr = cfr + y_rep_B;
                end
                YC(PV_i,B_i) = YC(PV_i,B_i) + costOeM / ((1 + r_int)^k);            % computing present values of O&M
            end
            YC(PV_i,B_i) = YC(PV_i,B_i) - costB * ( (cfr - VU) / y_rep_B ) / (1 + r_int)^(VU);  % salvage due to battery life
            YC(PV_i,B_i) = YC(PV_i,B_i) + costI / ((1 + r_int)^(VU / 2));                       % considering present cost of inverter
        end
    end

    % Computing Indicators
    NPC = IC + YC;                                                          % Net Present Cost 
    CRF = (r_int * ((1 + r_int)^VU)) / (((1 + r_int)^VU) - 1);              % LCoE coefficient
    LLP = LL / sum(Load, 2);                                                % Loss of Load Probability
    LCoE = (NPC * CRF)./((sum(Load, 2)).*(1 - LLP));                        % Levelized Cost of Energy
        
    save('results.mat')

    for a_x=1:length(x_llp) 
        LLP_target=x_llp(a_x)/100;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 5
        % LOOKING FOR THE OPTIMUM PLANT AS REGARDS THE TARGETED LLP

        LLP_var = 0.005;                                                                        % accepted error band near targeted LLP value
        [posPV, posB] = find( (LLP_target - LLP_var) < LLP & LLP < (LLP_target + LLP_var) );    % find systems with targeted LLP (within error band)
        NPC_opt = min( diag(NPC(posPV, posB)) );                                                % find minimum NPC for the target LLP;
        
        for i = 1 : size(posPV, 1)
            if NPC(posPV(i), posB(i)) == NPC_opt
                PV_opt = posPV(i);
                B_opt = posB(i);
            end
        end

        kW_opt = (PV_opt-1) * step_PV + min_PV;
        kWh_opt = (B_opt-1) * step_B + min_B;
        LLP_opt = LLP(PV_opt, B_opt)
        LCoE_opt = LCoE(PV_opt, B_opt);
        IC_opt=IC(PV_opt, B_opt);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 6
        % PLOTTING

        if makePlot==1
            figure(1);
            mesh(min_B : step_B : max_B , min_PV : step_PV : max_PV , NPC);
            title('Net Present Cost');
            set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
            xlabel('Battery Bank size [kWh]');
            ylabel('PV array size [kW]');

            figure(2);
            mesh(min_B : step_B : max_B , min_PV : step_PV : max_PV , LLP);
            title('Loss of Load Probability');
            set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
            xlabel('Battery Bank size [kWh]');
            ylabel('PV array size [kW]');

            figure(3);
            mesh(min_B : step_B : max_B , min_PV : step_PV : max_PV , LCoE);
            title('Levelized Cost of Eenergy');
            set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
            xlabel('Battery Bank size [kWh]');
            ylabel('PV array size [kW]');

            figure(4);
            mesh(min_B : step_B : max_B , min_PV : step_PV : max_PV , num_B);
            title('Num. of battery employed due to lifetime limit');
            set(gca,'FontSize',12,'FontName','Times New Roman','fontWeight','bold')
            xlabel('Battery Bank size [kWh]');
            ylabel('PV array size [kW]');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% PART 7
        % Make optimal-solution matrix

        if isempty(NPC_opt)==1
            NPC_opt=NaN;
        end

        opt_sol=[LLP_opt NPC_opt kW_opt kWh_opt LCoE_opt IC_opt];

        MA_opt_norm_bhut_jun15_20_10(a_x,((6*counter-5):6*counter))=opt_sol;
    end
end

save('MA_opt_norm_bhut_jun15_20_10.mat','MA_opt_norm_bhut_jun15_20_10')

toc
