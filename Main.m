%% LMP & DLMP main code
%  12-18-2024

clc
clear all

%% 定义选项
mpopt = mpoption( 'out.all', 0);%, 'opf.ac.solver', 'ipopt'
% mpopt = mpoption;
define_constants;

%% 导入系统数据
case_name = {'T-1';'T-2';'T-3';'T-4';'T-5';'T-6';'T-7';'T-8';'T-9';...
    'T-10';'T-11';'T-12';'T-13';'T-14';'T-15'};

mpc = {loadcase('trans_opf_case9_ieee'); loadcase('trans_opf_case14_ieee'); loadcase('trans_opf_case57_ieee'); ...% 39 57 *0.5; 118 *1.2
    loadcase('trans_opf_case_rts96'); loadcase('trans_opf_case118_ieee');loadcase('case162');...
    loadcase('trans_opf_case_rts79');loadcase('trans_opf_case793_goc');loadcase('trans_opf_case300_ieee');...
    loadcase('grid_IEEE123'); loadcase('dist_opf_case136_ieee'); loadcase('dist_opf_case51_ieee');...
    loadcase('dist_opf_case69_ieee'); loadcase('dist_opf_case74_ieee');loadcase('dist_opf_case141_ieee')};%case4 case9 case11 case15

%% 不同模型
error = {};
for c = 1
    %% MATPOWER ACOPF
    Res_BaseAC{c} = runopf(mpc{c}, mpopt);
    %     Pbus = real(makeSbus(Res_BaseAC{c}.baseMVA, Res_BaseAC{c}.bus, Res_BaseAC{c}.gen));
    %     Qbus = imag(makeSbus(Res_BaseAC{c}.baseMVA, Res_BaseAC{c}.bus, Res_BaseAC{c}.gen));
    %       Vm0 = Res_BaseAC{c}.bus(:,8);
    %       Va0 = pf_results.bus(:,9)*pi/180;
    Res_BaseAC{c}.f
    Ploss_base = real(get_losses(Res_BaseAC{c}))/ mpc{c}.baseMVA;
    Qloss_base = imag(get_losses(Res_BaseAC{c}))/ mpc{c}.baseMVA;
    
    % Ysc = 1 ./ (Res_BaseAC{c}.branch(:, BR_R) - 1j * Res_BaseAC{c}.branch(:, BR_X));
    % f = Res_BaseAC{c}.branch(:, F_BUS);
    % t = Res_BaseAC{c}.branch(:, T_BUS);
    % Ploss = Res_BaseAC{c}.baseMVA*real(Ysc).* (Res_BaseAC{c}.bus(f,8).^2+Res_BaseAC{c}.bus(t,8).^2-2*Res_BaseAC{c}.bus(f,8).*Res_BaseAC{c}.bus(t,8).*cos(Res_BaseAC{c}.bus(f,9)* pi/180-Res_BaseAC{c}.bus(t,9)* pi/180));
    % Qloss = Res_BaseAC{c}.baseMVA*imag(Ysc).* (Res_BaseAC{c}.bus(f,8).^2+Res_BaseAC{c}.bus(t,8).^2-2*Res_BaseAC{c}.bus(f,8).*Res_BaseAC{c}.bus(t,8).*cos(Res_BaseAC{c}.bus(f,9)* pi/180-Res_BaseAC{c}.bus(t,9)* pi/180));
    %
    %     p_branch = (Res_BaseAC{c}.branch(:,14) - Res_BaseAC{c}.branch(:,16))/2/ Res_BaseAC{c}.baseMVA;
    %     r = Res_BaseAC{c}.branch(:,3);
    %     Ploss0 = p_branch.^2 .* r *Res_BaseAC{c}.baseMVA;
    %     %Ploss0 = sum(real(get_losses(pf_results)))/ mpc{c}.baseMVA;
    %
    % for i = 1:size(Res_BaseAC{c}.bus(:,1),1)
    %     Mi = find(f == i | t == i);
    %     Pequ(i,1) = 0.5 * sum(Ploss(Mi));
    %     Qequ(i,1) = 0.5 * sum(Qloss(Mi));
    % end
    % Cf = full(sparse((1:size(Res_BaseAC{c}.branch, 1))', f, 1, size(Res_BaseAC{c}.branch, 1), size(Res_BaseAC{c}.bus(:,1),1)));
    % Ct = full(sparse((1:size(Res_BaseAC{c}.branch, 1))', t, 1, size(Res_BaseAC{c}.branch, 1), size(Res_BaseAC{c}.bus(:,1),1)));
    % Pequ2 = (Cf+Ct)' * (0.5 * Ploss);
    % Qequ2 = (Cf+Ct)' * (0.5 * Qloss);
    %     disp(Res_BaseAC{c}.mu);%shadow prices
    disp(Res_BaseAC{c}.bus(:, LAM_P));%LMP_Pd
    disp(Res_BaseAC{c}.bus(:, LAM_Q));%LMP_Qd
    disp(Res_BaseAC{c}.bus(:, MU_VMAX));%LMP_Vmax
    disp(Res_BaseAC{c}.bus(:, MU_VMIN));%LMP_Vmin
        
       
    %% Lossless Linear ACOPF
    LF_P = zeros(size(mpc{c}.bus(:,1),1),1);
    DF_P = ones(size(mpc{c}.bus(:,1),1),1);
    LF_Q = zeros(size(mpc{c}.bus(:,1),1),1);
    DF_Q = ones(size(mpc{c}.bus(:,1),1),1);
    %     tuni_less1 = clock;
    %     unit_less1 = cputime;
    [Res_LL{c},s,LMP_AC_less{c}] = Lossless_Linear_opf_LMP(mpc{c},LF_P, LF_Q, mpopt);
            %     tuni_less2 = clock;
            %     unit_less2 = cputime;
            %     tuni_less = etime(tuni_less2,tuni_less1);
            %     unit_less = unit_less2 - unit_less1;
        %
     %% iteration Lossy ACOPF model 1
    Ploss0 = zeros(size(mpc{c}.branch(:,1),1),1);
    Qloss0 = zeros(size(mpc{c}.branch(:,1),1),1);
    %mpc{c}.gen(:,2) = Res_LL{c}.gen(:,2);
    pf_results = runpf(Res_LL{c}, mpopt);
    pf_results = ext2int(pf_results);
    p_branch = (pf_results.branch(:,14) - pf_results.branch(:,16))/2/ mpc{c}.baseMVA;
    q_branch = (pf_results.branch(:,15) - pf_results.branch(:,17))/2/ mpc{c}.baseMVA;
    r = mpc{c}.branch(:,3);
    x = mpc{c}.branch(:,4);
    Ploss_temp = Ploss0;
    Ploss0 = p_branch.^2 .* r;
    %Ploss0 = real(get_losses(pf_results))/ mpc{c}.baseMVA;
    Qloss0 = q_branch.^2 .* x;
    %Qloss0 = imag(get_losses(pf_results))/ mpc{c}.baseMVA;
    while abs(sum(Ploss_temp-Ploss0))>=1e-2
        Ploss_temp = Ploss0;
        if (c<=9)
            %for Transmission systems
            [Res_LLL_IT_C{c},success,LMP_LLL_IT_C{c}, GSF, C, Ploss, Qloss] = Lossy_opf_LMP_IT_C(mpc{c}, Ploss0, Qloss0, mpopt);
        else
            %for Distribution networks
            [Res_LLL_IT_C{c},success,LMP_LLL_IT_C{c}, GSF, C, Ploss, Qloss] = Lossy_opf_DLMP_IT_C(mpc{c}, Ploss0, Qloss0, mpopt);
        end
        pf_results = runpf(Res_LLL_IT_C{c}, mpopt);
        pf_results = ext2int(pf_results);
        p_branch = (pf_results.branch(:,14) - pf_results.branch(:,16))/2/ mpc{c}.baseMVA;
        q_branch = (pf_results.branch(:,15) - pf_results.branch(:,17))/2/ mpc{c}.baseMVA;
        Ploss0 = p_branch.^2 .* r;
        Qloss0 = q_branch.^2 .* x;
    end
    Ploss_LLL_IT = Ploss0 / mpc{c}.baseMVA;
    Qloss_LLL_IT = Qloss0 / mpc{c}.baseMVA;
    
    %% iteration Lossy ACOPF model 2
    Ploss0 = zeros(size(mpc{c}.branch(:,1),1),1);
    Qloss0 = zeros(size(mpc{c}.branch(:,1),1),1);
    %mpc{c}.gen(:,2) = Res_LL{c}.gen(:,2);
    pf_results = runpf(Res_LL{c}, mpopt);
    pf_results = ext2int(pf_results);
    p_branch = (pf_results.branch(:,14) - pf_results.branch(:,16))/2/ mpc{c}.baseMVA;
    q_branch = (pf_results.branch(:,15) - pf_results.branch(:,17))/2/ mpc{c}.baseMVA;
    r = mpc{c}.branch(:,3);
    x = mpc{c}.branch(:,4);
    Ploss_temp = Ploss0;
    Ploss0 = p_branch.^2 .* r;
    Qloss0 = q_branch.^2 .* x;
    Pbus0 = real(makeSbus(pf_results.baseMVA, pf_results.bus, pf_results.gen));
    Vm0 = pf_results.bus(:,8);
    Va0 = pf_results.bus(:,9)*pi/180;
    %        T = makePTDF(mpc{c});
    %         LF_PP = 2*(r.*p_branch)' * T;
    %         LF_PP = 0;
    %         LF_PQ = 0;
    %         LF_P_Q = 0;
    %         DF_P = 1 - LF_P;
    Qbus0 = imag(makeSbus(pf_results.baseMVA,pf_results.bus, pf_results.gen));
    %LF_QQ = 2*(x.*q_branch)' * T;
    %         LF_QQ = 0;
    %         LF_QP = 0;
    %         LF_Q_P = 0;
    %         DF_Q = 1 - LF_Q;
    %         P_offset = LF_P_Q' * Pbus0 * Qbus0;
    %         Q_offset = LF_Q_P' * Qbus0 * Pbus0;
    it = 0;
    while (abs(sum(Ploss_temp-Ploss0))>=1e-2&&it<=4)
        Ploss_temp = Ploss0;
        if (c<=9)
            %for Transmission systems
            [Res_LLL_IT1{c},success,LMP_LLL_IT1{c}, GSF, C, Ploss_1, Qloss_1] = Lossy_opf_LMP_IT(c,mpc{c},Va0, Vm0, Pbus0, Ploss0, Qbus0, Qloss0, mpopt);
        else
            %for Distribution networks
            [Res_LLL_IT1{c},success,LMP_LLL_IT1{c}, GSF, C, Ploss_1, Qloss_1] = Lossy_opf_DLMP_IT(mpc{c},Va0, Vm0, Pbus0, Ploss0, Qbus0, Qloss0, mpopt);
        end
        % mpc{c}.gen(:,2) = Res_LLL_IT{c}.gen(:,2);
        pf_results = runpf(Res_LLL_IT1{c}, mpopt);
        pf_results = ext2int(pf_results);
        p_branch = (pf_results.branch(:,14) - pf_results.branch(:,16))/2/ mpc{c}.baseMVA;
        q_branch = (pf_results.branch(:,15) - pf_results.branch(:,17))/2/ mpc{c}.baseMVA;
        r = mpc{c}.branch(:,3);
        x = mpc{c}.branch(:,4);
        Ploss0 = p_branch.^2 .* r;
        Qloss0 = q_branch.^2 .* x;
        %             LF_P = 2*(r.*p_branch)' * T;%GSF.GSF_PP_T
        %             LF_P = LF_P';
        %             DF_P = 1 - LF_P;
        %         pf_results = ext2int(pf_results);
        Pbus0 = real(makeSbus(pf_results.baseMVA, pf_results.bus, pf_results.gen));
        Vm0 = pf_results.bus(:,8);
        Va0 = pf_results.bus(:,9)*pi/180;
        %             LF_Q = 2*(x.*q_branch)' * T;%* GSF.GSF_PQ_T
        %             LF_Q = LF_Q';
        %             DF_Q = 1 - LF_Q;
        Qbus0 = imag(makeSbus(pf_results.baseMVA, pf_results.bus, pf_results.gen));
        it = it + 1;
    end
    Ploss_LLL_IT = Ploss0 / mpc{c}.baseMVA;
    Qloss_LLL_IT = Qloss0 / mpc{c}.baseMVA;
      
    %%    Our non-iteration Lossy ACOPF
    if (c<=9)
        %   for Transmission systems
        [Res_LLL_non{c},success,LMP_LLL_non{c}, Ploss, Qloss] = Lossy_SOCP_opf_LMP(mpc{c}, mpopt);
    else
        %   for Distribution networks
        [Res_LLL_non{c},success,LMP_LLL_non{c}, Ploss, Qloss] = Lossy_SOCP_opf_DLMP(mpc{c}, mpopt);
    end
    %
     P_err_1 = abs(Ploss_base -  Ploss_1);
     P_err_2 = abs(Ploss_base -  Ploss);
     Q_err_1 = abs(Qloss_base -  Qloss_1);
     Q_err_2 = abs(Qloss_base -  Qloss);
     
    %% 记录LMP误差
    % Our Lossy ACOPF
    error{c}.LLL_non_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_non{c}.Total_Active_Price'),Inf);
    error{c}.LLL_non_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_non{c}.Total_Active_Price').^2));
    error{c}.LLL_non_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_non{c}.Total_Active_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    error{c}.LLL_non_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_non{c}.Total_Active_Price')));
    error{c}.LLL_non_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_non{c}.Total_Active_Price')./Res_BaseAC{c}.bus(:, LAM_P)')));
    error{c}.LLL_non_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_non{c}.Total_Active_Price').^2)));
    error{c}.LLL_non_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_non{c}.Total_Active_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    %     % It Lossy ACOPF Var M3
    %     error{c}.ACM3_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT2{c}.Total_Active_Price'),Inf);
    %     error{c}.ACM3_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT2{c}.Total_Active_Price').^2));
    %     error{c}.ACM3_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT2{c}.Total_Active_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    %     error{c}.ACM3_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT2{c}.Total_Active_Price')));
    %     error{c}.ACM3_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT2{c}.Total_Active_Price')./Res_BaseAC{c}.bus(:, LAM_P)')));
    %     error{c}.ACM3_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT2{c}.Total_Active_Price').^2)));
    %     error{c}.ACM3_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT2{c}.Total_Active_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    % It Lossy ACOPF Var M2
    error{c}.ACM2_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT1{c}.Total_Active_Price'),Inf);
    error{c}.ACM2_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT1{c}.Total_Active_Price').^2));
    error{c}.ACM2_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT1{c}.Total_Active_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    error{c}.ACM2_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT1{c}.Total_Active_Price')));
    error{c}.ACM2_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT1{c}.Total_Active_Price')./Res_BaseAC{c}.bus(:, LAM_P)')));
    error{c}.ACM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT1{c}.Total_Active_Price').^2)));
    error{c}.ACM2_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT1{c}.Total_Active_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    % It Lossy ACOPF Con M1
    error{c}.ACM1_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT_C{c}.Total_Active_Price'),Inf);
    error{c}.ACM1_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT_C{c}.Total_Active_Price').^2));
    error{c}.ACM1_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT_C{c}.Total_Active_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    error{c}.ACM1_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT_C{c}.Total_Active_Price')));
    error{c}.ACM1_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT_C{c}.Total_Active_Price')./Res_BaseAC{c}.bus(:, LAM_P)')));
    error{c}.ACM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT_C{c}.Total_Active_Price').^2)));
    error{c}.ACM1_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_LLL_IT_C{c}.Total_Active_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    % Lossless ACOPF
    error{c}.ACless_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_AC_less{c}.Total_Active_Price'),Inf);
    error{c}.ACless_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_AC_less{c}.Total_Active_Price').^2));
    error{c}.ACless_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_AC_less{c}.Total_Active_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    error{c}.ACless_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP_AC_less{c}.Total_Active_Price')));
    error{c}.ACless_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_AC_less{c}.Total_Active_Price')./Res_BaseAC{c}.bus(:, LAM_P)')));
    error{c}.ACless_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_AC_less{c}.Total_Active_Price').^2)));
    error{c}.ACless_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP_AC_less{c}.Total_Active_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    %     % It Lossy DCOPF Var M2
    %     error{c}.DCM2_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total'),Inf);
    %     error{c}.DCM2_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2));
    %     error{c}.DCM2_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    %     error{c}.DCM2_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total')));
    %     error{c}.DCM2_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total')./Res_BaseAC{c}.bus(:, LAM_P)')));
    %     error{c}.DCM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2)));
    %     error{c}.DCM2_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    %     % It Lossy DCOPF Con M1
    %     error{c}.DCM1_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total'),Inf);
    %     error{c}.DCM1_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2));
    %     error{c}.DCM1_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    %     error{c}.DCM1_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total')));
    %     error{c}.DCM1_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total')./Res_BaseAC{c}.bus(:, LAM_P)')));
    %     error{c}.DCM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2)));
    %     error{c}.DCM1_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    %     % Lossless DCOPF
    %     error{c}.DCless_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total'),Inf);
    %     error{c}.DCless_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2));
    %     error{c}.DCless_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    %     error{c}.DCless_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total')));
    %     error{c}.DCless_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total')./Res_BaseAC{c}.bus(:, LAM_P)')));
    %     error{c}.DCless_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2)));
    %     error{c}.DCless_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    
    %% 记录QLMP误差
    % Our Lossy ACOPF
    error{c}.QLLL_non_maxe = norm((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_non{c}.Total_Reactive_Price'),Inf);
    error{c}.QLLL_non_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_non{c}.Total_Reactive_Price').^2));
    error{c}.QLLL_non_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_non{c}.Total_Reactive_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_Q),1));
    error{c}.QLLL_non_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_non{c}.Total_Reactive_Price')));
    error{c}.QLLL_non_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_non{c}.Total_Reactive_Price')./Res_BaseAC{c}.bus(:, LAM_Q)')));
    error{c}.QLLL_non_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_non{c}.Total_Reactive_Price').^2)));
    error{c}.QLLL_non_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_non{c}.Total_Reactive_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_Q)' - mean(Res_BaseAC{c}.bus(:, LAM_Q)')).^2));
    %     % It Lossy ACOPF Var M3
    %     error{c}.QACM3_maxe = norm((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT2{c}.Total_Reactive_Price'),Inf);
    %     error{c}.QACM3_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT2{c}.Total_Reactive_Price').^2));
    %     error{c}.QACM3_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT2{c}.Total_Reactive_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_Q),1));
    %     error{c}.QACM3_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT2{c}.Total_Reactive_Price')));
    %     error{c}.QACM3_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT2{c}.Total_Reactive_Price')./Res_BaseAC{c}.bus(:, LAM_Q)')));
    %     error{c}.QACM3_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT2{c}.Total_Reactive_Price').^2)));
    %     error{c}.QACM3_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT2{c}.Total_Reactive_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_Q)' - mean(Res_BaseAC{c}.bus(:, LAM_Q)')).^2));
    % It Lossy ACOPF Var M2
    error{c}.QACM2_maxe = norm((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT1{c}.Total_Reactive_Price'),Inf);
    error{c}.QACM2_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT1{c}.Total_Reactive_Price').^2));
    error{c}.QACM2_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT1{c}.Total_Reactive_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_Q),1));
    error{c}.QACM2_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT1{c}.Total_Reactive_Price')));
    error{c}.QACM2_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT1{c}.Total_Reactive_Price')./Res_BaseAC{c}.bus(:, LAM_Q)')));
    error{c}.QACM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT1{c}.Total_Reactive_Price').^2)));
    error{c}.QACM2_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT1{c}.Total_Reactive_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_Q)' - mean(Res_BaseAC{c}.bus(:, LAM_Q)')).^2));
    % It Lossy ACOPF Con M1
    error{c}.QACM1_maxe = norm((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT_C{c}.Total_Reactive_Price'),Inf);
    error{c}.QACM1_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT_C{c}.Total_Reactive_Price').^2));
    error{c}.QACM1_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT_C{c}.Total_Reactive_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_Q),1));
    error{c}.QACM1_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT_C{c}.Total_Reactive_Price')));
    error{c}.QACM1_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT_C{c}.Total_Reactive_Price')./Res_BaseAC{c}.bus(:, LAM_Q)')));
    error{c}.QACM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT_C{c}.Total_Reactive_Price').^2)));
    error{c}.QACM1_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_LLL_IT_C{c}.Total_Reactive_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_Q)' - mean(Res_BaseAC{c}.bus(:, LAM_Q)')).^2));
    % Lossless ACOPF
    error{c}.QACless_maxe = norm((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_AC_less{c}.Total_Reactive_Price'),Inf);
    error{c}.QACless_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_AC_less{c}.Total_Reactive_Price').^2));
    error{c}.QACless_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_AC_less{c}.Total_Reactive_Price').^2))/size(Res_BaseAC{c}.bus(:, LAM_Q),1));
    error{c}.QACless_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_AC_less{c}.Total_Reactive_Price')));
    error{c}.QACless_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_AC_less{c}.Total_Reactive_Price')./Res_BaseAC{c}.bus(:, LAM_Q)')));
    error{c}.QACless_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_AC_less{c}.Total_Reactive_Price').^2)));
    error{c}.QACless_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_Q)' - LMP_AC_less{c}.Total_Reactive_Price').^2) / sum((Res_BaseAC{c}.bus(:, LAM_Q)' - mean(Res_BaseAC{c}.bus(:, LAM_Q)')).^2));
    %     % It Lossy DCOPF Var M2
    %     error{c}.DCM2_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total'),Inf);
    %     error{c}.DCM2_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2));
    %     error{c}.DCM2_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    %     error{c}.DCM2_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total')));
    %     error{c}.DCM2_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total')./Res_BaseAC{c}.bus(:, LAM_P)')));
    %     error{c}.DCM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2)));
    %     error{c}.DCM2_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP2{c}.total').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    %     % It Lossy DCOPF Con M1
    %     error{c}.DCM1_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total'),Inf);
    %     error{c}.DCM1_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2));
    %     error{c}.DCM1_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    %     error{c}.DCM1_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total')));
    %     error{c}.DCM1_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total')./Res_BaseAC{c}.bus(:, LAM_P)')));
    %     error{c}.DCM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2)));
    %     error{c}.DCM1_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP1{c}.total').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    %     % Lossless DCOPF
    %     error{c}.DCless_maxe = norm((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total'),Inf);
    %     error{c}.DCless_sse = mean(sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2));
    %     error{c}.DCless_mse = mean((sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2))/size(Res_BaseAC{c}.bus(:, LAM_P),1));
    %     error{c}.DCless_mae = mean(mean(abs(Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total')));
    %     error{c}.DCless_mape = mean(mean(abs((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total')./Res_BaseAC{c}.bus(:, LAM_P)')));
    %     error{c}.DCless_rmse = mean(sqrt(mean((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2)));
    %     error{c}.DCless_r2 = 1 - (sum((Res_BaseAC{c}.bus(:, LAM_P)' - LMP0{c}.total').^2) / sum((Res_BaseAC{c}.bus(:, LAM_P)' - mean(Res_BaseAC{c}.bus(:, LAM_P)')).^2));
    
    %% 记录成本误差
    error{c}.Cost_LLL_non_maxe = norm((Res_BaseAC{c}.f - Res_LLL_non{c}.f),Inf);
    error{c}.Cost_LLL_non_sse = mean(sum((Res_BaseAC{c}.f - Res_LLL_non{c}.f).^2));
    error{c}.Cost_LLL_non_mse = mean((sum((Res_BaseAC{c}.f - Res_LLL_non{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    error{c}.Cost_LLL_non_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LLL_non{c}.f)));
    error{c}.Cost_LLL_non_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LLL_non{c}.f)./Res_BaseAC{c}.f)));
    error{c}.Cost_LLL_non_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LLL_non{c}.f).^2)));
    error{c}.Cost_LLL_non_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LLL_non{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    %     % It Lossy ACOPF Var M3
    %     error{c}.Cost_ACM3_maxe = norm((Res_BaseAC{c}.f - Res_LLL_IT2{c}.f),Inf);
    %     error{c}.Cost_ACM3_sse = mean(sum((Res_BaseAC{c}.f - Res_LLL_IT2{c}.f).^2));
    %     error{c}.Cost_ACM3_mse = mean((sum((Res_BaseAC{c}.f - Res_LLL_IT2{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    %     error{c}.Cost_ACM3_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LLL_IT2{c}.f)));
    %     error{c}.Cost_ACM3_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LLL_IT2{c}.f)./Res_BaseAC{c}.f)));
    %     error{c}.Cost_ACM3_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LLL_IT2{c}.f).^2)));
    %     error{c}.Cost_ACM3_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LLL_IT2{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    % It Lossy ACOPF Var M2
    error{c}.Cost_ACM2_maxe = norm((Res_BaseAC{c}.f - Res_LLL_IT1{c}.f),Inf);
    error{c}.Cost_ACM2_sse = mean(sum((Res_BaseAC{c}.f - Res_LLL_IT1{c}.f).^2));
    error{c}.Cost_ACM2_mse = mean((sum((Res_BaseAC{c}.f - Res_LLL_IT1{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    error{c}.Cost_ACM2_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LLL_IT1{c}.f)));
    error{c}.Cost_ACM2_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LLL_IT1{c}.f)./Res_BaseAC{c}.f)));
    error{c}.Cost_ACM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LLL_IT1{c}.f).^2)));
    error{c}.Cost_ACM2_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LLL_IT1{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    % It Lossy ACOPF Con M1
    error{c}.Cost_ACM1_maxe = norm((Res_BaseAC{c}.f - Res_LLL_IT_C{c}.f),Inf);
    error{c}.Cost_ACM1_sse = mean(sum((Res_BaseAC{c}.f - Res_LLL_IT_C{c}.f).^2));
    error{c}.Cost_ACM1_mse = mean((sum((Res_BaseAC{c}.f - Res_LLL_IT_C{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    error{c}.Cost_ACM1_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LLL_IT_C{c}.f)));
    error{c}.Cost_ACM1_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LLL_IT_C{c}.f)./Res_BaseAC{c}.f)));
    error{c}.Cost_ACM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LLL_IT_C{c}.f).^2)));
    error{c}.Cost_ACM1_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LLL_IT_C{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    % Lossless ACOPF
    error{c}.Cost_ACless_maxe = norm((Res_BaseAC{c}.f - Res_LL{c}.f),Inf);
    error{c}.Cost_ACless_sse = mean(sum((Res_BaseAC{c}.f - Res_LL{c}.f).^2));
    error{c}.Cost_ACless_mse = mean((sum((Res_BaseAC{c}.f - Res_LL{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    error{c}.Cost_ACless_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LL{c}.f)));
    error{c}.Cost_ACless_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LL{c}.f)./Res_BaseAC{c}.f)));
    error{c}.Cost_ACless_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LL{c}.f).^2)));
    error{c}.Cost_ACless_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LL{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    % It Lossy DCOPF Var M2
    %     error{c}.Cost_DCM2_maxe = norm((Res_BaseAC{c}.f - Res_LDC2{c}.f),Inf);
    %     error{c}.Cost_DCM2_sse = mean(sum((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2));
    %     error{c}.Cost_DCM2_mse = mean((sum((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    %     error{c}.Cost_DCM2_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LDC2{c}.f)));
    %     error{c}.Cost_DCM2_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LDC2{c}.f)./Res_BaseAC{c}.f)));
    %     error{c}.Cost_DCM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2)));
    %     error{c}.Cost_DCM2_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    %     % It Lossy DCOPF Con M1
    %     error{c}.Cost_DCM1_maxe = norm((Res_BaseAC{c}.f - Res_LDC1{c}.f),Inf);
    %     error{c}.Cost_DCM1_sse = mean(sum((Res_BaseAC{c}.f - Res_LDC1{c}.f).^2));
    %     error{c}.Cost_DCM1_mse = mean((sum((Res_BaseAC{c}.f - Res_LDC1{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    %     error{c}.Cost_DCM1_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LDC1{c}.f)));
    %     error{c}.Cost_DCM1_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LDC1{c}.f)./Res_BaseAC{c}.f)));
    %     error{c}.Cost_DCM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2)));
    %     error{c}.Cost_DCM1_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LDC1{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    %     % Lossless DCOPF
    %     error{c}.Cost_DCless_maxe = norm((Res_BaseAC{c}.f - Res_LDC0{c}.f),Inf);
    %     error{c}.Cost_DCless_sse = mean(sum((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2));
    %     error{c}.Cost_DCless_mse = mean((sum((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    %     error{c}.Cost_DCless_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LDC0{c}.f)));
    %     error{c}.Cost_DCless_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LDC0{c}.f)./Res_BaseAC{c}.f)));
    %     error{c}.Cost_DCless_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2)));
    %     error{c}.Cost_DCless_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    
    %% 记录线路潮流误差
    error{c}.Ploss_non_ae = abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_non{c}.branch(:,[14])+Res_LLL_non{c}.branch(:,[16])));
    error{c}.Ploss_non_mae = mean(abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_non{c}.branch(:,[14])+Res_LLL_non{c}.branch(:,[16]))));
    error{c}.Ploss_non_rmse = sqrt(mean(error{c}.Ploss_non_ae.^2));
    error{c}.Pf_LLL_non_maxe = norm((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_non{c}.branch(:,[14])),Inf);
    error{c}.Pf_LLL_non_sse = mean(sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_non{c}.branch(:,[14])).^2));
    error{c}.Pf_LLL_non_mse = mean((sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_non{c}.branch(:,[14])).^2))/size(Res_BaseAC{c}.branch(:,[14]),1));
    error{c}.Pf_LLL_non_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[14]) - Res_LLL_non{c}.branch(:,[14]))));
    error{c}.Pf_LLL_non_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_non{c}.branch(:,[14]))./Res_BaseAC{c}.branch(:,[14]))));
    error{c}.Pf_LLL_non_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_non{c}.branch(:,[14])).^2)));
    error{c}.Pf_LLL_non_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_non{c}.branch(:,[14])).^2) / sum((Res_BaseAC{c}.branch(:,[14]) - mean(Res_BaseAC{c}.branch(:,[14]))).^2));
    error{c}.Pt_LLL_non_maxe = norm((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_non{c}.branch(:,[16])),Inf);
    error{c}.Pt_LLL_non_sse = mean(sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_non{c}.branch(:,[16])).^2));
    error{c}.Pt_LLL_non_mse = mean((sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_non{c}.branch(:,[16])).^2))/size(Res_BaseAC{c}.branch(:,[16]),1));
    error{c}.Pt_LLL_non_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[16]) - Res_LLL_non{c}.branch(:,[16]))));
    error{c}.Pt_LLL_non_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_non{c}.branch(:,[16]))./Res_BaseAC{c}.branch(:,[16]))));
    error{c}.Pt_LLL_non_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_non{c}.branch(:,[16])).^2)));
    error{c}.Pt_LLL_non_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_non{c}.branch(:,[16])).^2) / sum((Res_BaseAC{c}.branch(:,[16]) - mean(Res_BaseAC{c}.branch(:,[16]))).^2));
    error{c}.Qloss_non_ae = abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_non{c}.branch(:,[15])+Res_LLL_non{c}.branch(:,[17])));
    error{c}.Qloss_non_mae = mean(abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_non{c}.branch(:,[15])+Res_LLL_non{c}.branch(:,[17]))));
    error{c}.Qloss_non_rmse = sqrt(mean(error{c}.Qloss_non_ae.^2));
    error{c}.Qf_LLL_non_maxe = norm((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_non{c}.branch(:,[15])),Inf);
    error{c}.Qf_LLL_non_sse = mean(sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_non{c}.branch(:,[15])).^2));
    error{c}.Qf_LLL_non_mse = mean((sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_non{c}.branch(:,[15])).^2))/size(Res_BaseAC{c}.branch(:,[15]),1));
    error{c}.Qf_LLL_non_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[15]) - Res_LLL_non{c}.branch(:,[15]))));
    error{c}.Qf_LLL_non_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_non{c}.branch(:,[15]))./Res_BaseAC{c}.branch(:,[15]))));
    error{c}.Qf_LLL_non_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_non{c}.branch(:,[15])).^2)));
    error{c}.Qf_LLL_non_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_non{c}.branch(:,[15])).^2) / sum((Res_BaseAC{c}.branch(:,[15]) - mean(Res_BaseAC{c}.branch(:,[15]))).^2));
    error{c}.Qt_LLL_non_maxe = norm((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_non{c}.branch(:,[17])),Inf);
    error{c}.Qt_LLL_non_sse = mean(sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_non{c}.branch(:,[17])).^2));
    error{c}.Qt_LLL_non_mse = mean((sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_non{c}.branch(:,[17])).^2))/size(Res_BaseAC{c}.branch(:,[17]),1));
    error{c}.Qt_LLL_non_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[17]) - Res_LLL_non{c}.branch(:,[17]))));
    error{c}.Qt_LLL_non_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_non{c}.branch(:,[17]))./Res_BaseAC{c}.branch(:,[17]))));
    error{c}.Qt_LLL_non_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_non{c}.branch(:,[17])).^2)));
    error{c}.Qt_LLL_non_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_non{c}.branch(:,[17])).^2) / sum((Res_BaseAC{c}.branch(:,[17]) - mean(Res_BaseAC{c}.branch(:,[17]))).^2));
    
    % It Lossy ACOPF Var M3
    %         error{c}.Ploss_ACM3_ae = abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_IT2{c}.branch(:,[14])+Res_LLL_IT2{c}.branch(:,[16])));
    %         error{c}.Ploss_ACM3_mae = mean(abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_IT2{c}.branch(:,[14])+Res_LLL_IT2{c}.branch(:,[16]))));
    %         error{c}.Ploss_ACM3_rmse = sqrt(mean(error{c}.Ploss_ACM3_ae.^2));
    %         error{c}.Pf_ACM3_ae = abs(Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14]));
    %         error{c}.Pf_ACM3_maxe = norm((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14])),Inf);
    %         error{c}.Pf_ACM3_sse = mean(sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14])).^2));
    %         error{c}.Pf_ACM3_mse = mean((sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14])).^2))/size(Res_BaseAC{c}.branch(:,[14]),1));
    %         error{c}.Pf_ACM3_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14]))));
    %         error{c}.Pf_ACM3_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14]))./Res_BaseAC{c}.branch(:,[14]))));
    %         error{c}.Pf_ACM3_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14])).^2)));
    %         error{c}.Pf_ACM3_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT2{c}.branch(:,[14])).^2) / sum((Res_BaseAC{c}.branch(:,[14]) - mean(Res_BaseAC{c}.branch(:,[14]))).^2));
    %         error{c}.Pt_ACM3_ae = abs(Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16]));
    %         error{c}.Pt_ACM3_maxe = norm((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16])),Inf);
    %         error{c}.Pt_ACM3_sse = mean(sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16])).^2));
    %         error{c}.Pt_ACM3_mse = mean((sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16])).^2))/size(Res_BaseAC{c}.branch(:,[16]),1));
    %         error{c}.Pt_ACM3_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16]))));
    %         error{c}.Pt_ACM3_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16]))./Res_BaseAC{c}.branch(:,[16]))));
    %         error{c}.Pt_ACM3_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16])).^2)));
    %         error{c}.Pt_ACM3_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT2{c}.branch(:,[16])).^2) / sum((Res_BaseAC{c}.branch(:,[16]) - mean(Res_BaseAC{c}.branch(:,[16]))).^2));
    %         error{c}.Qloss_ACM3_ae = abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_IT2{c}.branch(:,[15])+Res_LLL_IT2{c}.branch(:,[17])));
    %         error{c}.Qloss_ACM3_mae = mean(abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_IT2{c}.branch(:,[15])+Res_LLL_IT2{c}.branch(:,[17]))));
    %         error{c}.Qloss_ACM3_rmse = sqrt(mean(error{c}.Qloss_ACM3_ae.^2));
    %         error{c}.Qf_ACM3_ae = abs(Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15]));
    %         error{c}.Qf_ACM3_maxe = norm((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15])),Inf);
    %         error{c}.Qf_ACM3_sse = mean(sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15])).^2));
    %         error{c}.Qf_ACM3_mse = mean((sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15])).^2))/size(Res_BaseAC{c}.branch(:,[15]),1));
    %         error{c}.Qf_ACM3_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15]))));
    %         error{c}.Qf_ACM3_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15]))./Res_BaseAC{c}.branch(:,[15]))));
    %         error{c}.Qf_ACM3_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15])).^2)));
    %         error{c}.Qf_ACM3_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT2{c}.branch(:,[15])).^2) / sum((Res_BaseAC{c}.branch(:,[15]) - mean(Res_BaseAC{c}.branch(:,[15]))).^2));
    %         error{c}.Qt_ACM3_ae = abs(Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17]));
    %         error{c}.Qt_ACM3_maxe = norm((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17])),Inf);
    %         error{c}.Qt_ACM3_sse = mean(sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17])).^2));
    %         error{c}.Qt_ACM3_mse = mean((sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17])).^2))/size(Res_BaseAC{c}.branch(:,[17]),1));
    %         error{c}.Qt_ACM3_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17]))));
    %         error{c}.Qt_ACM3_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17]))./Res_BaseAC{c}.branch(:,[17]))));
    %         error{c}.Qt_ACM3_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17])).^2)));
    %         error{c}.Qt_ACM3_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT2{c}.branch(:,[17])).^2) / sum((Res_BaseAC{c}.branch(:,[17]) - mean(Res_BaseAC{c}.branch(:,[17]))).^2));
    %         error{c}.Branch_ACM3_mae = ((error{c}.Pf_ACM3_mae + error{c}.Pt_ACM3_mae) / 2 + (error{c}.Qf_ACM3_mae + error{c}.Qt_ACM3_mae) / 2) / 2;
    %         error{c}.Branch_ACM3_rmse = ((error{c}.Pf_ACM3_rmse + error{c}.Pt_ACM3_rmse) / 2 + (error{c}.Qf_ACM3_rmse + error{c}.Qt_ACM3_rmse) / 2) / 2;
    %         error{c}.Branch_ACM3_r2 = ((error{c}.Pf_ACM3_r2 + error{c}.Pt_ACM3_r2) / 2 + (error{c}.Qf_ACM3_r2 + error{c}.Qt_ACM3_r2) / 2) / 2;
    
    % It Lossy ACOPF Var M2
    error{c}.Ploss_ACM2_ae = abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_IT1{c}.branch(:,[14])+Res_LLL_IT1{c}.branch(:,[16])));
    error{c}.Ploss_ACM2_mae = mean(abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_IT1{c}.branch(:,[14])+Res_LLL_IT1{c}.branch(:,[16]))));
    error{c}.Ploss_ACM2_rmse = sqrt(mean(error{c}.Ploss_ACM2_ae.^2));
    error{c}.Pf_ACM2_ae = abs(Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14]));
    error{c}.Pf_ACM2_maxe = norm((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14])),Inf);
    error{c}.Pf_ACM2_sse = mean(sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14])).^2));
    error{c}.Pf_ACM2_mse = mean((sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14])).^2))/size(Res_BaseAC{c}.branch(:,[14]),1));
    error{c}.Pf_ACM2_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14]))));
    error{c}.Pf_ACM2_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14]))./Res_BaseAC{c}.branch(:,[14]))));
    error{c}.Pf_ACM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14])).^2)));
    error{c}.Pf_ACM2_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT1{c}.branch(:,[14])).^2) / sum((Res_BaseAC{c}.branch(:,[14]) - mean(Res_BaseAC{c}.branch(:,[14]))).^2));
    error{c}.Pt_ACM2_ae = abs(Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16]));
    error{c}.Pt_ACM2_maxe = norm((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16])),Inf);
    error{c}.Pt_ACM2_sse = mean(sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16])).^2));
    error{c}.Pt_ACM2_mse = mean((sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16])).^2))/size(Res_BaseAC{c}.branch(:,[16]),1));
    error{c}.Pt_ACM2_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16]))));
    error{c}.Pt_ACM2_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16]))./Res_BaseAC{c}.branch(:,[16]))));
    error{c}.Pt_ACM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16])).^2)));
    error{c}.Pt_ACM2_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT1{c}.branch(:,[16])).^2) / sum((Res_BaseAC{c}.branch(:,[16]) - mean(Res_BaseAC{c}.branch(:,[16]))).^2));
    error{c}.Qloss_ACM2_ae = abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_IT1{c}.branch(:,[15])+Res_LLL_IT1{c}.branch(:,[17])));
    error{c}.Qloss_ACM2_mae = mean(abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_IT1{c}.branch(:,[15])+Res_LLL_IT1{c}.branch(:,[17]))));
    error{c}.Qloss_ACM2_rmse = sqrt(mean(error{c}.Qloss_ACM2_ae.^2));
    error{c}.Qf_ACM2_ae = abs(Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15]));
    error{c}.Qf_ACM2_maxe = norm((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15])),Inf);
    error{c}.Qf_ACM2_sse = mean(sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15])).^2));
    error{c}.Qf_ACM2_mse = mean((sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15])).^2))/size(Res_BaseAC{c}.branch(:,[15]),1));
    error{c}.Qf_ACM2_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15]))));
    error{c}.Qf_ACM2_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15]))./Res_BaseAC{c}.branch(:,[15]))));
    error{c}.Qf_ACM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15])).^2)));
    error{c}.Qf_ACM2_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT1{c}.branch(:,[15])).^2) / sum((Res_BaseAC{c}.branch(:,[15]) - mean(Res_BaseAC{c}.branch(:,[15]))).^2));
    error{c}.Qt_ACM2_ae = abs(Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17]));
    error{c}.Qt_ACM2_maxe = norm((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17])),Inf);
    error{c}.Qt_ACM2_sse = mean(sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17])).^2));
    error{c}.Qt_ACM2_mse = mean((sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17])).^2))/size(Res_BaseAC{c}.branch(:,[17]),1));
    error{c}.Qt_ACM2_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17]))));
    error{c}.Qt_ACM2_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17]))./Res_BaseAC{c}.branch(:,[17]))));
    error{c}.Qt_ACM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17])).^2)));
    error{c}.Qt_ACM2_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT1{c}.branch(:,[17])).^2) / sum((Res_BaseAC{c}.branch(:,[17]) - mean(Res_BaseAC{c}.branch(:,[17]))).^2));
    error{c}.Branch_ACM2_mae = ((error{c}.Pf_ACM2_mae + error{c}.Pt_ACM2_mae) / 2 + (error{c}.Qf_ACM2_mae + error{c}.Qt_ACM2_mae) / 2) / 2;
    error{c}.Branch_ACM2_rmse = ((error{c}.Pf_ACM2_rmse + error{c}.Pt_ACM2_rmse) / 2 + (error{c}.Qf_ACM2_rmse + error{c}.Qt_ACM2_rmse) / 2) / 2;
    error{c}.Branch_ACM2_r2 = ((error{c}.Pf_ACM2_r2 + error{c}.Pt_ACM2_r2) / 2 + (error{c}.Qf_ACM2_r2 + error{c}.Qt_ACM2_r2) / 2) / 2;
    
    % It Lossy ACOPF Con M1
    error{c}.Ploss_ACM1_ae = abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_IT_C{c}.branch(:,[14])+Res_LLL_IT_C{c}.branch(:,[16])));
    error{c}.Ploss_ACM1_mae = mean(abs((Res_BaseAC{c}.branch(:,[14])+Res_BaseAC{c}.branch(:,[16])) - (Res_LLL_IT_C{c}.branch(:,[14])+Res_LLL_IT_C{c}.branch(:,[16]))));
    error{c}.Ploss_ACM1_rmse = sqrt(mean(error{c}.Ploss_ACM1_ae.^2));
    error{c}.Pf_ACM1_ae = abs(Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14]));
    error{c}.Pf_ACM1_maxe = norm((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14])),Inf);
    error{c}.Pf_ACM1_sse = mean(sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14])).^2));
    error{c}.Pf_ACM1_mse = mean((sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14])).^2))/size(Res_BaseAC{c}.branch(:,[14]),1));
    error{c}.Pf_ACM1_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14]))));
    error{c}.Pf_ACM1_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14]))./Res_BaseAC{c}.branch(:,[14]))));
    error{c}.Pf_ACM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14])).^2)));
    error{c}.Pf_ACM1_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[14]) - Res_LLL_IT_C{c}.branch(:,[14])).^2) / sum((Res_BaseAC{c}.branch(:,[14]) - mean(Res_BaseAC{c}.branch(:,[14]))).^2));
    error{c}.Pt_ACM1_ae = abs(Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16]));
    error{c}.Pt_ACM1_maxe = norm((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16])),Inf);
    error{c}.Pt_ACM1_sse = mean(sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16])).^2));
    error{c}.Pt_ACM1_mse = mean((sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16])).^2))/size(Res_BaseAC{c}.branch(:,[16]),1));
    error{c}.Pt_ACM1_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16]))));
    error{c}.Pt_ACM1_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16]))./Res_BaseAC{c}.branch(:,[16]))));
    error{c}.Pt_ACM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16])).^2)));
    error{c}.Pt_ACM1_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[16]) - Res_LLL_IT_C{c}.branch(:,[16])).^2) / sum((Res_BaseAC{c}.branch(:,[16]) - mean(Res_BaseAC{c}.branch(:,[16]))).^2));
    error{c}.Qloss_ACM1_ae = abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_IT_C{c}.branch(:,[15])+Res_LLL_IT_C{c}.branch(:,[17])));
    error{c}.Qloss_ACM1_mae = mean(abs((Res_BaseAC{c}.branch(:,[15])+Res_BaseAC{c}.branch(:,[17])) - (Res_LLL_IT_C{c}.branch(:,[15])+Res_LLL_IT_C{c}.branch(:,[17]))));
    error{c}.Qloss_ACM1_rmse = sqrt(mean(error{c}.Qloss_ACM1_ae.^2));
    error{c}.Qf_ACM1_ae = abs(Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15]));
    error{c}.Qf_ACM1_maxe = norm((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15])),Inf);
    error{c}.Qf_ACM1_sse = mean(sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15])).^2));
    error{c}.Qf_ACM1_mse = mean((sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15])).^2))/size(Res_BaseAC{c}.branch(:,[15]),1));
    error{c}.Qf_ACM1_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15]))));
    error{c}.Qf_ACM1_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15]))./Res_BaseAC{c}.branch(:,[15]))));
    error{c}.Qf_ACM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15])).^2)));
    error{c}.Qf_ACM1_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[15]) - Res_LLL_IT_C{c}.branch(:,[15])).^2) / sum((Res_BaseAC{c}.branch(:,[15]) - mean(Res_BaseAC{c}.branch(:,[15]))).^2));
    error{c}.Qt_ACM1_ae = abs(Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17]));
    error{c}.Qt_ACM1_maxe = norm((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17])),Inf);
    error{c}.Qt_ACM1_sse = mean(sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17])).^2));
    error{c}.Qt_ACM1_mse = mean((sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17])).^2))/size(Res_BaseAC{c}.branch(:,[17]),1));
    error{c}.Qt_ACM1_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17]))));
    error{c}.Qt_ACM1_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17]))./Res_BaseAC{c}.branch(:,[17]))));
    error{c}.Qt_ACM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17])).^2)));
    error{c}.Qt_ACM1_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[17]) - Res_LLL_IT_C{c}.branch(:,[17])).^2) / sum((Res_BaseAC{c}.branch(:,[17]) - mean(Res_BaseAC{c}.branch(:,[17]))).^2));
    error{c}.Branch_ACM1_mae = ((error{c}.Pf_ACM1_mae + error{c}.Pt_ACM1_mae) / 2 + (error{c}.Qf_ACM1_mae + error{c}.Qt_ACM1_mae) / 2) / 2;
    error{c}.Branch_ACM1_rmse = ((error{c}.Pf_ACM1_rmse + error{c}.Pt_ACM1_rmse) / 2 + (error{c}.Qf_ACM1_rmse + error{c}.Qt_ACM1_rmse) / 2) / 2;
    error{c}.Branch_ACM1_r2 = ((error{c}.Pf_ACM1_r2 + error{c}.Pt_ACM1_r2) / 2 + (error{c}.Qf_ACM1_r2 + error{c}.Qt_ACM1_r2) / 2) / 2;
    
    % Lossless ACOPF
    error{c}.Pf_ACless_ae = abs(Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14]));
    error{c}.Pf_ACless_maxe = norm((Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14])),Inf);
    error{c}.Pf_ACless_sse = mean(sum((Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14])).^2));
    error{c}.Pf_ACless_mse = mean((sum((Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14])).^2))/size(Res_BaseAC{c}.branch(:,[14]),1));
    error{c}.Pf_ACless_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14]))));
    error{c}.Pf_ACless_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14]))./Res_BaseAC{c}.branch(:,[14]))));
    error{c}.Pf_ACless_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14])).^2)));
    error{c}.Pf_ACless_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[14]) - Res_LL{c}.branch(:,[14])).^2) / sum((Res_BaseAC{c}.branch(:,[14]) - mean(Res_BaseAC{c}.branch(:,[14]))).^2));
    error{c}.Pf_ACless_ae = abs(Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16]));
    error{c}.Pt_ACless_maxe = norm((Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16])),Inf);
    error{c}.Pt_ACless_sse = mean(sum((Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16])).^2));
    error{c}.Pt_ACless_mse = mean((sum((Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16])).^2))/size(Res_BaseAC{c}.branch(:,[16]),1));
    error{c}.Pt_ACless_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16]))));
    error{c}.Pt_ACless_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16]))./Res_BaseAC{c}.branch(:,[16]))));
    error{c}.Pt_ACless_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16])).^2)));
    error{c}.Pt_ACless_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[16]) - Res_LL{c}.branch(:,[16])).^2) / sum((Res_BaseAC{c}.branch(:,[16]) - mean(Res_BaseAC{c}.branch(:,[16]))).^2));
    error{c}.Qf_ACless_ae = abs(Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15]));
    error{c}.Qf_ACless_maxe = norm((Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15])),Inf);
    error{c}.Qf_ACless_sse = mean(sum((Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15])).^2));
    error{c}.Qf_ACless_mse = mean((sum((Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15])).^2))/size(Res_BaseAC{c}.branch(:,[15]),1));
    error{c}.Qf_ACless_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15]))));
    error{c}.Qf_ACless_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15]))./Res_BaseAC{c}.branch(:,[15]))));
    error{c}.Qf_ACless_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15])).^2)));
    error{c}.Qf_ACless_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[15]) - Res_LL{c}.branch(:,[15])).^2) / sum((Res_BaseAC{c}.branch(:,[15]) - mean(Res_BaseAC{c}.branch(:,[15]))).^2));
    error{c}.Qt_ACless_ae = abs(Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17]));
    error{c}.Qt_ACless_maxe = norm((Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17])),Inf);
    error{c}.Qt_ACless_sse = mean(sum((Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17])).^2));
    error{c}.Qt_ACless_mse = mean((sum((Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17])).^2))/size(Res_BaseAC{c}.branch(:,[17]),1));
    error{c}.Qt_ACless_mae = mean(mean(abs(Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17]))));
    error{c}.Qt_ACless_mape = mean(mean(abs((Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17]))./Res_BaseAC{c}.branch(:,[17]))));
    error{c}.Qt_ACless_rmse = mean(sqrt(mean((Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17])).^2)));
    error{c}.Qt_ACless_r2 = 1 - (sum((Res_BaseAC{c}.branch(:,[17]) - Res_LL{c}.branch(:,[17])).^2) / sum((Res_BaseAC{c}.branch(:,[17]) - mean(Res_BaseAC{c}.branch(:,[17]))).^2));
    error{c}.Branch_ACless_mae = ((error{c}.Pf_ACless_mae + error{c}.Pt_ACless_mae) / 2 + (error{c}.Qf_ACless_mae + error{c}.Qt_ACless_mae) / 2) / 2;
    error{c}.Branch_ACless_rmse = ((error{c}.Pf_ACless_rmse + error{c}.Pt_ACless_rmse) / 2 + (error{c}.Qf_ACless_rmse + error{c}.Qt_ACless_rmse) / 2) / 2;
    error{c}.Branch_ACless_r2 = ((error{c}.Pf_ACless_r2 + error{c}.Pt_ACless_r2) / 2 + (error{c}.Qf_ACless_r2 + error{c}.Qt_ACless_r2) / 2) / 2;
    
    %         % It Lossy DCOPF Var M2
    %         error.Cost_DCM2_maxe = norm((Res_BaseAC{c}.f - Res_LDC2{c}.f),Inf);
    %         error.Cost_DCM2_sse = mean(sum((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2));
    %         error.Cost_DCM2_mse = mean((sum((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    %         error.Cost_DCM2_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LDC2{c}.f)));
    %         error.Cost_DCM2_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LDC2{c}.f)./Res_BaseAC{c}.f)));
    %         error.Cost_DCM2_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2)));
    %         error.Cost_DCM2_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    %         % It Lossy DCOPF Con M1
    %         error.Cost_DCM1_maxe = norm((Res_BaseAC{c}.f - Res_LDC1{c}.f),Inf);
    %         error.Cost_DCM1_sse = mean(sum((Res_BaseAC{c}.f - Res_LDC1{c}.f).^2));
    %         error.Cost_DCM1_mse = mean((sum((Res_BaseAC{c}.f - Res_LDC1{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    %         error.Cost_DCM1_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LDC1{c}.f)));
    %         error.Cost_DCM1_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LDC1{c}.f)./Res_BaseAC{c}.f)));
    %         error.Cost_DCM1_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LDC2{c}.f).^2)));
    %         error.Cost_DCM1_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LDC1{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
    %         % Lossless DCOPF
    %         error.Cost_DCless_maxe = norm((Res_BaseAC{c}.f - Res_LDC0{c}.f),Inf);
    %         error.Cost_DCless_sse = mean(sum((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2));
    %         error.Cost_DCless_mse = mean((sum((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2))/size(Res_BaseAC{c}.f,1));
    %         error.Cost_DCless_mae = mean(mean(abs(Res_BaseAC{c}.f - Res_LDC0{c}.f)));
    %         error.Cost_DCless_mape = mean(mean(abs((Res_BaseAC{c}.f - Res_LDC0{c}.f)./Res_BaseAC{c}.f)));
    %         error.Cost_DCless_rmse = mean(sqrt(mean((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2)));
    %         error.Cost_DCless_r2 = 1 - (sum((Res_BaseAC{c}.f - Res_LDC0{c}.f).^2) / sum((Res_BaseAC{c}.f - mean(Res_BaseAC{c}.f)).^2));
end