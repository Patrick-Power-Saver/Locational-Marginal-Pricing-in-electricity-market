function [MVAbase, bus, gen, gencost, branch, f, success, et] = Lossless_Linear_opf_LMP(mpc, LF_P, LF_Q, varargin)
% Lossless Linear ACOPF
% 11-05-2022

%% 导入系统数据
define_constants;
mpopt = mpoption;
[baseMVA, bus, gen, branch, gencost] = loadcase(mpc);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[ref, pv, pq] = bustypes(bus, gen); % reference bus index
t0 = clock;

%% 参数设置
branch(:, BR_STATUS) = 1;
br = find(branch(:, BR_STATUS));            %% in-service branches
no =bus(:,1);
nb = size(bus, 1);          % number of buses
nl = size(branch, 1);       % number of lines
ng = size(gen, 1);          % number of gens
wbuses = find((bus(:,2) == PV |bus(:,2) == PQ ) & (bus(:,PD) > 0));%bus(:,2) == REF |

Pd = bus(:, PD);              % bus load (MVA)
Qd = bus(:, QD);
Fmax = branch(:, RATE_A);                       % branch flow limit (MVA)

Vmin = bus(:, VMIN);                            % minimum voltage magnitude (p.u.)
Vmax = bus(:, VMAX);                            % maximum voltage magnitude (p.u.)
Sgmin = gen(:, PMIN) + 1i * gen(:, QMIN);       % gen max. power output (MVA)
Sgmax = gen(:, PMAX) + 1i * gen(:, QMAX);       % gen max. power output (MVA)
C2 = gencost(:, COST);                          % gen injection cost ($/MWh^2)
C1 = gencost(:, COST + 1);                      % gen injection cost ($/MWh)
C0 = gencost(:, COST + 2);                      % gen injection cost ($/h)

% 节点-支路关联矩阵和节点导纳矩阵
gbus = gen(:, GEN_BUS);
Cg = full(sparse(gbus, (1:ng)', 1, nb, ng));
Cw = sparse(wbuses, (1:length(wbuses))', 1, nb, length(wbuses));
f = branch(:, F_BUS);
t = branch(:, T_BUS);
i = [(1:nl)'; (1:nl)'];             % element k, j is -1 if branch k connects "to" bus j
Cf = full(sparse((1:nl)', f, 1, nl, nb));
Ct = full(sparse((1:nl)', t, 1, nl, nb));
Cft = full(sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb));
Ctf = sparse(i, [t; f], [ones(nl, 1); -ones(nl, 1)], nl, nb);

% H = makePTDF(mpc);
% Hg = H(find(branch(:,RATE_A) > 0), gbus);
% Hd = H(find(branch(:,RATE_A) > 0),find(bus(:,PD) > 0));
% Pd_hat = Pd(find(bus(:,PD) > 0));Qd_hat = Qd(find(bus(:,PD) > 0));

[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA;
Ysc = 1 ./ (branch(:, BR_R) - 1j * branch(:, BR_X));
Gbus = real(Ybus);
Bbus = imag(Ybus);
Gsh = real(Ysh);
Bsh = imag(Ysh);

GP = Gbus; %导纳矩阵的实部
BD = diag(sum(Bbus));
B_mod= Bbus-BD;% B'
BP = -B_mod; %BP=-B'
BQ = -Bbus; %BQ=-B
GQ = -Gbus;                       %GQ approxiately equals -Gbus

Xp_dlpf = [BP GP];
Xq_dlpf = [GQ BQ];
CJ = full([Xp_dlpf;Xq_dlpf]);
% valid_idx = setdiff(1:2*nb, ref);
% C_prime = CJ(valid_idx, valid_idx);
% CJ = pinv(C_prime);
% C_J = zeros(2*nb,2*nb);
% C_J(valid_idx,valid_idx) = CJ;
% CJ = C_J;
CJ(ref,ref) = 0;
CJ = pinv(CJ);
C_PVm = CJ(nb+1:end,1:nb);
C_QVm = CJ(nb+1:end,nb+1:end);
C_PVa = CJ(1:nb,1:nb);
C_QVa = CJ(1:nb,nb+1:end);

B11 = BP([pv;pq],[pv;pq]);  %delta的行对应B11的列 (PV+PQ) 没有slackbus
G12 = GP([pv;pq],pq);%G   V的行对应G12的列（PQ）
G21 = GQ(pq,[pv;pq]);%-G  delta的行对应G21的列 （PV+PQ）
B22 = BQ(pq,pq); %B22  V的行对应B22的列 （PQ）
t1 =G12/B22; % Ninv(L) = G12*inv(B22)  G12的列与B22的行一致 （PQ）
t2 =G21/B11; % Minv(H) = G21*inv(B11)  G21的列与B11的行一致 （PV+PQ）
B11m = B11 - t1*G21;
B22m = B22 - t2*G12;
XT_dlpf = [inv(B11m) -inv(B11m)*G12*inv(B22)];
XV_dlpf = [-inv(B22m)*G21*inv(B11) inv(B22m)];


% Xp_pqpv = [B11 G12];
% Xq_pqpv = [G21 B22];
% CJ_pqpv = [Xp_pqpv;Xq_pqpv];

%单位换算
%-------- unit conversion --------- %
Pd = Pd / baseMVA;
Qd = Qd / baseMVA;
Fmax = Fmax / baseMVA;          % p.u.
% FLmax = sparse([f;t],[t,f],[Fmax;-Fmax],nb,nb);
Sgmin = Sgmin / baseMVA;        % p.u.
Sgmax = Sgmax / baseMVA;        % p.u.
C2 = C2 * (baseMVA^2);          % $/h
C1 = C1 * baseMVA;              % $/h

%% 定义变量
Vm = sdpvar(nb, 1, 'full');
Vm_2 = sdpvar(nb, 1, 'full');
Vmij_2 = sdpvar(nl,1, 'full');
abs_Vmij_2 = sdpvar(nl,1, 'full');
Vaij_2 = sdpvar(nl,1, 'full');
abs_Vaacij_2 = sdpvar(nl,1, 'full');
Va = sdpvar(nb,1, 'full');
Pbus = sdpvar(nb,1, 'full');
Qbus = sdpvar(nb,1, 'full');
Pm1 = sdpvar(size(bus([pv;pq]), 1),1, 'full');
Qm1 = sdpvar(size(bus([pq]), 1),1, 'full');
Pequ = sdpvar(nb,1, 'full');
Qequ = sdpvar(nb,1, 'full');
Pf = sdpvar(nl,1, 'full');
Qf = sdpvar(nl,1, 'full');
Pt = sdpvar(nl,1, 'full');
Qt = sdpvar(nl,1, 'full');
Ploss = sdpvar(nl,1, 'full');
Qloss = sdpvar(nl,1, 'full');
Pg = sdpvar(ng, 1, 'full');
Qg = sdpvar(ng, 1, 'full');
Pres_d = sdpvar(ng, 1, 'full');
Pres_u = sdpvar(ng, 1, 'full');

% 为节点电压赋初值
% -------- assign initial value -------- %
% Matpower solution
% mpopt = mpoption( 'out.all', 0);
% result = runpf(mpc, mpopt);
% assign(Vm_2, result.bus(:,VM).^2);
% assign(Va, result.bus(:,VA)*pi/180);
% Sbus = makeSbus(result.baseMVA, result.bus, result.gen);
% assign(Pbus, real(Sbus));
% assign(Qbus, imag(Sbus));
% assign(Pg, result.gen(:,PG));
% assign(Qg, result.gen(:,QG));

%% 目标函数
%obj = C1' * Pg + C1' * (Qg);         % linear cost function ($/h)
obj = C2' * (Pg.^2) + C1' * Pg + ones(ng, 1)' * C0;    % + C2' * (Qg.^2) + C1' * Qg + ones(ng, 1)' * C0

%% 约束条件
Constraints = [];

%% A 节点功率平衡等式约束
Constraints = [Constraints; (Pbus + Pd - Cg * Pg == 0):'Balance_P'];
Constraints = [Constraints; (Qbus + Qd - Cg * Qg == 0):'Balance_Q'];

%% B 节点注入与节点电压的线性等式约束（B,C二选一）
% for i = 1:nb
%     Fi = f(find(f == i));
%     Ti = t(find(f == i));
%     Constraints=[Constraints;Pbus(i,1) == BP(i,Ti) * (Va(Fi)-Va(Ti)) + GP(i,Ti) * (0.5 * (Vm_2(Fi) + Vm_2(Ti)))];
%     Constraints=[Constraints;Qbus(i,1) == BQ(i,Ti) * (0.5 * (Vm_2(Fi) + Vm_2(Ti))) + GQ(i,Ti) * (Va(Fi)-Va(Ti))];
% end
% Constraints=[Constraints;Pbus == Xp * [Va(f)-Va(t);0.5*(Vm_2(f)+Vm_2(t))]];
% Constraints=[Constraints;Qbus == Xq * [Va(f)-Va(t);0.5*(Vm_2(f)+Vm_2(t))]];
Constraints=[Constraints;Pbus == Xp_dlpf * [Va;Vm_2]];
Constraints=[Constraints;Qbus == Xq_dlpf * [Va;Vm_2]];
%Constraints = [Constraints; [Pbus;Qbus] == CJ * [Va;Vm]];
% 节点注入功率P和Q的元素和维度是已知的
% Constraints = [Constraints;Pm1 == Pbus([pv;pq]) - GP([pv;pq],ref)*bus(ref,VM) - GP([pv;pq],pv)*bus(pv,VM)];%节点注入有功功率
% Constraints = [Constraints;Qm1 == Qbus(pq) - BQ(pq,ref)*bus(ref,VM) - BQ(pq,pv)*bus(pv,VM)];%节点注入有功功率
% Constraints = [Constraints; Va([pv;pq]) == XT_dlpf * [Pm1;Qm1]];
% Constraints = [Constraints; Vm_2(pq) == XV_dlpf * [Pm1;Qm1]];

%Constraints = [Constraints; [Va;Vm] == inv(CJ) * [(Cg * Pg- Pd);(Cg * Qg- Qd)]];
%Constraints = [Constraints; [Va;Vm] == CJ \ [(Cg * Pg- Pd);(Cg * Qg- Qd)]];
% Constraints = [Constraints; Vm == CJ(:,nb+1:end) \ [Pbus;Qbus]];
% Constraints = [Constraints; Va == CJ(:,1:nb) \ [Pbus;Qbus]];

% Constraints = [Constraints; (Vm == CJ(nb+1:end,1:nb) * Pbus + CJ(nb+1:end,nb+1:end) * Qbus):'Voltage'];
% Constraints = [Constraints; (Va == CJ(1:nb,1:nb) * Pbus + CJ(1:nb,nb+1:end) * Qbus):'Angle'];

% Constraints = [Constraints; Vm == C.C_PVm * Pbus + C.C_QVm * Qbus];
% Constraints = [Constraints; Va == C.C_PVa * Pbus + C.C_QVa * Qbus];

% Constraints = [Constraints; (Vm_2 == C_PVm * Pbus + C_QVm * Qbus):'V'];
% Constraints = [Constraints; (Va == C_PVa * Pbus + C_QVa * Qbus):'VA'];

%% C 节点注入与支路功率的线性等式约束（B,C二选一）
% Constraints=[Constraints,Pbus == Cft' * Pf + sum(Gbus,2).*Vm_2 ];
% Constraints=[Constraints,Qbus == Cft' * Qf - sum(Bbus,2).*Vm_2 ];
% Constraints=[Constraints,Pbus == Ctf' * Pt + sum(Gbus,2).*Vm_2 ];
% Constraints=[Constraints,Qbus == Ctf' * Qt - sum(Bbus,2).*Vm_2 ];

%% D 支路功率与节点电压的等式约束（D1,D2二选一,注意Vm或Vm_2要与B、C一致）
% D1 (Vi - Vj)
% Constraints=[Constraints; Pf == (Vm(f)-Vm(t)).* real(Ysc) - (Va(f)-Va(t)).*imag(Ysc) ];
% Constraints=[Constraints; Qf == -(Vm(f)-Vm(t)).* imag(Ysc) - (Va(f)-Va(t)).*real(Ysc)];
% Constraints=[Constraints; Pt == (Vm(t)-Vm(f)).* real(Ysc) - (Va(t)-Va(f)).*imag(Ysc)];
% Constraints=[Constraints; Qt == -(Vm(t)-Vm(f)).* imag(Ysc) - (Va(t)-Va(f)).*real(Ysc)];
% 
GSF_PP_F = (C_PVm(f,:) - C_PVm(t,:)) .* real(Ysc) - (C_PVa(f,:) - C_PVa(t,:)) .* imag(Ysc);
GSF_PQ_F = (C_QVm(f,:) - C_QVm(t,:)) .* real(Ysc) - (C_QVa(f,:) - C_QVa(t,:)) .* imag(Ysc);
GSF_QP_F = -(C_PVa(f,:) - C_PVa(t,:)) .* real(Ysc) - (C_PVm(f,:) - C_PVm(t,:)) .* imag(Ysc); 
GSF_QQ_F = -(C_QVa(f,:) - C_QVa(t,:)) .* real(Ysc) - (C_QVm(f,:) - C_QVm(t,:)) .* imag(Ysc);
GSF_PP_T = -GSF_PP_F;
GSF_PQ_T = -GSF_PQ_F;
GSF_QP_T = -GSF_QP_F;
GSF_QQ_T = -GSF_QQ_F;
Constraints=[Constraints; Pf == (GSF_PP_F * Pbus + GSF_PQ_F * Qbus)];
Constraints=[Constraints; Qf == (GSF_QP_F * Pbus + GSF_QQ_F * Qbus)];
Constraints=[Constraints; Pt == (GSF_PP_T * Pbus + GSF_PQ_T * Qbus)];
Constraints=[Constraints; Qt == (GSF_QP_T * Pbus + GSF_QQ_T * Qbus)];

% Constraints=[Constraints; Pf == ((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* real(Ysc) - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; Qf == -((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* imag(Ysc) - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*real(Ysc)];
% Constraints=[Constraints; Pt == ((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* real(Ysc) - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; Qt == -((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* imag(Ysc) - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*real(Ysc)];

% D2 (Vi.^2 - Vj.^2)/2
% Constraints=[Constraints; Pf  == (Vm_2(f)-Vm_2(t)).* real(Ysc) / 2 - (Va(f)-Va(t)).*imag(Ysc)];
% Constraints=[Constraints; Qf  == -(Vm_2(f)-Vm_2(t)).* imag(Ysc)/ 2 - (Va(f)-Va(t)).*real(Ysc)];
% Constraints=[Constraints; Pt  == (Vm_2(t)-Vm_2(f)).* real(Ysc)/ 2 - (Va(t)-Va(f)).*imag(Ysc)];
% Constraints=[Constraints; Qt  == -(Vm_2(t)-Vm_2(f)).* imag(Ysc)/ 2 - (Va(t)-Va(f)).*real(Ysc)];

% Constraints=[Constraints; -Pf == ((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* real(Ysc)  / 2 - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; -Qf == -((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* imag(Ysc)  / 2 - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*real(Ysc)];
% Constraints=[Constraints; -Pt == ((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* real(Ysc)  / 2 - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; -Qt == -((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* imag(Ysc)  / 2 - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*real(Ysc)];

% %% 网损约束
% Constraints=[Constraints;Pequ == 0.5 * Cf' * Ploss];
% Constraints=[Constraints;Qequ == 0.5 * Cf' * Qloss];
% Constraints=[Constraints; Ploss >= 0];
% Constraints=[Constraints; Qloss >= 0];
% Constraints=[Constraints; Ploss == (Vmij_2+Vaij_2).* real(Ysc)];
% Constraints=[Constraints; Qloss == (Vmij_2+Vaij_2).* imag(Ysc)];
% % M1
% for i = 1:nb
%     Constraints=[Constraints; (Vm_2(f(i)) - Vm_2(t(i)))^2 <= Vmij_2(i)];
% end
% for i = 1:nb
%     Constraints=[Constraints; (Va(f(i)) - Va(t(i)))^2 <= Vaij_2(i)];
% end
% % M2
% Constraints=[Constraints; (Vm_2(f) - Vm_2(t)).^2 <= Vmij_2];
% Constraints=[Constraints; (Va(f) - Va(t))^2 <= Vaij_2];
% % M3
% Constraints=[Constraints; diag((Vm_2(f) - Vm_2(t))) * (Vm_2(f) - Vm_2(t)) <= Vmij_2];
% Constraints=[Constraints; diag((Va(f) - Va(t))) * (Va(f) - Va(t)) <= Vaij_2];

%% 电压上下限不等式约束（Vm或Vm_2与BCD保持一致）
% Constraints = [Constraints; (Vmin <= Vm):'Voltage_low'];         % voltage magnitude (p.u.)
% Constraints = [Constraints; (Vm <= Vmax):'Voltage_high'];         % voltage magnitude (p.u.)
%Constraints=[Constraints, (Vmin <= Vm <= Vmax):'Voltage'];
Constraints=[Constraints, Vm_2([ref], 1) == bus([ref],8).^2];                                         
Constraints=[Constraints, (Vmin.^2 <= Vm_2):'Voltage_low'];
Constraints=[Constraints, (Vm_2 <= Vmax.^2):'Voltage_high'];

%% 机组出力上下限不等式约束
% Constraints = [Constraints; real(Sgmin) <= Pg <= real(Sgmax)];                                  % complex power output of generator (p.u.)
% Constraints = [Constraints; imag(Sgmin) <= Qg <= imag(Sgmax)];

Constraints = [Constraints; (real(Sgmin) <= Pg):'Pg_low'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (Pg <= real(Sgmax)):'Pg_high'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (imag(Sgmin) <= Qg):'Qg_low'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (Qg <= imag(Sgmax)):'Qg_high'];   

% Constraints = [Constraints; (real(Sgmin) + Pres_d <= Pg):'Pg_min'];
% Constraints = [Constraints; (real(Sgmax) - Pres_u >= Pg):'Pg_max']; 
% Constraints = [Constraints; (Pres_d <= 0.1 * real(Sgmax)):'Res_d_max'];
% Constraints = [Constraints; (Pres_u <= 0.1 * real(Sgmax)):'Res_u_max'];
% Constraints = [Constraints; (Pres_d >= 0):'Res_d_min'];
% Constraints = [Constraints; (Pres_u >= 0):'Res_u_min'];
% Constraints = [Constraints; (imag(Sgmin) <= Qg):'Qg_min'];% + var.Qres_d
% Constraints = [Constraints; (imag(Sgmax) >= Qg):'Qg_max']; % - var.Qres_u

%% 平衡节点电压相角为0
Constraints = [Constraints; Va(ref) == 0];                                  % reference bus angle (rad)

%% 线路容量不等式约束（配电网不需要此约束）
t_poly=10;
if (nb <= nl)
%     Constraints = [Constraints,  Pf.^2 + Qf.^2 <= Fmax.^2];                  % branch flow limit at "from" end (p.u.)
%     Constraints = [Constraints,  Pt.^2 + Qt.^2 <= Fmax.^2];
    for i = 1:size(Pf,1)
    Constraints = [Constraints,  (Pf(i)^2 + Qf(i)^2 <= Fmax(i)^2):'Con_Sf'];                  % branch flow limit at "from" end (p.u.)
    Constraints = [Constraints,  (Pt(i)^2 + Qt(i)^2 <= Fmax(i)^2):'Con_St'];
    end
% 
%      Constraints = [Constraints,  (-Fmax <= (GSF_PP_F * Pbus + GSF_PQ_F * Qbus)):'Con_PF_low'];                  % branch flow limit at "from" end (p.u.)
%      Constraints = [Constraints,  ((GSF_PP_F * Pbus + GSF_PQ_F * Qbus) <= Fmax):'Con_PF_high'];                  % branch flow limit at "from" end (p.u.)
%      Constraints = [Constraints,  (-Fmax <= (GSF_QP_F * Pbus + GSF_QQ_F * Qbus)):'Con_QF_low'];
%      Constraints = [Constraints,  ((GSF_QP_F * Pbus + GSF_QQ_F * Qbus) <= Fmax):'Con_QF_high'];
% 
%      Constraints = [Constraints,  (Fmax >= (GSF_PP_T * Pbus + GSF_PQ_T * Qbus)):'Con_PT_low'];                  % branch flow limit at "from" end (p.u.)
%      Constraints = [Constraints,  ((GSF_PP_T * Pbus + GSF_PQ_T * Qbus) >= -Fmax):'Con_PT_high'];                  % branch flow limit at "from" end (p.u.)
%      Constraints = [Constraints,  (Fmax >= (GSF_QP_T * Pbus + GSF_QQ_T * Qbus)):'Con_QT_low'];
%      Constraints = [Constraints,  ((GSF_QP_T * Pbus + GSF_QQ_T * Qbus) >= -Fmax):'Con_QT_high'];

%  for i=1:t_poly
%      theta_poly = pi * (i-1) / (2*t_poly);
%         Constraints = [Constraints,(-Fmax<= cos(theta_poly) * Pf + sin(theta_poly) * Qf <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * Pf + cos(theta_poly) * Qf <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= cos(theta_poly) * Pt + sin(theta_poly) * Qt <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * Pt + cos(theta_poly) * Qt <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * Pt + cos(theta_poly) * Qt <= Fmax):['Con_T_high',num2str(i)]];
%     
%         Constraints = [Constraints,(-Fmax<= cos(theta_poly) * Pf + sin(theta_poly) * Qf)&(cos(theta_poly) * Pf + sin(theta_poly) * Qf <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * Pf + cos(theta_poly) * Qf)&(-sin(theta_poly) * Pf + cos(theta_poly) * Qf <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= cos(theta_poly) * Pt + sin(theta_poly) * Qt)&(cos(theta_poly) * Pt + sin(theta_poly) * Qt <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * Pt + cos(theta_poly) * Qt)&(-sin(theta_poly) * Pt + cos(theta_poly) * Qt <= Fmax)];
% 
% %     Constraints = [Constraints,(-2*Fmax<= cos(theta_poly) * (GSF_PP * Pbus + GSF_PQ * Qbus) + sin(theta_poly) * (GSF_QP * Pbus + GSF_QQ * Qbus)...
% %         -sin(theta_poly) * (GSF_PP * Pbus + GSF_PQ * Qbus) + cos(theta_poly) * (GSF_QP * Pbus + GSF_QQ * Qbus)):['Con_F_low',num2str(i)]];
% %     Constraints = [Constraints,(cos(theta_poly) * (GSF_PP * Pbus + GSF_PQ * Qbus) + sin(theta_poly) * (GSF_QP * Pbus + GSF_QQ * Qbus)...
% %         -sin(theta_poly) * (GSF_PP * Pbus + GSF_PQ * Qbus) + cos(theta_poly) * (GSF_QP * Pbus + GSF_QQ * Qbus)<= 2*Fmax):['Con_F_high',num2str(i)]];
% %     
% %     Constraints = [Constraints,(-2*Fmax<= cos(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + sin(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus)...
% %         -sin(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + cos(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus)):['Con_T_low',num2str(i)]];
% %     Constraints = [Constraints,(cos(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + sin(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus)...
% %         -sin(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + cos(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus)<= 2*Fmax):['Con_T_high',num2str(i)]];
% %     
% %             Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * (GSF_PP * Pbus + GSF_PQ * Qbus) + cos(theta_poly) * (GSF_QP * Pbus + GSF_QQ * Qbus) <= Fmax)];
% %             Constraints = [Constraints,(-Fmax<= cos(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + sin(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus) <= Fmax)];
% %             Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + cos(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus) <= Fmax)];
% %             Constraints = [Constraints,(-Fmax<= cos(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + sin(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus) <= Fmax)];
% %             Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * (-GSF_PP * Pbus - GSF_PQ * Qbus) + cos(theta_poly) * (-GSF_QP * Pbus - GSF_QQ * Qbus) <= Fmax)];
% end
else
end

%% 求解问题
opt = sdpsettings('solver','+gurobi','verbose',0,'gurobi.NumericFocus',3);% ipopt gurobi cplex 
opt.gurobi.QCPDual = 1;
varargout=optimize(Constraints,obj,opt);
varargout.info

%% 检查约束违反
check(Constraints);

%% 输出结果
% -------- post-processing -------- %
%res.Ploss = sum(value(Ploss)) / baseMVA;
res.Vm = sqrt(value(Vm_2));
res.Va = value(Va)/pi*180;
res.Pbus = value(Pbus);
res.Qbus = value(Qbus);
res.Pf = value(Pf)* baseMVA;
res.Qf = value(Qf)* baseMVA;
res.Pt = value(Pt)* baseMVA;
res.Qt = value(Qt)* baseMVA;
res.Pg = value(Pg)* baseMVA;         % gen real power output (MW)
res.Qg = value(Qg)* baseMVA;
res.cost = value(obj);      % gen real power output cost ($/h)
% LMP：Energy
res.lmp_energy_P = dual(Constraints('Balance_P'))/ baseMVA;   % lmp($/MWh)
res.lmp_energy_Q = dual(Constraints('Balance_Q'))/ baseMVA;   % lmp($/MWh)
LMP.PEnergy_Price = res.lmp_energy_P;
LMP.QEnergy_Price = res.lmp_energy_Q;
% LMP：V
% res.lmp_V = dual(Constraints('V'));   % lmp($/MWh)
% res.lmp_VA = dual(Constraints('VA'));   % lmp($/MWh)
res.lmp_V_mul = dual(Constraints('Voltage_low'));   % lmp($/MWh)
res.lmp_V_muh = dual(Constraints('Voltage_high'));   % lmp($/MWh)
LMP.PVl_Price = C_PVm * res.lmp_V_mul;
LMP.PVh_Price = C_PVm * res.lmp_V_muh;
LMP.QVl_Price = C_QVm * res.lmp_V_mul;
LMP.QVh_Price = C_QVm * res.lmp_V_muh;
LMP.PV_Price = LMP.PVh_Price - LMP.PVl_Price;
LMP.QV_Price = LMP.QVh_Price - LMP.QVl_Price;

%res.lmp_mul_F = 0;res.lmp_muh_F = 0;res.lmp_mul_T = 0;res.lmp_muh_T = 0;
% for i = 1:t_poly
%     res.lmp_mul_F = res.lmp_mul_F + dual(Constraints(['Con_F_low', num2str(i)]));
%     res.lmp_muh_F = res.lmp_muh_F + dual(Constraints(['Con_F_high', num2str(i)]));
%     res.lmp_mul_T = res.lmp_mul_T + dual(Constraints(['Con_T_low', num2str(i)]));
%     res.lmp_muh_T = res.lmp_muh_T + dual(Constraints(['Con_T_high', num2str(i)]));
% end

% LMP：Loss
% res.lmp_loss_P = -0.5 * Cf * res.lmp_energy_P;   %dual(Constraints('Loss_P'))/ baseMVA -0.5* Cf*res.lmp_energy_P lmp($/MWh)
% res.lmp_loss_Q = -0.5 * Cf * res.lmp_energy_Q;   % dual(Constraints('Loss_Q'))/ baseMVA -0.5 * Cf * res.lmp_energy_Q lmp($/MWh)
% LF_P = value(2 * (C_PVm(f,:) - C_PVm(t,:)).^2 * Pbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + 2 * (C_PVa(f,:) - C_PVa(t,:)).^2 * Pbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* real(Ysc);
% LF_Q = value(2 * (C_QVm(f,:) - C_QVm(t,:)).^2 * Qbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + 2 * (C_QVa(f,:) - C_QVa(t,:)).^2 * Qbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* imag(Ysc);
% LMP.PLoss_Price = -Cf' * LF_P .* res.lmp_energy_P;%-Cf' * LF_P .* res.lmp_energy_P (-Cf'* LF_P).*(-2*pinv(Cf)*res.lmp_loss_P)
% LMP.QLoss_Price = -Cf' * LF_Q .* res.lmp_energy_Q;%-Cf' * LF_Q .* res.lmp_energy_Q (-Cf'* LF_Q).*(-2*pinv(Cf)*res.lmp_loss_Q)


% LMP：Con

% if (nb <= nl)
res.lmp_mu_SF = dual(Constraints('Con_Sf'))/ baseMVA;   % lmp($/MWh)
res.lmp_mu_ST = dual(Constraints('Con_St'))/ baseMVA;   % lmp($/MWh)
LMP.PfCon_Price = Cf' * (value(GSF_PP_F.^2 * Pbus + 2 * GSF_PP_F .* GSF_PQ_F * Qbus + GSF_QP_F.^2 * Pbus + 2 * GSF_QP_F .* GSF_QQ_F * Qbus) .* res.lmp_mu_SF);
LMP.PtCon_Price = Ct' * (value(GSF_PP_T.^2 * Pbus + 2 * GSF_PP_T .* GSF_PQ_T * Qbus + GSF_QP_T.^2 * Pbus + 2 * GSF_QP_T .* GSF_QQ_T * Qbus) .* res.lmp_mu_ST);
LMP.QfCon_Price = Cf' * (value(GSF_PQ_F.^2 * Qbus + 2 * GSF_PP_F .* GSF_PQ_F * Pbus + GSF_QQ_F.^2 * Qbus + 2 * GSF_QP_F .* GSF_QQ_F * Pbus) .* res.lmp_mu_SF);
LMP.QtCon_Price = Ct' * (value(GSF_PQ_T.^2 * Qbus + 2 * GSF_PP_T .* GSF_PQ_T * Pbus + GSF_QQ_T.^2 * Qbus + 2 * GSF_QP_T .* GSF_QQ_T * Pbus) .* res.lmp_mu_ST);
LMP.PL = (LMP.PfCon_Price + LMP.PtCon_Price) / 2;
LMP.QL = (LMP.QfCon_Price + LMP.QtCon_Price) / 2;

%     res.lmp_mul_PF = dual(Constraints('Con_PF_low'))/ baseMVA;   % lmp($/MWh)
%     res.lmp_muh_PF = dual(Constraints('Con_PF_high'))/ baseMVA;   % lmp($/MWh)
%     res.lmp_mul_PT = dual(Constraints('Con_PT_low'))/ baseMVA;   % lmp($/MWh)
%     res.lmp_muh_PT = dual(Constraints('Con_PT_high'))/ baseMVA;   % lmp($/MWh)
%     res.lmp_mul_QF = dual(Constraints('Con_QF_low'))/ baseMVA;   % lmp($/MWh)
%     res.lmp_muh_QF = dual(Constraints('Con_QF_high'))/ baseMVA;   % lmp($/MWh)
%     res.lmp_mul_QT = dual(Constraints('Con_QT_low'))/ baseMVA;   % lmp($/MWh)
%     res.lmp_muh_QT = dual(Constraints('Con_QT_high'))/ baseMVA;   % lmp($/MWh)
%     LMP.PfCon_Price = GSF_PP_F' * res.lmp_mul_PF + GSF_PP_F' * res.lmp_muh_PF;
%     LMP.PtCon_Price = GSF_PP_T' * res.lmp_mul_PT + GSF_PP_T' * res.lmp_muh_PT;
%     LMP.QfCon_Price = GSF_QP_F' * res.lmp_mul_QF + GSF_QP_F' * res.lmp_muh_QF;
%     LMP.QtCon_Price = GSF_QP_T' * res.lmp_mul_QT + GSF_QP_T' * res.lmp_muh_QT;
% end

    LMP.Total_Active_Price = LMP.PEnergy_Price + LMP.PL + LMP.PV_Price;
    LMP.Total_Reactive_Price = LMP.QEnergy_Price + LMP.QL + LMP.QV_Price;

% LMP.Total_Active_Price = res.lmp_energy_P + (LMP.PVl_Price + LMP.PVh_Price) + (LMP.PfCon_Price + LMP.PtCon_Price) / 2;
% LMP.Total_Reactive_Price = res.lmp_energy_Q + (LMP.QVl_Price + LMP.QVh_Price) + (LMP.QfCon_Price + LMP.QtCon_Price) / 2;

% res.lmp_mul_Pg = dual(Constraints('Pg_low'))/ baseMVA; 
% res.lmp_muh_Pg = dual(Constraints('Pg_high'))/ baseMVA;
% res.lmp_mul_Qg = dual(Constraints('Qg_low'))/ baseMVA; 
% res.lmp_muh_Qg = dual(Constraints('Qg_high'))/ baseMVA;

mpc.bus(:,VM) = res.Vm;
mpc.bus(:,VA) = res.Va;
mpc.gen(:,PG) = res.Pg;
mpc.gen(:,QG) = res.Qg;
mpc.branch(:,PF) = res.Pf;
mpc.branch(:,PT) = res.Pt;
mpc.branch(:,QF) = res.Qf;
mpc.branch(:,QT) = res.Qt;
% mpc.bus(:, 14) = LMP.PEnergy_Price;%LMP_Pd
% mpc.bus(:, 15) = LMP.QEnergy_Price;%LMP_Qd
% mpc.bus(:, 16) = LMP.Vh_Price;%LMP_Vmax
% mpc.bus(:, 17) = LMP.Vl_Price;%LMP_Vmin
% mpc.branch(:, 18) = LMP.PfCon_Price + 1i * LMP.QfCon_Price;%LMP_Sf
% mpc.branch(:, 19) = LMP.PtCon_Price + 1i * LMP.QtCon_Price;%LMP_St
% mpc.gen(:, 22) = ;%LMP_Pgmax
% mpc.gen(:, 23);%LMP_Pgmin
% mpc.gen(:, 24);%LMP_Qgmax
% mpc.gen(:, 25);%LMP_Qgmin

V = mpc.bus(:,VM) .* exp(sqrt(-1) * pi/180 * mpc.bus(:,VA));
[mpc.bus, mpc.gen, mpc.branch] = pfsoln(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
[i2eac, mpc.bus, mpc.gen, mpc.branch] = ext2int(mpc.bus, mpc.gen, mpc.branch);
success = 1;
results.success = success;
results.et = etime(clock, t0);
[mpc.bus, mpc.gen, mpc.branch, mpc.f] = deal(mpc.bus, mpc.gen, mpc.branch, res.cost);
results = int2ext(mpc);

%% zero out result fields of out-of-service gens & branches
if ~isempty(results.order.gen.status.off)
    results.gen(results.order.gen.status.off, [PG QG]) = 0;
end
if ~isempty(results.order.branch.status.off)
    results.branch(results.order.branch.status.off, [PF QF PT QT]) = 0;
end

if nargout == 1 || nargout == 2 || nargout == 3
    MVAbase = results;
    bus = success;
    gen = LMP;
elseif nargout > 3
    [MVAbase, bus, gen, branch, f, et] = ...
        deal(results.baseMVA, results.bus, results.gen, results.branch, results.f, results.et);
    % else  %% don't define MVAbase, so it doesn't print anything
end
disp('done');

% %% 绘制功率圆
% % 参数设置
% R = 1; % 圆的半径
% num_squares = t_poly; % 旋转正方形的数量
% 
% % 圆的边界点
% theta_circle = linspace(0, 2 * pi, 100);
% x_circle = R * cos(theta_circle);
% y_circle = R * sin(theta_circle);
% 
% % 创建图形窗口
% figure;
% plot(x_circle, y_circle, 'b-', 'LineWidth', 1.5); % 绘制圆的边界
% hold on;
% 
% % 绘制多个旋转正方形
% for i = 1:num_squares
%     % 计算旋转角度
%     theta_ro = pi * (i-1) / (2*num_squares)
% 
%     % 定义正方形的四个顶点，未旋转时位于 x 和 y 轴方向
%     square_x = R * [-1, 1, 1, -1, -1];
%     square_y = R * [-1, -1, 1, 1, -1];
% 
%     % 旋转正方形
%     x_rotated = cos(theta_ro) * square_x - sin(theta_ro) * square_y;
%     y_rotated = sin(theta_ro) * square_x + cos(theta_ro) * square_y;
% 
%     % 绘制旋转后的正方形
%     plot(x_rotated, y_rotated, 'r-', 'LineWidth', 1.2);
% end
% 
% % 设置图形属性
% axis equal;
% xlim([-1.5, 1.5]);
% ylim([-1.5, 1.5]);
% xlabel('Rx');
% ylabel('Ry');
% title('Polygonal Approximation of Circle with Rotated Squares');
% legend('x^2 + y^2 = R^2', 'Rotated squares');
% grid on;
% hold off;
end