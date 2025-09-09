function [MVAbase, bus, gen, gencost, branch, f, success, et] = Lossy_opf_DLMP_IT_C(mpc, Ploss0,  Qloss0, varargin)
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
CJ = [Xp_dlpf;Xq_dlpf];
CJ = full([Xp_dlpf;Xq_dlpf]);
CJ(ref,ref) = 0;
CJ = pinv(CJ);
%CJ = pinv(CJ);%配网
C_PVm = CJ(nb+1:end,1:nb);
C_QVm = CJ(nb+1:end,nb+1:end);
C_PVa = CJ(1:nb,1:nb);
C_QVa = CJ(1:nb,nb+1:end);

% B11 = BP([pv;pq],[pv;pq]);  %delta的行对应B11的列 (PV+PQ) 没有slackbus
% G12 = GP([pv;pq],pq);%G   V的行对应G12的列（PQ）
% G21 = GQ(pq,[pv;pq]);%-G  delta的行对应G21的列 （PV+PQ）
% B22 = BQ(pq,pq); %B22  V的行对应B22的列 （PQ）
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
z = binvar(nl, 1, 'full'); % 引入二进制变量 
M = 1e6; % 一个足够大的数，用于表示正负无穷大
C = sdpvar(nl,1);

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
Constraints = [Constraints; (Pbus + Pd + (0.5 * (Cf+Ct)' * (0.5 * Ploss0)) - Cg * Pg == 0):'Balance_P'];
Constraints = [Constraints; (Qbus + Qd + (0.5 * (Cf+Ct)' * (0.5 * Qloss0)) - Cg * Qg == 0):'Balance_Q'];
% Constraints = [Constraints; sum(Pd) - sum(Pg)  == 0];   % nodal real power balance (MW)
% Constraints = [Constraints; sum(Qd) - sum(Qg)  == 0];   % nodal real power balance (MW)

%% B 节点注入与节点电压的线性等式约束（B,C二选一）
Constraints=[Constraints;Pbus == Xp_dlpf * [Va;Vm_2]];
Constraints=[Constraints;Qbus == Xq_dlpf * [Va;Vm_2]];
%Constraints = [Constraints; [Pbus;Qbus] == CJ * [Va;Vm]];
%Constraints = [Constraints; [Va;Vm] == CJ \ [Pbus;Qbus]];
%Constraints = [Constraints; [Va;Vm] == inv(CJ) * [(Cg * Pg- Pd);(Cg * Qg- Qd)]];
%Constraints = [Constraints; [Va;Vm] == CJ \ [(Cg * Pg- Pd);(Cg * Qg- Qd)]];
% Constraints = [Constraints; Vm == CJ(:,nb+1:end) \ [Pbus;Qbus]];
% Constraints = [Constraints; Va == CJ(:,1:nb) \ [Pbus;Qbus]];

% Constraints = [Constraints; (Vm == CJ(nb+1:end,1:nb) * Pbus + CJ(nb+1:end,nb+1:end) * Qbus):'Voltage'];
% Constraints = [Constraints; (Va == CJ(1:nb,1:nb) * Pbus + CJ(1:nb,nb+1:end) * Qbus):'Angle'];

% Constraints = [Constraints; Vm == C_PVm * Pbus + C_QVm * Qbus];
% Constraints = [Constraints; Va == C_PVa * Pbus + C_QVa * Qbus];

% Constraints = [Constraints; (Vm_2 == C_PVm * Pbus + C_QVm * Qbus):'V'];
% Constraints = [Constraints; Va == C_PVa * Pbus + C_QVa * Qbus];

%% C 节点注入与支路功率的线性等式约束（B,C二选一）
% Constraints=[Constraints,Pbus == Cft' * Pf + sum(Gbus,2).*Vm_2 ];
% Constraints=[Constraints,Qbus == Cft' * Qf - sum(Bbus,2).*Vm_2 ];
% Constraints=[Constraints,Pbus == Ctf' * Pt + sum(Gbus,2).*Vm_2 ];
% Constraints=[Constraints,Qbus == Ctf' * Qt - sum(Bbus,2).*Vm_2 ];

%% D 支路功率与节点电压的等式约束（D1,D2二选一,注意Vm或Vm_2要与B、C一致）
% D1 (Vi - Vj)
% Constraints=[Constraints; Pf == (Vm(f)-Vm(t)).* real(Ysc) - (Va(f)-Va(t)).*imag(Ysc) ];
% Constraints=[Constraints; Qf == -(Vm(f)-Vm(t)).* imag(Ysc) - (Va(f)-Va(t)).*real(Ysc)];
% Constraints=[Constraints; -Pt == (Vm(t)-Vm(f)).* real(Ysc) - (Va(t)-Va(f)).*imag(Ysc)];
% Constraints=[Constraints; -Qt == -(Vm(t)-Vm(f)).* imag(Ysc) - (Va(t)-Va(f)).*real(Ysc)];

GSF.GSF_PP_F = (C_PVm(f,:) - C_PVm(t,:)) .* real(Ysc) - (C_PVa(f,:) - C_PVa(t,:)) .* imag(Ysc);
GSF.GSF_PQ_F = (C_QVm(f,:) - C_QVm(t,:)) .* real(Ysc) - (C_QVa(f,:) - C_QVa(t,:)) .* imag(Ysc);
GSF.GSF_QP_F = -(C_PVa(f,:) - C_PVa(t,:)) .* real(Ysc) - (C_PVm(f,:) - C_PVm(t,:)) .* imag(Ysc); 
GSF.GSF_QQ_F = -(C_QVa(f,:) - C_QVa(t,:)) .* real(Ysc) - (C_QVm(f,:) - C_QVm(t,:)) .* imag(Ysc);
GSF.GSF_PP_T = -GSF.GSF_PP_F;
GSF.GSF_PQ_T = -GSF.GSF_PQ_F;
GSF.GSF_QP_T = -GSF.GSF_QP_F;
GSF.GSF_QQ_T = -GSF.GSF_QQ_F;
Constraints=[Constraints; Pf == (GSF.GSF_PP_F * (Pbus) + GSF.GSF_PQ_F * (Qbus))+ 0.5 * Ploss0];
Constraints=[Constraints; Qf == (GSF.GSF_QP_F * (Pbus) + GSF.GSF_QQ_F * (Qbus))+ 0.5 * Qloss0];
Constraints=[Constraints; Pt == (GSF.GSF_PP_T * (Pbus) + GSF.GSF_PQ_T * (Qbus))+ 0.5 * Ploss0];
Constraints=[Constraints; Qt == (GSF.GSF_QP_T * (Pbus) + GSF.GSF_QQ_T * (Qbus))+ 0.5 * Qloss0];

%% 使用McCormick包络法进行凸松弛
% % 设置上下界
% x_min = -Fmax;
% x_max = Fmax;
% y_min = -Fmax;
% y_max = Fmax;
% % 定义双线性项
% PfPt = sdpvar(nl,1);
% QfQt = sdpvar(nl,1);
% % McCormick包络法的四个不等式
% for i = 1:nl
% Constraints = [Constraints, PfPt(i) >= x_min(i).*Pt(i) + y_min(i).*Pf(i) - x_min(i).*y_min(i)];
% Constraints = [Constraints, PfPt(i) >= x_max(i).*Pt(i) + y_max(i).*Pf(i) - x_max(i).*y_max(i)];
% Constraints = [Constraints, PfPt(i) <= x_min(i).*Pt(i) + y_max(i).*Pf(i) - x_min(i).*y_max(i)];
% Constraints = [Constraints, PfPt(i) <= x_max(i).*Pt(i) + y_min(i).*Pf(i) - x_max(i).*y_min(i)];
% Constraints = [Constraints, PfPt(i) <= 0];
% Constraints = [Constraints, QfQt(i) >= x_min(i).*Qt(i) + y_min(i).*Qf(i) - x_min(i).*y_min(i)];
% Constraints = [Constraints, QfQt(i) >= x_max(i).*Qt(i) + y_max(i).*Qf(i) - x_max(i).*y_max(i)];
% Constraints = [Constraints, QfQt(i) <= x_min(i).*Qt(i) + y_max(i).*Qf(i) - x_min(i).*y_max(i)];
% Constraints = [Constraints, QfQt(i) <= x_max(i).*Qt(i) + y_min(i).*Qf(i) - x_max(i).*y_min(i)];
% Constraints = [Constraints, QfQt(i) <= 0];
% end

% Constraints=[Constraints; Pf == ((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* real(Ysc) - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; Qf == -((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* imag(Ysc) - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*real(Ysc)];
% Constraints=[Constraints; Pt == ((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* real(Ysc) - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; Qt == -((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* imag(Ysc) - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*real(Ysc)];

% D2 (Vi.^2 - Vj.^2)/2
% Constraints=[Constraints; Pf  == (Vm_2(f)-Vm_2(t)).* real(Ysc) / 2 - (Va(f)-Va(t)).*imag(Ysc) + 0.5 * Ploss];
% Constraints=[Constraints; Qf  == -(Vm_2(f)-Vm_2(t)).* imag(Ysc)/ 2 - (Va(f)-Va(t)).*real(Ysc) +  0.5 * Qloss];
% Constraints=[Constraints; Pt  == (Vm_2(t)-Vm_2(f)).* real(Ysc)/ 2 - (Va(t)-Va(f)).*imag(Ysc) +  0.5 * Ploss];
% Constraints=[Constraints; Qt  == -(Vm_2(t)-Vm_2(f)).* imag(Ysc)/ 2 - (Va(t)-Va(f)).*real(Ysc) + 0.5 * Qloss];

% Constraints=[Constraints; -Pf == ((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* real(Ysc)  / 2 - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; -Qf == -((C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus).* imag(Ysc)  / 2 - ((C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus).*real(Ysc)];
% Constraints=[Constraints; -Pt == ((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* real(Ysc)  / 2 - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*imag(Ysc)];
% Constraints=[Constraints; -Qt == -((C_PVm(t,:) - C_PVm(f,:)) * Pbus + (C_QVm(t,:) - C_QVm(f,:)) * Qbus).* imag(Ysc)  / 2 - ((C_PVa(t,:) - C_PVa(f,:)) * Pbus + (C_QVa(t,:) - C_QVa(f,:)) * Qbus).*real(Ysc)];

%% 网损约束
% LF_PP = ((C_PVm(f,:) - C_PVm(t,:)).^2 + (C_PVa(f,:) - C_PVa(t,:)).^2) .* real(Ysc);
% LF_PQ = ((C_QVm(f,:) - C_QVm(t,:)).^2 + (C_QVa(f,:) - C_QVa(t,:)).^2) .* real(Ysc);
% LF_P_Q = (2 * (C_PVm(f,:) - C_PVm(t,:)) .* (C_PVa(f,:) - C_PVa(t,:)) + 2 * (C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:))) .* real(Ysc);
% LF_QP = ((C_PVm(f,:) - C_PVm(t,:)).^2 + (C_PVa(f,:) - C_PVa(t,:)).^2) .* imag(Ysc);
% LF_QQ = ((C_QVm(f,:) - C_QVm(t,:)).^2 + (C_QVa(f,:) - C_QVa(t,:)).^2) .* imag(Ysc);
% LF_Q_P = (2 * (C_PVm(f,:) - C_PVm(t,:)) .* (C_PVa(f,:) - C_PVa(t,:)) + 2 * (C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:))) .* imag(Ysc);
% Constraints=[Constraints;(Ploss == LF_PP * (Pbus0 .* Pbus) + LF_PQ * (Qbus0 .* Qbus) + LF_P_Q * (Pbus0 .* Qbus0)):'Ploss'];%LF_P_Q' * Pbus0 * Qbus0 =  P_offset
% Constraints=[Constraints;(Qloss == LF_QP * (Pbus0 .* Pbus) + LF_QQ * (Qbus0 .* Qbus) + LF_Q_P * (Qbus0 .* Pbus0)):'Qloss'];%LF_Q_P' * Qbus0 * Pbus0 =  Q_offset

% LF_PVa = ((Va0(f,:) - Va0(t,:)) .* real(Ysc));
% LF_PVm_2 =0.25 * ((Vm0(f,:).^2 - Vm0(t,:).^2) .* real(Ysc));
% %offset_P = -((Va0(f,:) - Va0(t,:)).^2 + (Vm0(f,:) - Vm0(t,:)).^2) .* real(Ysc);
% 
% LF_QVa = ((Va0(f,:) - Va0(t,:)) .* imag(Ysc));
% LF_QVm_2 = 0.25 * ((Vm0(f,:).^2 - Vm0(t,:).^2) .* imag(Ysc));
% %offset_Q =  ((Va0(f,:) - Va0(t,:)).^2 + (Vm0(f,:) - Vm0(t,:)).^2) .* imag(Ysc);
% Constraints=[Constraints;(Ploss == LF_PVa .* (Va(f,:) - Va(t,:)) + LF_PVm_2 .* (Vm_2(f,:) - Vm_2(t,:))):'Ploss'];%LF_P_Q' * Pbus0 * Qbus0 =  P_offset
% Constraints=[Constraints;(Qloss == LF_QVa .* (Va(f,:) - Va(t,:)) + LF_QVm_2 .* (Vm_2(f,:) - Vm_2(t,:)) ):'Qloss'];%LF_Q_P' * Qbus0 * Pbus0 =  Q_offset

% Constraints=[Constraints; 0 <= Ploss <= 0.03*Fmax];
% Constraints=[Constraints; 0 <= Qloss <= 0.03*Fmax];
% Constraints=[Constraints; 0 <= Ploss <= abs(Pf)];
% Constraints=[Constraints; 0 <= Qloss <= abs(Qf)];

% Pf_abs = sdpvar(nl,1);
% Qf_abs = sdpvar(nl,1);
% % 原始约束的凸松弛 
% Constraints = [Constraints, 0 <= Ploss]; 
% Constraints = [Constraints, 0 <= Qloss]; 
% % 线性化绝对值 
% Constraints = [Constraints, Pf_abs >= Pf, Pf_abs >= -Pf]; 
% Constraints = [Constraints, Qf_abs >= Qf, Qf_abs >= -Qf];
% Vm_2ft = (C_PVm(f,:) - C_PVm(t,:)) * Pbus + (C_QVm(f,:) - C_QVm(t,:)) * Qbus;%Vm_2(f(i)) - Vm_2(t(i))
% Vaft = (C_PVa(f,:) - C_PVa(t,:)) * Pbus + (C_QVa(f,:) - C_QVa(t,:)) * Qbus;%Va(f(i)) - Va(t(i))
% 
% % for i = 1:nl
% % Constraints=[Constraints; (Ploss(i) >= ((Vm_2ft(i))^2+(Vaft(i))^2) * real(Ysc(i))):'Loss_P'];
% % Constraints=[Constraints; (Qloss(i) >= ((Vm_2ft(i))^2+(Vaft(i))^2) * imag(Ysc(i))):'Loss_Q'];
% % end
% 
% Constraints=[Constraints; (Ploss == (Vmij_2+Vaij_2).* real(Ysc)):'Loss_P'];
% Constraints=[Constraints; (Qloss == (Vmij_2+Vaij_2).* imag(Ysc)):'Loss_Q'];
% % M1
% for i = 1:nb
%     Constraints=[Constraints; ((Vm_2ft(i))^2 <= Vmij_2(i))];
% end
% for i = 1:nb
%     Constraints=[Constraints; ((Vaft(i))^2 <= Vaij_2(i))];
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
Constraints=[Constraints, (Vmin.^2 <= Vm_2):'Voltage_low'];
Constraints=[Constraints, (Vm_2 <= Vmax.^2):'Voltage_high'];

%% 机组出力上下限不等式约束
% Constraints = [Constraints; real(Sgmin) <= Pg <= real(Sgmax)];                                  % complex power output of generator (p.u.)
% Constraints = [Constraints; imag(Sgmin) <= Qg <= imag(Sgmax)];

Constraints = [Constraints; (real(Sgmin) <= Pg):'Pg_low'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (Pg <= real(Sgmax)):'Pg_high'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (imag(Sgmin) <= Qg):'Qg_low'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (Qg <= imag(Sgmax)):'Qg_high'];   

%% 平衡节点电压相角为0
Constraints = [Constraints; Va(ref) == 0];                                  % reference bus angle (rad)

%% 线路容量不等式约束（配电网不需要此约束）
%t_poly=10;
if (nb <= nl)
%     Constraints = [Constraints,  Pf.^2 + Qf.^2 <= Fmax.^2];                  % branch flow limit at "from" end (p.u.)
%     Constraints = [Constraints,  Pt.^2 + Qt.^2 <= Fmax.^2];
%     for i = 1:size(Pf,1)
%     Constraints = [Constraints,  (Pf(i)^2 + Qf(i)^2 <= Fmax(i)^2):'Con_Sf'];                  % branch flow limit at "from" end (p.u.)
%     Constraints = [Constraints,  (Pt(i)^2 + Qt(i)^2 <= Fmax(i)^2):'Con_St'];
%     end
% 
     Constraints = [Constraints,  (-Fmax <= (GSF.GSF_PP_F * (Pbus) + GSF.GSF_PQ_F * (Qbus))):'Con_PF_low'];                  % branch flow limit at "from" end (p.u.)
     Constraints = [Constraints,  ((GSF.GSF_PP_F * (Pbus) + GSF.GSF_PQ_F * (Qbus)) <= Fmax):'Con_PF_high'];                  % branch flow limit at "from" end (p.u.)
     Constraints = [Constraints,  (-Fmax <= (GSF.GSF_QP_F * (Pbus) + GSF.GSF_QQ_F * (Qbus))):'Con_QF_low'];
     Constraints = [Constraints,  ((GSF.GSF_QP_F * (Pbus) + GSF.GSF_QQ_F * (Qbus)) <= Fmax):'Con_QF_high'];
      
     Constraints = [Constraints,  (Fmax >= (GSF.GSF_PP_T * (Pbus) + GSF.GSF_PQ_T * (Qbus))):'Con_PT_low'];                  % branch flow limit at "from" end (p.u.)
     Constraints = [Constraints,  ((GSF.GSF_PP_T * (Pbus) + GSF.GSF_PQ_T * (Qbus)) >= -Fmax):'Con_PT_high'];                  % branch flow limit at "from" end (p.u.)
     Constraints = [Constraints,  (Fmax >= (GSF.GSF_QP_T * (Pbus) + GSF.GSF_QQ_T * (Qbus))):'Con_QT_low'];
     Constraints = [Constraints,  ((GSF.GSF_QP_T * (Pbus) + GSF.GSF_QQ_T * (Qbus)) >= -Fmax):'Con_QT_high'];     

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
opt = sdpsettings('solver','+gurobi','verbose',1,'gurobi.NumericFocus',3);% ipopt gurobi cplex 
opt.gurobi.QCPDual = 1;
varargout=optimize(Constraints,obj,opt);
varargout.info

%% 检查约束违反
check(Constraints);

%% 输出结果
% -------- post-processing -------- %
%res.Ploss = sum(value(Ploss)) / baseMVA;
res.Vm = value(sqrt(Vm_2));
% res.Vmij_2 = value(Vmij_2);
% res.Vaij_2 = value(Vaij_2);
% res.Vm_2ft = value(Vm_2ft);
% res.Vaft = value(Vaft);
res.Va = value(Va)/pi*180;
res.Pbus = value(Pbus)* baseMVA;
res.Qbus = value(Qbus)* baseMVA;
res.Ploss = value(Ploss)* baseMVA;
res.Qloss = value(Qloss)* baseMVA;
res.Pf = value(Pf)* baseMVA;
res.Qf = value(Qf)* baseMVA;
res.Pt = value(Pt)* baseMVA;
res.Qt = value(Qt)* baseMVA;
res.Pg = value(Pg)* baseMVA;         % gen real power output (MW)
res.Qg = value(Qg)* baseMVA;
% res.PfPt = value(PfPt)* baseMVA;
% res.QfQt = value(QfQt)* baseMVA;
res.cost = value(obj);      % gen real power output cost ($/h)
% LMP：Energy
res.lmp_energy_P = dual(Constraints(1))/ baseMVA;   % lmp($/MWh)
res.lmp_energy_Q = dual(Constraints(2))/ baseMVA;   % lmp($/MWh)
LMP.PEnergy_Price = res.lmp_energy_P;
LMP.QEnergy_Price = res.lmp_energy_Q;
% LMP：Loss 
res.lmp_loss_P = dual(Constraints(1))/ baseMVA;   % lmp($/MWh)
res.lmp_loss_Q = dual(Constraints(2))/ baseMVA;   % lmp($/MWh)
LF_PP = ((C_PVm(f,:) - C_PVm(t,:)).^2 + (C_PVa(f,:) - C_PVa(t,:)).^2) .* real(Ysc);
LF_PQ = ((C_QVm(f,:) - C_QVm(t,:)).^2 + (C_QVa(f,:) - C_QVa(t,:)).^2) .* real(Ysc);
LF_P_Q = (2 * (C_PVm(f,:) - C_PVm(t,:)) .* (C_PVa(f,:) - C_PVa(t,:)) + 2 * (C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:))) .* real(Ysc);
LF_QP = ((C_PVm(f,:) - C_PVm(t,:)).^2 + (C_PVa(f,:) - C_PVa(t,:)).^2) .* imag(Ysc);
LF_QQ = ((C_QVm(f,:) - C_QVm(t,:)).^2 + (C_QVa(f,:) - C_QVa(t,:)).^2) .* imag(Ysc);
LF_Q_P = (2 * (C_PVm(f,:) - C_PVm(t,:)) .* (C_PVa(f,:) - C_PVa(t,:)) + 2 * (C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:))) .* imag(Ysc);
LMP.LF_P = value(2 * LF_PP * Pbus + LF_P_Q * Qbus);%LF_P_Q' * Pbus0 * Qbus0 =  P_offset
LMP.LF_Q = value(2 * LF_QQ * Qbus + LF_Q_P * Pbus);%LF_Q_P' * Qbus0 * Pbus0 =  Q_offset
% LF_P = value(2 * (C_PVm(f,:) - C_PVm(t,:)).^2 * Pbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + 2 * (C_PVa(f,:) - C_PVa(t,:)).^2 * Pbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* real(Ysc);
% LF_Q = value(2 * (C_QVm(f,:) - C_QVm(t,:)).^2 * Qbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + 2 * (C_QVa(f,:) - C_QVa(t,:)).^2 * Qbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* imag(Ysc);
LMP.PLoss_Price = -Cft' * LMP.LF_P .* res.lmp_energy_P;%
LMP.QLoss_Price = -Cft' * LMP.LF_Q .* res.lmp_energy_Q;%-Cf' * LF_Q .* res.lmp_energy_Q

% LMP：V
%res.lmp_V = dual(Constraints('V'));   % lmp($/MWh)
res.lmp_V_mul = dual(Constraints('Voltage_low'))/ baseMVA;   % lmp($/MWh)
res.lmp_V_muh = dual(Constraints('Voltage_high'))/ baseMVA;   % lmp($/MWh)
LMP.PVl_Price = C_PVm * res.lmp_V_mul;
LMP.PVh_Price = C_PVm * res.lmp_V_muh;
LMP.PV_Price = LMP.PVl_Price + LMP.PVh_Price;
LMP.QVl_Price = C_QVm * res.lmp_V_mul;
LMP.QVh_Price = C_QVm * res.lmp_V_muh;
LMP.QV_Price = LMP.QVl_Price + LMP.QVh_Price;

%res.lmp_mul_F = 0;res.lmp_muh_F = 0;res.lmp_mul_T = 0;res.lmp_muh_T = 0;
% for i = 1:t_poly
%     res.lmp_mul_F = res.lmp_mul_F + dual(Constraints(['Con_F_low', num2str(i)]));
%     res.lmp_muh_F = res.lmp_muh_F + dual(Constraints(['Con_F_high', num2str(i)]));
%     res.lmp_mul_T = res.lmp_mul_T + dual(Constraints(['Con_T_low', num2str(i)]));
%     res.lmp_muh_T = res.lmp_muh_T + dual(Constraints(['Con_T_high', num2str(i)]));
% end

% LMP：Con
% res.lmp_mu_SF = dual(Constraints('Con_Sf'))/ baseMVA;   % lmp($/MWh)
% res.lmp_mu_ST = dual(Constraints('Con_St'))/ baseMVA;   % lmp($/MWh)
% LMP.SfCon_Price = GSF_PP_F * res.lmp_mu_SF;
% LMP.StCon_Price = GSF_PP_F * res.lmp_mu_ST;
if (nb <= nl)
    res.lmp_mul_PF = dual(Constraints('Con_PF_low'))/ baseMVA;   % lmp($/MWh)
    res.lmp_muh_PF = dual(Constraints('Con_PF_high'))/ baseMVA;   % lmp($/MWh)
    res.lmp_mul_PT = dual(Constraints('Con_PT_low'))/ baseMVA;   % lmp($/MWh)
    res.lmp_muh_PT = dual(Constraints('Con_PT_high'))/ baseMVA;   % lmp($/MWh)
    res.lmp_mul_QF = dual(Constraints('Con_QF_low'))/ baseMVA;   % lmp($/MWh)
    res.lmp_muh_QF = dual(Constraints('Con_QF_high'))/ baseMVA;   % lmp($/MWh)
    res.lmp_mul_QT = dual(Constraints('Con_QT_low'))/ baseMVA;   % lmp($/MWh)
    res.lmp_muh_QT = dual(Constraints('Con_QT_high'))/ baseMVA;   % lmp($/MWh)
    LMP.PfCon_Price = GSF.GSF_PP_F' * res.lmp_mul_PF + GSF.GSF_PP_F' * res.lmp_muh_PF;
    LMP.PtCon_Price = GSF.GSF_PP_T' * res.lmp_mul_PT + GSF.GSF_PP_T' * res.lmp_muh_PT;
    LMP.QfCon_Price = GSF.GSF_QP_F' * res.lmp_mul_QF + GSF.GSF_QP_F' * res.lmp_muh_QF;
    LMP.QtCon_Price = GSF.GSF_QP_T' * res.lmp_mul_QT + GSF.GSF_QP_T' * res.lmp_muh_QT;
end

if (nb <= nl)
    LMP.Total_Active_Price = LMP.PEnergy_Price + LMP.PLoss_Price + (LMP.PfCon_Price + LMP.PtCon_Price) / 2;
    LMP.Total_Reactive_Price = LMP.QEnergy_Price + LMP.QLoss_Price + (LMP.QfCon_Price + LMP.QtCon_Price) / 2;
else
    LMP.Total_Active_Price = LMP.PEnergy_Price + LMP.PLoss_Price;
    LMP.Total_Reactive_Price = LMP.QEnergy_Price + LMP.QLoss_Price;
end
% LMP.Total_Active_Price = res.lmp_energy_P + LMP.PLoss_Price + (LMP.PVl_Price + LMP.PVh_Price) + (LMP.PfCon_Price + LMP.PtCon_Price) / 2;
% LMP.Total_Reactive_Price = res.lmp_energy_Q + LMP.QLoss_Price + (LMP.QVl_Price + LMP.QVh_Price) + (LMP.QfCon_Price + LMP.QtCon_Price) / 2;

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

if nargout == 1 || nargout == 2 || nargout == 3 || nargout == 4 || nargout == 5 || nargout == 6 || nargout == 7
    MVAbase = results;
    bus = success;
    gen = LMP;
    gencost = GSF;
    branch = C;
    f = res.Ploss;
    success = res.Qloss;
elseif nargout > 7
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