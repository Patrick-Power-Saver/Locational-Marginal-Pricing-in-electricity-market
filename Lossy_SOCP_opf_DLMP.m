function [MVAbase, bus, gen, gencost, branch, f, success, et] = Lossy_SOCP_opf_DLMP(mpc, varargin)
% Lossy Linear AC OPF
% 11-06-2023

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

Pd = bus(:, PD);              % bus load (MVA)
Qd = bus(:, QD);  
% P0 = real(makeSbus(baseMVA, bus, gen));
% Q0 = imag(makeSbus(baseMVA, bus, gen));
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
Cg = sparse(gbus, (1:ng)', 1, nb, ng);

f = branch(:, F_BUS);
t = branch(:, T_BUS);
Cf = sparse((1:nl)', f, 1, nl, nb);
Ct = sparse((1:nl)', t, 1, nl, nb);
i = [(1:nl)'; (1:nl)'];             % element k, j is -1 if branch k connects "to" bus j
Cft = sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb);
Ctf = sparse(i, [t; f], [ones(nl, 1); -ones(nl, 1)], nl, nb);
F = full(sparse(br,f,1,nl,nb));
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Ysc = 1 ./ (branch(:, BR_R) - 1j * branch(:, BR_X));
Gbus = real(Ybus);
Bbus = imag(Ybus);

GP = Gbus; %导纳矩阵的实部
BD = diag(sum(Bbus));
B_mod= Bbus-BD;% B'
BP = -B_mod; %BP=-B' 
BQ = -Bbus; %BQ=-B
GQ = -Gbus;                       %GQ approxiately equals -Gbus
Xp_dlpf = [BP GP];
Xq_dlpf = [GQ BQ];
CJ = full([Xp_dlpf;Xq_dlpf]);
CJ(ref,ref) = 0;
CJ = pinv(CJ);
C_PVm = CJ(nb+1:end,1:nb);
C_QVm = CJ(nb+1:end,nb+1:end);
C_PVa = CJ(1:nb,1:nb);
C_QVa = CJ(1:nb,nb+1:end);

Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA;
tap = ones(nl, 1);  
xfmr = find(branch(:, TAP));                    %% indices of transformers
tap(xfmr) = branch(xfmr, TAP);                  %% include transformer tap ratios
tap = tap .* exp(1j*pi/180 * branch(:, SHIFT)); %% add phase shifters

% Gsh = real(Ysh);%并联导纳为0
% Bsh = imag(Ysh);
% temp1=zeros(size(branch,1),4);
% for i=1:length(temp1)
%     temp1(i,1)=Gbus(branch(i,F_BUS),branch(i,T_BUS));
%     temp1(i,2)=Bbus(branch(i,F_BUS),branch(i,T_BUS));
%     temp1(i,3)=Gbus(branch(i,T_BUS),branch(i,F_BUS));
%     temp1(i,4)=Bbus(branch(i,T_BUS),branch(i,F_BUS));
% end

% g = 1 ./ branch(:, BR_R);
% b = 1 ./ branch(:, BR_X);
% Gf = spdiags(g, 0, nl, nl) * Cft;   % Gf * Va is the vector of real branch power (p.u.)
% Gt = spdiags(g, 0, nl, nl) * Ctf;
% Bf = spdiags(b, 0, nl, nl) * Cft;   % Bf * Va is the vector of real branch power (p.u.)
% Bt = spdiags(b, 0, nl, nl) * Ctf;
% Gbus_ft = Cft' * Gf;                   % Gbus * Va is the vector of nodal real power injection (p.u.)
% Gbus_tf = Ctf' * Gt;
% Bbus_ft = Cft' * Bf;                   % Bbus * Va is the vector of nodal real power injection (p.u.)
% Bbus_tf = Ctf' * Bt;

% 单位换算
% -------- unit conversion --------- %
Pd = Pd / baseMVA;
Qd = Qd / baseMVA;
% P0 = P0 / baseMVA;
% Q0 = Q0 / baseMVA;
Fmax = Fmax / baseMVA;          % p.u.
FLmax = sparse([f;t],[t,f],[Fmax;-Fmax],nb,nb);
Sgmin = Sgmin / baseMVA;        % p.u.
Sgmax = Sgmax / baseMVA;        % p.u.
C2 = C2 * (baseMVA^2);          % $/h
C1 = C1 * baseMVA;              % $/h

%% 定义变量
Vm_2 = sdpvar(nb, 1, 'full');
Vm = sdpvar(nb,1, 'full');
Vmij_2 = sdpvar(nl,1, 'full');
abs_Vmij_2 = sdpvar(nl,1, 'full');
Va = sdpvar(nb,1, 'full');
Vaij_2 = sdpvar(nl,1, 'full');
abs_Vaij_2 = sdpvar(nl,1, 'full');
Pbus = sdpvar(nb,1, 'full');
Qbus = sdpvar(nb,1, 'full');
Pf = sdpvar(nl,1, 'full');
Qf = sdpvar(nl,1, 'full');
Pt = sdpvar(nl,1, 'full');
Qt = sdpvar(nl,1, 'full');
% PL = sparse([f;t],[t;f],[Pf;Pt],nb,nb);
% QL = sparse([f;t],[t;f],[Qf;Qt],nb,nb);
Pg = sdpvar(ng, 1, 'full');
Qg = sdpvar(ng, 1, 'full');
%Loss= sdpvar(nl,1, 'full');
Ploss = sdpvar(nl,1, 'full');
Qloss = sdpvar(nl,1, 'full');
Pequ = sdpvar(nb, 1, 'full');
Qequ = sdpvar(nb, 1, 'full');

% 为节点电压赋初值
% -------- assign initial value -------- %
% Matpower solution
% mpopt = mpoption( 'out.all', 0);
% result = runopf(mpc, mpopt);
% assign(Vm, result.bus(:,VM));
% assign(Va, result.bus(:,VA)*pi/180);
% Sbus = makeSbus(result.baseMVA, result.bus, result.gen);
% assign(Pbus, real(Sbus));
% assign(Qbus, imag(Sbus));
% assign(Pg, result.gen(:,PG));
% assign(Qg, result.gen(:,QG));

%% 目标函数
obj = C2' * (Pg.^2) + C1' * Pg + ones(ng, 1)' * C0;         % linear cost function ($/h)

%% 约束条件
Constraints = [];

%% 节点功率平衡等式约束
Constraints = [Constraints; (Pbus + Pd + Pequ - Cg * Pg == 0):'Balance_P'];
Constraints = [Constraints; (Qbus + Qd + Qequ - Cg * Qg == 0):'Balance_Q'];

%% B 节点注入与节点电压的线性等式约束（B,C二选一）
Constraints=[Constraints;Pbus == Xp_dlpf * [Va;Vm_2]];
Constraints=[Constraints;Qbus == Xq_dlpf * [Va;Vm_2]];
% Constraints = [Constraints; (Vm_2 == C_PVm * Pbus + C_QVm * Qbus):'V'];
% Constraints = [Constraints; Va == C_PVa * Pbus + C_QVa * Qbus];

%******************配网******************
% Constraints = [Constraints; Ploss == ((Vmij_2) .* real(Ysc))];
% Constraints = [Constraints; Qloss == ((Vmij_2) .* imag(Ysc))];
% Constraints = [Constraints;(Vm_2(f)-Vm_2(t)).^2 <= Vmij_2];
Constraints=[Constraints; (Ploss == (Vmij_2+Vaij_2).* real(Ysc)):'Loss_P'];
Constraints=[Constraints; (Qloss == (Vmij_2+Vaij_2).* imag(Ysc)):'Loss_Q'];
% M1
for i = 1:size(f,1)
    Constraints=[Constraints; (Vm_2(f(i)) - Vm_2(t(i)))^2 <= Vmij_2(i)];
end
for i = 1:size(f,1)
    Constraints=[Constraints; (Va(f(i)) - Va(t(i)))^2 <= Vaij_2(i)];
end
%Constraints = [Constraints;abs_Vmij_2 >= Vmij_2];
% Constraints = [Constraints;abs_Vmij_2 >= -Vmij_2];
% Constraints = [Constraints;(Va(f)-Va(t)).^2 <= Vaij_2];
% Constraints = [Constraints;abs_Vaij_2 >= Vaij_2];
% Constraints = [Constraints;abs_Vaij_2 >= -Vaij_2];
% 
% Constraints = [Constraints; Pequ == 0.5 * Cf' * (Vmij_2 + Ploss)];
% Constraints = [Constraints; Qequ == 0.5 * Cf' * (Vmij_2 + Qloss)];
Constraints = [Constraints; Pequ == 0.5 * (Cf'+Ct') * (Ploss)];
Constraints = [Constraints; Qequ == 0.5 * (Cf'+Ct') * (Qloss)];
% Constraints=[Constraints;Pequ == 0.5 * Cf' * Ploss];
% Constraints=[Constraints;Qequ == 0.5 * Cf' * Qloss];
% Constraints=[Constraints; Ploss >= 0];
% Constraints=[Constraints; Qloss >= 0];

% 节点注入与节点电压的线性等式约束（不需要）
% Constraints=[Constraints;Pm1 - Pbus([pq;pv]) + GP([pq;pv],ref)*bus(ref,VM) + GP([pq;pv],pv)*Vm(pv) == 0];  
% Constraints=[Constraints;Qm1 - Qbus(pq) + BQ(pq,ref)*bus(ref,VM) + BQ(pq,pv)*Vm(pv) == 0];
% Constraints=[Constraints;Pm2 - Pm1 + t1*Qm1 == 0];
% Constraints=[Constraints;Qm2 - Qm1 + t2*Pm1 == 0];
% Constraints=[Constraints;Vm(pq) == B22m\Qm2];
% Constraints=[Constraints;Va([pq;pv]) == B11m\Pm2];
% % Constraints=[Constraints;Vm([ref]) == result.bus([ref],VM)];
% Constraints=[Constraints;Va([ref]) == 0];

%% 节点注入与支路功率的线性等式约束
% for i = 1:nb %节点数量
%     Corrlbranchij=SearchNodeConnection(f,t,no(i,1)); %获取每个节点对应的支路信息
%     SumCorrlBranchACP(i,1) = sum(PL(bus(i,1),Corrlbranchij(:,2)));
%     SumCorrlBranchACQ(i,1) = sum(QL(bus(i,1),Corrlbranchij(:,2)));
% %     SumCorrlPLoss(i,1)  = sum(PLoss(bus(i,1),Corrlbranchij(:,2)));
% %     SumCorrlQLoss(i,1)  = sum(QLoss(bus(i,1),Corrlbranchij(:,2)));
%     SumGii(1,no(i,1))=sum(Gbus(no(i,1),:));
%     SumBii(1,no(i,1))=-sum(Bbus(no(i,1),:));
% end
% for i = 1:nb %KCL_
%     Constraints=[Constraints,Pbus(no(i,1)) == SumCorrlBranchACP(no(i,1))+Vm(no(i,1))*SumGii(no(i,1)) ];%
%     Constraints=[Constraints,Qbus(no(i,1))== SumCorrlBranchACQ(no(i,1))+Vm(no(i,1))*SumBii(no(i,1)) ];%
% end

% Constraints=[Constraints,Pbus == Cft' * Pf];
% Constraints=[Constraints,Qbus == Cft' * Qf];
% Constraints=[Constraints,Pbus == Ctf' * Pt];
% Constraints=[Constraints,Qbus == Ctf' * Qt];

% Constraints=[Constraints,Pbus - 0.5 * F' * Ploss / baseMVA  == Cft' * Pf + sum(Gbus,2).*Vm];
% Constraints=[Constraints,Qbus - 0.5 * F' * Qloss / baseMVA  == Cft' * Qf - sum(Bbus,2).*Vm];
% Constraints=[Constraints,Pbus - 0.5 * F' * Ploss/ baseMVA  == Ctf' * Pt + sum(Gbus,2).*Vm];
% Constraints=[Constraints,Qbus - 0.5 * F' * Qloss / baseMVA == Ctf' * Qt - sum(Bbus,2).*Vm];

%% 支路功率与节点电压的等式约束
% (Vi - Vj)
% Constraints=[Constraints; Pf == ((Vm(branch(br,F_BUS)) - Vm(branch(br,T_BUS))).*temp1(br,1) ...
%     -(Va(branch(br,F_BUS))-Va(branch(br,T_BUS))).*temp1(br,2))]; % the sitas in equation (27) is 0,so neglected.
% Constraints=[Constraints; Qf == -((Vm(branch(br,F_BUS)) - Vm(branch(br,T_BUS))).*temp1(br,2) ...
%     -(Va(branch(br,F_BUS))-Va(branch(br,T_BUS))).*temp1(br,1))];  % the sitas in equation (27) is 0,so neglected.
% Constraints=[Constraints; Pt == ((Vm(branch(br,T_BUS)) - Vm(branch(br,F_BUS))).*temp1(br,3) ...
%     -(Va(branch(br,T_BUS))-Va(branch(br,F_BUS))).*temp1(br,4))]; % the sitas in equation (27) is 0,so neglected.
% Constraints=[Constraints; Qt == -((Vm(branch(br,T_BUS)) - Vm(branch(br,F_BUS))).*temp1(br,4) ...
%     -(Va(branch(br,T_BUS))-Va(branch(br,F_BUS))).*temp1(br,3))];  % the sitas in equation (27) is 0,so neglected.
% 
% Constraints=[Constraints; Pf == (Vm(f)-Vm(t)).* real(Ysc) - (Va(f)-Va(t)).*imag(Ysc)];
% Constraints=[Constraints; Qf == -(Vm(f)-Vm(t)).* imag(Ysc) - (Va(f)-Va(t)).*real(Ysc)];
% Constraints=[Constraints; Pt == (Vm(t)-Vm(f)).* real(Ysc) - (Va(t)-Va(f)).*imag(Ysc)];
% Constraints=[Constraints; Qt == -(Vm(t)-Vm(f)).* imag(Ysc) - (Va(t)-Va(f)).*real(Ysc)];

% (Vi.^2 - Vj.^2)/2
% Constraints=[Constraints; Pf == ((Vm_2(branch(br,F_BUS)) - Vm_2(branch(br,T_BUS))).*temp1(br,1) / 2 ...
%     -(Va(branch(br,F_BUS))-Va(branch(br,T_BUS))).*temp1(br,2))]; % the sitas in equation (27) is 0,so neglected.
% Constraints=[Constraints; Qf == -((Vm_2(branch(br,F_BUS)) - Vm_2(branch(br,T_BUS))).*temp1(br,2) / 2 ...
%     -(Va(branch(br,F_BUS))-Va(branch(br,T_BUS))).*temp1(br,1))];  % the sitas in equation (27) is 0,so neglected.
% Constraints=[Constraints; Pt == ((Vm_2(branch(br,T_BUS)) - Vm_2(branch(br,F_BUS))).*temp1(br,3) / 2 ...
%     -(Va(branch(br,T_BUS))-Va(branch(br,F_BUS))).*temp1(br,4))]; % the sitas in equation (27) is 0,so neglected.
% Constraints=[Constraints; Qt == -((Vm_2(branch(br,T_BUS)) - Vm_2(branch(br,F_BUS))).*temp1(br,4) / 2 ...
%     -(Va(branch(br,T_BUS))-Va(branch(br,F_BUS))).*temp1(br,3))];  % the sitas in equation (27) is 0,so neglected.
GSF_PP_F = (C_PVm(f,:) - C_PVm(t,:)) .* real(Ysc) - (C_PVa(f,:) - C_PVa(t,:)) .* imag(Ysc);
GSF_PQ_F = (C_QVm(f,:) - C_QVm(t,:)) .* real(Ysc) - (C_QVa(f,:) - C_QVa(t,:)) .* imag(Ysc);
GSF_QP_F = -(C_PVa(f,:) - C_PVa(t,:)) .* real(Ysc) - (C_PVm(f,:) - C_PVm(t,:)) .* imag(Ysc); 
GSF_QQ_F = -(C_QVa(f,:) - C_QVa(t,:)) .* real(Ysc) - (C_QVm(f,:) - C_QVm(t,:)) .* imag(Ysc);
GSF_PP_T = -GSF_PP_F;
GSF_PQ_T = -GSF_PQ_F;
GSF_QP_T = -GSF_QP_F;
GSF_QQ_T = -GSF_QQ_F;
Constraints=[Constraints; Pf == (GSF_PP_F * Pbus + GSF_PQ_F * Qbus) + 0.5 * Ploss];
Constraints=[Constraints; Qf == (GSF_QP_F * Pbus + GSF_QQ_F * Qbus) + 0.5 * Qloss];
Constraints=[Constraints; Pt == (GSF_PP_T * Pbus + GSF_PQ_T * Qbus) + 0.5 * Ploss];
Constraints=[Constraints; Qt == (GSF_QP_T * Pbus + GSF_QQ_T * Qbus) + 0.5 * Qloss];

% Constraints=[Constraints; Pf  == (Vm_2(f)-Vm_2(t)).* real(Ysc) / 2 - (Va(f)-Va(t)).*imag(Ysc) + 0.5 * Ploss];
% Constraints=[Constraints; Qf  == -(Vm_2(f)-Vm_2(t)).* imag(Ysc)/ 2 - (Va(f)-Va(t)).*real(Ysc) + 0.5 * Qloss];
% Constraints=[Constraints; Pt  == (Vm_2(t)-Vm_2(f)).* real(Ysc)/ 2 - (Va(t)-Va(f)).*imag(Ysc) + 0.5 * Ploss];
% Constraints=[Constraints; Qt  == -(Vm_2(t)-Vm_2(f)).* imag(Ysc)/ 2 - (Va(t)-Va(f)).*real(Ysc) + 0.5 * Qloss];

% % 使用McCormick包络法进行凸松弛
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

% Constraints=[Constraints; Pf + Pt == Ploss];
% Constraints=[Constraints; Qf + Qt == Qloss];

% Constraints=[Constraints; Gf * Vm - Bf * Va == Pf];
% Constraints=[Constraints; Gt * Vm - Bt * Va == Pt];
% Constraints=[Constraints; -Bf * Vm - Gf * Va == Qf];
% Constraints=[Constraints; -Bt * Vm - Gt * Va == Qt];

%% 电压上下限不等式约束
%Constraints = [Constraints; Vmin <= Vm <= Vmax];         % voltage magnitude (p.u.)
Constraints=[Constraints, (Vmin.^2 <= Vm_2):'Voltage_low'];
Constraints=[Constraints, (Vm_2 <= Vmax.^2):'Voltage_high'];
%Constraints = [Constraints; -pi/6 <= Va <= pi/6];
% Constraints = [Constraints;Voltage.min<=Vol(pq,1)<=Voltage.max];

%% 机组出力上下限不等式约束
Constraints = [Constraints; (real(Sgmin) <= Pg):'Pg_low'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (Pg <= real(Sgmax)):'Pg_high'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (imag(Sgmin) <= Qg):'Qg_low'];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; (Qg <= imag(Sgmax)):'Qg_high'];   

%% 平衡节点电压相角为0
Constraints = [Constraints; Va(ref) == 0];                                  % reference bus angle (rad)

%% 线路容量不等式约束（配电网不需要此约束）
if (nb <= nl)
%     Constraints = [Constraints,  Pf.^2 + Qf.^2 <= Fmax.^2];                  % branch flow limit at "from" end (p.u.)
%     Constraints = [Constraints,  Pt.^2 + Qt.^2 <= Fmax.^2];
    %Constraints = [Constraints,  PL.^2 + QL.^2 <= FLmax.^2];
    Constraints = [Constraints,  (-Fmax <= (GSF_PP_F * Pbus + GSF_PQ_F * Qbus)):'Con_PF_low'];                  % branch flow limit at "from" end (p.u.)
    Constraints = [Constraints,  ((GSF_PP_F * Pbus + GSF_PQ_F * Qbus) <= Fmax):'Con_PF_high'];                  % branch flow limit at "from" end (p.u.)
    Constraints = [Constraints,  (-Fmax <= (GSF_QP_F * Pbus + GSF_QQ_F * Qbus)):'Con_QF_low'];
    Constraints = [Constraints,  ((GSF_QP_F * Pbus + GSF_QQ_F * Qbus) <= Fmax):'Con_QF_high'];
    
    Constraints = [Constraints,  (Fmax >= (GSF_PP_T * Pbus + GSF_PQ_T * Qbus)):'Con_PT_low'];                  % branch flow limit at "from" end (p.u.)
    Constraints = [Constraints,  ((GSF_PP_T * Pbus + GSF_PQ_T * Qbus) >= -Fmax):'Con_PT_high'];                  % branch flow limit at "from" end (p.u.)
    Constraints = [Constraints,  (Fmax >= (GSF_QP_T * Pbus + GSF_QQ_T * Qbus)):'Con_QT_low'];
    Constraints = [Constraints,  ((GSF_QP_T * Pbus + GSF_QQ_T * Qbus) >= -Fmax):'Con_QT_high'];
else
end

% %% 交流线路潮流上下限不等式约束 
% Npart=30; %将上、下半圆截为20份线段
% alpha=pi/30; %6度
% beta=(pi-2*alpha)/Npart;
% M=Npart;  %上半圆份数
% N=Npart;  %下半圆份数
% Smax=1;
% %**********上半圆**********%
% KAU=zeros(Npart,1);KBU=zeros(Npart,1);XPAU=zeros(Npart,1);
% YPAU=zeros(Npart,1);XPBU=zeros(Npart,1);YPBU=zeros(Npart,1);
% %**********下半圆**********%
% KAD=zeros(Npart,1); KBD=zeros(Npart,1);XPAD=zeros(Npart,1);
% YPAD=zeros(Npart,1);XPBD=zeros(Npart,1);YPBD=zeros(Npart,1);
% 
% %**********上半圆**********%  
% for i=1:M
%     KAU(i)=tan(alpha+(i-1)*beta);
%     KBU(i)=tan(alpha+i*beta);
%     XPAU(i)=1/( sqrt( 1+( 1/(KAU(i))^2 ) )*KAU(i) ) * Smax;
%     YPAU(i)=1/( sqrt( 1+( 1/(KAU(i))^2 ) ) ) * Smax;
%     XPBU(i)=1/( sqrt( 1+( 1/(KBU(i))^2 ) )*KBU(i) ) * Smax;
%     YPBU(i)=1/( sqrt( 1+( 1/(KBU(i))^2 ) ) ) * Smax;
%     Constraints=[Constraints, ( ( YPBU(i)-YPAU(i) )/( XPBU(i)-XPAU(i) )*(QL(f,t) - XPAU(i) )+YPAU(i)-PL(f,t) )>=0,... %交流支路
%     ];  
% end
% 
% %**********下半圆**********%
% 
% for i=1:N
%     KAD(i)=tan(-alpha-(i-1)*beta);
%     KBD(i)=tan(-alpha-i*beta);
%     XPAD(i)=-1/( sqrt( 1+( 1/(KAD(i))^2 ) )*KAD(i) ) * Smax;
%     YPAD(i)=-1/( sqrt( 1+( 1/(KAD(i))^2 ) ) ) * Smax;
%     XPBD(i)=-1/( sqrt( 1+( 1/(KBD(i))^2 ) )*KBD(i) ) * Smax;
%     YPBD(i)=-1/( sqrt( 1+( 1/(KBD(i))^2 ) ) ) * Smax;
%     Constraints=[Constraints, ( ( YPBD(i)-YPAD(i) )/( XPBD(i)-XPAD(i) )*(QL(f,t) - XPAD(i) )+YPAD(i)-PL(f,t) )<=0,...%交流支路
%     ]; %交流支路
% end

%% 求解问题
opt = sdpsettings('solver','gurobi','verbose',0);%gurobi cplex
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
% res.PfPt = value(PfPt)* baseMVA;
% res.QfQt = value(QfQt)* baseMVA;
res.Va = value(Va)/pi*180;
res.Pbus = value(Pbus)* baseMVA;
res.Qbus = value(Qbus)* baseMVA;
res.Ploss = value(Ploss);
res.Qloss = value(Qloss);
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
% LMP：Loss
% res.lmp_loss_P = -0.5 * (Cf + Ct) * res.lmp_energy_P;   %dual(Constraints('Loss_P'))/ baseMVA -0.5* Cf*res.lmp_energy_P lmp($/MWh)
% res.lmp_loss_Q = -0.5 * (Cf + Ct) * res.lmp_energy_Q;   % dual(Constraints('Loss_Q'))/ baseMVA -0.5 * Cf * res.lmp_energy_Q lmp($/MWh)
% LF_P = value((C_PVm(f,:) - C_PVm(t,:)).^2 * Pbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + (C_PVa(f,:) - C_PVa(t,:)).^2 * Pbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* real(Ysc);
% LF_Q = value((C_QVm(f,:) - C_QVm(t,:)).^2 * Qbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + (C_QVa(f,:) - C_QVa(t,:)).^2 * Qbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* imag(Ysc);
% LMP.PLoss_Price = -(Cf + Ct)' * LF_P .* res.lmp_energy_P;%-Cf' * LF_P .* res.lmp_energy_P (-Cf'* LF_P).*(-2*pinv(Cf)*res.lmp_loss_P)
% LMP.QLoss_Price = -(Cf + Ct)' * LF_Q .* res.lmp_energy_Q;%-Cf' * LF_Q .* res.lmp_energy_Q (-Cf'* LF_Q).*(-2*pinv(Cf)*res.lmp_loss_Q)

% res.lmp_loss_P = -0.5 * (Cf + Ct) * res.lmp_energy_P;   %dual(Constraints('Loss_P'))/ baseMVA -0.5* Cf*res.lmp_energy_P lmp($/MWh)
% res.lmp_loss_Q = -0.5 * (Cf + Ct) * res.lmp_energy_Q;   % dual(Constraints('Loss_Q'))/ baseMVA -0.5 * Cf * res.lmp_energy_Q lmp($/MWh)
% LMP.LF_P = value((C_PVm(f,:) - C_PVm(t,:)).^2 * Pbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + (C_PVa(f,:) - C_PVa(t,:)).^2 * Pbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* real(Ysc);
% LMP.LF_Q = value((C_QVm(f,:) - C_QVm(t,:)).^2 * Qbus + 2 * ((C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) * Qbus) ...
%     + (C_QVa(f,:) - C_QVa(t,:)).^2 * Qbus + 2 * ((C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:)) * Qbus)) .* imag(Ysc);
% LMP.PLoss_Price = (-(Cf+Ct)'* LMP.LF_P).*(-2*pinv(Cf + Ct + lambda * eye(size(Cf)))*res.lmp_loss_P);%-Cf' * LF_P .* res.lmp_energy_P
% LMP.QLoss_Price = (-(Cf+Ct)'* LMP.LF_Q).*(-2*pinv(Cf + Ct + lambda * eye(size(Cf)))*res.lmp_loss_Q);%-Cf' * LF_Q .* res.lmp_energy_Q

lambda = 1e-10; % 正则化参数 
LMP.SP_loss_P = 2 * pinv(Cf + Ct + lambda * eye(size(Cf))) * dual(Constraints('Loss_P'))/ baseMVA;   % lmp($/MWh)
LMP.SP_loss_Q = 2 * pinv(Cf + Ct + lambda * eye(size(Cf))) * dual(Constraints('Loss_Q'))/ baseMVA;   % lmp($/MWh)
LF_PP = ((C_PVm(f,:) - C_PVm(t,:)).^2 + (C_PVa(f,:) - C_PVa(t,:)).^2) .* real(Ysc);
LF_PQ = ((C_QVm(f,:) - C_QVm(t,:)).^2 + (C_QVa(f,:) - C_QVa(t,:)).^2) .* real(Ysc);
LF_P_Q = (2 * (C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) + 2 * (C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:))) .* real(Ysc);
LF_QP = ((C_PVm(f,:) - C_PVm(t,:)).^2 + (C_PVa(f,:) - C_PVa(t,:)).^2) .* imag(Ysc);
LF_QQ = ((C_QVm(f,:) - C_QVm(t,:)).^2 + (C_QVa(f,:) - C_QVa(t,:)).^2) .* imag(Ysc);
LF_Q_P = (2 * (C_PVm(f,:) - C_PVm(t,:)) .* (C_QVm(f,:) - C_QVm(t,:)) + 2 * (C_PVa(f,:) - C_PVa(t,:)) .* (C_QVa(f,:) - C_QVa(t,:))) .* imag(Ysc);
LMP.LF_P = value(2 * LF_PP * Pbus + LF_P_Q * Qbus);%LF_P_Q' * Pbus0 * Qbus0 =  P_offset
LMP.LF_Q = value(2 * LF_QQ * Qbus + LF_Q_P * Pbus);%LF_Q_P' * Qbus0 * Pbus0 =  Q_offset
LMP.PLoss_Price = (-(Cft)'* LMP.LF_P).*(-2 * LMP.SP_loss_P);%-Cf' * LF_P .* res.lmp_energy_P
LMP.QLoss_Price = (-(Cft)'* LMP.LF_Q).*(-2 * LMP.SP_loss_Q);%-Cf' * LF_Q .* res.lmp_energy_Q

% p_branch = (res.Pf - res.Pt)/2/ baseMVA;
% r = (branch(:, BR_R));
% T = makePTDF(mpc);
% LF = 2*(r.*p_branch)' * T;
% LF = LF';
% LD = 1-LF;
% (res.lmp_loss_P' * T)' .*-LF

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
    LMP.PfCon_Price = GSF_PP_F' * res.lmp_mul_PF + GSF_PP_F' * res.lmp_muh_PF;
    LMP.PtCon_Price = GSF_PP_T' * res.lmp_mul_PT + GSF_PP_T' * res.lmp_muh_PT;
    LMP.QfCon_Price = GSF_QP_F' * res.lmp_mul_QF + GSF_QP_F' * res.lmp_muh_QF;
    LMP.QtCon_Price = GSF_QP_T' * res.lmp_mul_QT + GSF_QP_T' * res.lmp_muh_QT;
end

if (nb <= nl)
    LMP.Total_Active_Price = LMP.PEnergy_Price + (LMP.PfCon_Price + LMP.PtCon_Price) / 2;
    LMP.Total_Reactive_Price = LMP.QEnergy_Price + (LMP.QfCon_Price + LMP.QtCon_Price) / 2;
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

if nargout == 1 || nargout == 2 || nargout == 3 || nargout == 4 || nargout == 5
    MVAbase = results;
    bus = success;
    gen = LMP;
    gencost = res.Ploss;
    branch = res.Qloss;
elseif nargout > 5
    [MVAbase, bus, gen, branch, f, et] = ...
        deal(results.baseMVA, results.bus, results.gen, results.branch, results.f, results.et);
    % else  %% don't define MVAbase, so it doesn't print anything
end
disp('done');
end