addpath('../../matlab')
addpath('./scs_matlab')

cd '../../matlab'; make_scs; cd '../examples/matlab';

randn('seed',0);rand('seed',0);

%% generate random cone problem (solution NOT necessariy unique):

%%% types of probs to solve:
gen_feasible = true;
gen_infeasible = true;
gen_unbounded = true;
%%% solvers to test:
run_indirect = false;
run_direct = true;
run_cvx = false; % won't work if ep or ed > 0
cvx_solver = 'sdpt3';
run_scs_matlab = true; % SCS implemented in MATLAB

% set cone sizes (ep = ed = 0 if you want to compare against cvx):
K = struct('f',50,'l',60,'q',[0;1;0;2;100;20],'s',[0;2;1;2;20],'ep',10,'ed',15,'p',[-0.3, 0.25, 0.75, -0.4]);

%K = struct('f',1,'l',0,'q',0,'s',0,'ep',0,'ed',0,'p',0);

density = 0.1; % A matrix density

m = getConeDims(K)
n = round(m/3);
%n = 2;
%%
params = struct('eps', 1e-3, 'normalize', 0, 'scale', 1,...
    'cg_rate',2, 'max_iters', 59, 'alpha', 1, 'rho_x',1,'use_indirect',false);
%% generate primal-dual feasible cone prob for linesearch:
randn('seed',0);rand('seed',0);
% Ax + s = b, s \in K, A'y + c = 0, y \in K*, s'y = 0

z = randn(m,1);
y = proj_dual_cone(z,K); % y = s - z;
s = y - z; %s = proj_cone(z,K);

%A = [0,0;0,0; 0,0; 1 ,1];
A=sprandn(m,n,density);
%x = [1;1];
x=randn(n,1);
c = -A'*y;
b = A*x + s;

data.A = A;
data.b = b;
data.c = c;

%cd '../../matlab'; write_scs_data(data,K,params,'randomConeFeasible'); cd '../examples/matlab';

params.line_search = 1;
[xi_m,yi_m,si_m,infoi_m,statsls] = scs_matlab(data,K,params);
params.line_search = 0;
[xd_m,yd_m,sd_m,infod_m,stats] = scs_matlab(data,K,params);
[c'*xi_m,(c'*xi_m - c'*x) / (c'*x), b'*yi_m, (b'*yi_m - b'*y) / (b'*y)]
[c'*xd_m,(c'*xd_m - c'*x) / (c'*x), b'*yd_m, (b'*yd_m - b'*y) / (b'*y)]