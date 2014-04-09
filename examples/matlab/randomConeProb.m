function randomConeProb()
clear all;close all;
addpath('../../matlab')
cd '../../matlab'; make_scs; cd '../examples/matlab';
%% generate random cone problem (solution NOT necessariy unique):

%%% types of probs to solve:
gen_feasible = true;
gen_infeasible = true;
gen_unbounded = true;
%%% solvers to test:
run_direct = true;
run_indirect = true;
run_cvx = true; % won't work if ep or ed > 0
cvx_solver = 'sdpt3';

% set cone sizes (ep = ed = 0 if you want to compare against cvx):
K = struct('f',100,'l',150,'q',[2;3;4;5;6;7;8;9;10;5;6;100;0;1],'s',[5;5;0;1;50],'ep',5,'ed',5)

density = 0.1; % A matrix density

m = getConeDims(K);
n = round(m/3);
params = struct('EPS',1e-4, 'NORMALIZE',1,'SCALE',5,'CG_RATE',1.5);

%% generate primal-dual feasible cone prob:
% Ax + s = b, s \in K, A'y + c = 0, y \in K*, s'y = 0
if (gen_feasible)
    z = randn(m,1);
    z = symmetrizeSDP(z,K); % for SD cones
    y = proj_dual_cone(z,K); % y = s - z;
    s = y - z; %s = proj_cone(z,K);
    
    
    A = sprandn(m,n,density);
    x = randn(n,1);
    c = -A'*y;
    b = A*x + s;
    nnz(A)
    
    data.A = A;
    data.b = b;
    data.c = c;
    
    %cd '../../matlab'; write_scs_data(data,K,params,'randomConeFeasible'); cd '../examples/matlab';

    %indirect
    if (run_indirect)
        [xi,yi,si,infoi] = scs_indirect(data,K,params);
        c'*x
        (c'*xi - c'*x) / (c'*x)
        b'*y
        (b'*yi - b'*y) / (b'*y)
    end
    if (run_direct)
        % direct:
        [xd,yd,sd,infod] = scs_direct(data,K,params);
        c'*x
        (c'*xd - c'*x) / (c'*x)
        b'*y
        (b'*yd - b'*y) / (b'*y)
    end
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
end

%% generate infeasible (NOT SPARSE,SLOW AS A RESULT)
% A'y = 0, y \in K*, b'*y = -1
if (gen_infeasible)
    z = randn(m,1);
    z = symmetrizeSDP(z,K); % for SD cones
    y = proj_dual_cone(z,K); % y = s - z;
    A = randn(m,n);
    
    A = A - ((A'*y)*y'/norm(y)^2)'; % dense...
    
    b = randn(m,1);
    b = -b / (b'*y);
    
    data.A = sparse(A);
    data.b = b;
    data.c = randn(n,1);
    
    params.SCALE = 0.5;
    %indirect
    if(run_indirect) [xi,yi,si,infoi] = scs_indirect(data,K,params); end
    % direct:
    if (run_direct) [xd,yd,sd,infod] = scs_direct(data,K,params); end
    
    % cvx:
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
    
end
%% generate unbounded (NOT SPARSE,SLOW AS A RESULT)
% Ax + s = 0, s \in K, c'*x = -1
if(gen_unbounded)
    z = randn(m,1);
    z = symmetrizeSDP(z,K); % for SD cones
    s = proj_cone(z,K);
    A = randn(m,n);
    x = randn(n,1);
    A = A - (s + A*x)*x'/(norm(x)^2); % dense...
    c = randn(n,1);
    c = - c / (c'*x);
    
    data.A = sparse(A);
    data.b = randn(m,1);
    data.c = c;
    
    params.SCALE = 0.5;
    %indirect
    if(run_indirect) [xi,yi,si,infoi] = scs_indirect(data,K,params); end
    % direct:
    if (run_direct) [xd,yd,sd,infod] = scs_direct(data,K,params); end
    
    if (run_cvx) [xc,yc,sc] = solveConeCvx(data,K,cvx_solver); end
end

end

function [x,y,s]=solveConeCvx(data,K,cs)
if (K.ep>0 || K.ed>0)
    x = nan;
    y = nan;
    s = nan;
    disp('cvx cannot solve EXPs');
    return;
end
n = length(data.c);
m = length(data.b);
% can NOT solve EXPs
cvx_begin
cvx_solver(cs)
variables x(n) s(m)
dual variable y
minimize(data.c'*x)
y:data.A*x + s == data.b
s(1:K.f)==0
l = K.f;
s(l+1:l+K.l) >= 0
l = l + K.l;
for i=1:length(K.q)
    s(l+1) >= norm(s(l+2:l + K.q(i)));
    l = l + K.q(i);
end
for i=1:length(K.s)
    reshape(s(l+1:l + (K.s(i))^2),K.s(i),K.s(i)) == semidefinite(K.s(i));
    l = l + (K.s(i))^2;
end
cvx_end
end

function l = getConeDims(K)
l = K.f + K.l;
for i=1:length(K.q)
    l = l + K.q(i);
end
for i=1:length(K.s)
    l = l + (K.s(i))^2;
end
l = l + K.ep*3;
l = l + K.ed*3;
end

function z = proj_dual_cone(z,c)
z = z + proj_cone(-z,c);
end

function z = proj_cone(z,c)
free_len = c.f;
lp_len = c.l;
k_soc = length(c.q);
q = c.q;
s = c.s;
ssize = length(c.s);
% free/zero cone
z(1:free_len) = 0;
% lp cone
z(free_len+1:lp_len+free_len) = pos(z(free_len+1:lp_len+free_len));
% SOCs
idx=lp_len+free_len;
for i=1:k_soc
    z(idx+1:idx+q(i)) = proj_soc(z(idx+1:idx+q(i)));
    idx=idx+q(i);
end
%SDCs
for i=1:ssize
    z(idx+1:idx+s(i)^2) = proj_sdp(z(idx+1:idx+s(i)^2),s(i));
    idx=idx+s(i)^2;
end
%Exp primal
for i=1:c.ep
    z(idx+1:idx+3) = project_exp_bisection(z(idx+1:idx+3));
    idx = idx+3;
end
%Exp dual
for i=1:c.ed
    z(idx+1:idx+3) = z(idx+1:idx+3) + project_exp_bisection(-z(idx+1:idx+3));
    idx = idx+3;
end

end

function z = proj_soc(tt)
if isempty(tt)
    z=[];
    return;
elseif length(tt)==1
    z = pos(tt);
    return;
end
v1=tt(1);v2=tt(2:end);
if norm(v2)<=-v1
    v2=zeros(length(v2),1);
    v1=0;
elseif norm(v2) > abs(v1)
    v2=0.5*(1+v1/norm(v2))*v2;
    v1=norm(v2);
end
z=[v1;v2];
end

function z = proj_sdp(z,n)
if isempty(z)
    z=[];
    return;
elseif length(z)==1
    z = pos(z);
    return;
end

z = reshape(z,n,n);
zs=(z+z')/2;

%ii
[V,S] = eig(zs);
S = diag(S);
num_pos = sum(S>0);
num_neg = sum(S<0);
if (num_pos < num_neg)
    positive = true;
    idx = find(S>0);
    V = V(:,idx);
    S = S(idx);
else
    positive = false;
    idx = find(S<0);
    V = V(:,idx);
    S = S(idx);
end

if (positive)
    T = S;
    T(T<0) = 0;
    z = V*diag(T)*V';
else
    T = S;
    T(T>0) = 0;
    z = zs - V*diag(T)*V';
end
z = z(:);
end

function x = project_exp_bisection(v)
r = v(1); s = v(2); t = v(3);
% v in cl(Kexp)
if( (s.*exp(r./s) <= t && s > 0) || (r <= 0 && s == 0 && t >= 0) );
    x = v;
    return
end
x = zeros(3,1);
% -v in Kexp^*
if ( (-r < 0 && r.*exp(s./r) <= -exp(1).*t) || (-r == 0 && -s >= 0 && -t >= 0) );
    return
end

% special case with analytical solution
if(r < 0 && s < 0);
    x = v;
    x(2) = 0;
    x(3) = max(v(3),0);
    return
end

x = v;
[ub,lb] = getRhoUb(v);
for iter=1:1e2;
    rho = (ub + lb)/2;
    [g,x] = calcGrad(v,rho,x);
    if (g > 0)
        lb = rho;
    else
        ub = rho;
    end
    if (ub - lb < 1e-6)
        break
    end
end
end

function [ub,lb] = getRhoUb(v)
lb = 0;
rho = 2^(-3);
[g,z] = calcGrad(v,rho,v);
while (g>0)
    lb = rho;
    rho = rho*2;
    [g,z] = calcGrad(v,rho,z);
end
ub = rho;
end

function [g,x] = calcGrad(v,rho,warm_start)
x = solve_with_rho(v,rho,warm_start(3));
if (x(2)==0)
    g = x(1);
else
    g = (x(1) + x(2)*log(x(2)/x(3)));
end
end


function x = solve_with_rho(v,rho,w)
x = zeros(3,1);
x(3) = newton_exp_onz(rho,v(2),v(3),w);
x(2) = (1/rho)*(x(3) - v(3))*x(3);
x(1) = v(1) - rho;
end


function z = newton_exp_onz(rho, y_hat, z_hat,w)
t = max(max(w - z_hat, -z_hat),1e-6);
for iter=1:1e2;
    f = (1/rho^2)*t*(t + z_hat) - y_hat/rho + log(t/rho) + 1;
    fp = (1/rho^2)*(2*t + z_hat) + 1/t;
    
    t = t - f/fp;
    if (t <= -z_hat)
        t = -z_hat;
        break;
    elseif (t <= 0)
        t = 0;
        break;
    elseif (abs(f)<1e-6)
        break;
    end
end
z = t + z_hat;
end

function z = symmetrizeSDP(z,K)
l = K.f + K.l;
for i=1:length(K.q)
    l = l + K.q(i);
end
for i=1:length(K.s)
    V = reshape(z(l+1:l+K.s(i)^2),K.s(i),K.s(i));
    V = (V+V')/2;
    z(l+1:l+K.s(i)^2) = V(:);
    l=l+K.s(i)^2;
end
end