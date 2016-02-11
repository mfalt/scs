h = 0.1;
A0 = [0.39 0.0709 0.61 0.0223; -8.01 0.293 8.01 0.5693; 0.672 0.025 0.328 0.0667; 8.68 0.638 -8.68 0.206];
A1 = [0.979 0.0928 0.0207 0.000669; -0.402 0.853 0.402 0.0195; 0.0229 0.00075 0.977 0.0907; 0.440 0.0218 -0.44 0.812];
B0 = [0.00182; 0.031; 0.000296; 0.0109];
B1 = [0.00208; 0.0405; 0.00000833; 0.000327];
Bb = [0; 1; 0; 0];
Qx = h*diag([100 0 100 0]);
QNx= eye(4);
Qu = h*diag([1e-3 0 0]);
QNu= Qu;
Qw = 0;
QNw= 0;

T = 20;
x0 = [1 0 1 0];
M = 1000;
m = -1000;

%12 elements [x1 x2 x3 x4 zx1 zx2 zx3 zx4 zu u1 u2 u3]
% Au2 x + Bu2 u1 + Bb u3 = A0 x + B0 u1 + (A1-A0)zx + (B1-B0) zu + Bb u3
EqBlock = [-eye(4) zeros(4,8) A0 (A1-A0) (B1-B0) B0 zeros(4,1) Bb];
AEq = zeros(4*(T+1),12*T);
AEq(1:4,1:4) = eye(4);
for i = 1:T
    AEq((4*i+1):4*(i+1),(12*(i-1)+1):(12*i+12)) = EqBlock;
end
bEq = zeros(4*(T+1),1);
bEq(1:4,1) = x0;

Aineq = zeros(0,12);
bineq = zeros(0,1);
Aineq = [Aineq; eye(4) zeros(4,8); -eye(4) zeros(4,8)];
bineq = [bineq; 10; 100; 10; 100; 10; 100; 10; 100];
Aineq = [Aineq; zeros(1,9) 1 0 0; zeros(1,9) -1 0 0];
bineq = [bineq; 1; 1];

Aineq = [Aineq; zeros(1,11) 1; zeros(1,11) -1];
bineq = [bineq; 1; 0];

%Fix zx
Aineq = [Aineq; zeros(4,4) eye(4) zeros(4,2) -M*ones(4,1) zeros(4,1)];
bineq = [bineq; zeros(4,1)];
Aineq = [Aineq; zeros(4,4) -eye(4) zeros(4,2) m*ones(4,1) zeros(4,1)];
bineq = [bineq; zeros(4,1)];
Aineq = [Aineq; eye(4) -eye(4) zeros(4,2) M*ones(4,1) zeros(4,1)];
bineq = [bineq; M*ones(4,1)];
Aineq = [Aineq; -eye(4) eye(4) zeros(4,2) -m*ones(4,1) zeros(4,1)];
bineq = [bineq; -m*ones(4,1)];

%TODO FIx this
%Fix zu
Aineq = [Aineq; zeros(4,4) eye(4) zeros(4,2) -M*ones(4,1) zeros(4,1)];
bineq = [bineq; zeros(4,1)];
Aineq = [Aineq; zeros(4,4) -eye(4) zeros(4,2) m*ones(4,1) zeros(4,1)];
bineq = [bineq; zeros(4,1)];
Aineq = [Aineq; eye(4) -eye(4) zeros(4,2) M*ones(4,1) zeros(4,1)];
bineq = [bineq; M*ones(4,1)];
Aineq = [Aineq; -eye(4) eye(4) zeros(4,2) -m*ones(4,1) zeros(4,1)];
bineq = [bineq; -m*ones(4,1)];

ineqs = size(Aineq,1);

CA = zeros(ineqs*(T+1),12*(T+1));
for i = 1:T+1
    CA((ineqs*(i-1)+1):ineqs*i,(12*(i-1)+1):12*i) = Aineq;
end
CB = repmat(bineq,T+1,1);

Qblock = blkdiag(Qx, zeros(5), Qu);
QNblock = blkdiag(QNx, zeros(5), QNu);
Q = zeros(12*(T+1));
for i = 1:T
    Q(12*(i-1)+1:12*i,12*(i-1)+1:12*i) = Qblock;
end
s = 12*(T+1);

fixAblock = [zeros(1,11) 1];
fixBblock = [0];
fixA = zeros(T+1,12);
fixB = zeros(T+1,1);
for i = 1:T+1
    fixA(i,12*(i-1)+1:12*i) = fixAblock;
    fixB(i,1) = fixBblock;
end

fixA2block = [zeros(1,10) 1 0];
fixB2block = [0];
fixA2 = zeros(T+1,12);
fixB2 = zeros(T+1,1);
for i = 1:T+1
    fixA2(i,12*(i-1)+1:12*i) = fixA2block;
    fixB2(i,1) = fixB2block;
end

cvx_begin
    variable x(s)
    minimize(quad_form(x,Q))
    subject to
        AEq * x == bEq
        fixA * x == fixB
        fixA2 * x == fixB2
        CA*x <= CB
cvx_end
