n = 200;

sol = randn(n,1);   %One feasible solution
Q = sprandn(n,n,.1); Q = Q'*Q;
A = randn(100,n);
b = A*sol;

for doLinesearch = 0:1
    cvx_begin
        variable x(n)
        cvx_solver('scs_matlab')
        cvx_solver_settings('line_search',doLinesearch,'eps',1e-4,'alpha',1.5)
        minimize(quad_form(x,Q)+norm(x,1))
        subject to
            A * x == b
    cvx_end
    sprintf('Cost: %5.2e\nFeasiblility %5.2e \n', x'*Q*x+norm(x,1), norm(A*x-b))
end
