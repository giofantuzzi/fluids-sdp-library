%% couette_bounds.m
%
% Compute bounds on energy dissipation for 3D plane Couette flow using the
% background method. For a detailed description of the problem, see
%
% http://quinopt.readthedocs.io/04_examples/planeCouetteBF.html
%
% G. Fantuzzi, 3 June 2023

%% Initialization
%
% First, we remove any existing variables to avoid unexpected dependencies, and
% clear YALMIP's and QUINOPT's internals
clear;
yalmip clear;
quinopt clear;

% Then,we set some problem parameters: the Reynolds number, the period Lambda_y
% in the y direction, the degree of the linear term PHI in the auxiliary
% functional, and the maximum horizontal wavenumber to test:
ReVals   = 500;%ceil(logspace(2,3,100));
Lambda_y = 4*pi;
PHI_DEG  = 50;
k_max    = 10;
savedata = false;

for ii = 1:length(ReVals)

    % Reset
    Re = ReVals(ii);
    fprintf('Re = %6.4e\n',Re)
    yalmip clear;
    quinopt clear;

    % Finally, we define the integration variables, the flow variables, and the
    % boundary conditions
    z = indvar(0,1);
    [u,w] = depvar(z);
    BC = [u(0); u(1); w(0); w(1); w(0,1); w(1,1)];

    %% Setting up the optimization variables
    % The optimization variables are:
    % 1) a, a scalar known in the literature as "balance parameter"
    % 2) U, the upper bound on the time-average bulk dissipation
    % 3) PHI, a polynomial of degree PHI_DEG satisfying the conditions PHI(0)=0 and
    %    PHI(1)=0. We set up PHI as a polynomial in the Legendre basis, which is how
    %    QUINOPT represents variable internally.
    parameters a U
    C_DPHI = sdpvar(PHI_DEG-1,1);
    D1PHI = legpoly(z,PHI_DEG-1,[0; C_DPHI]);
    PHI = legpolyint(D1PHI,z);

    % Finally, we compute the second derivative of PHI
    D2PHI = jacobian(D1PHI,z);

    %% Setting up the inequality constraints
    % To set up the integral inequality constraints, we construct a vector EXPR
    % containing the integrand of each inequality.
    n = 0;
    k = 0;
    EXPR(1) = (a-1)*u(z,1)^2 + D2PHI*u(z) + U-1;
    while k<k_max
        n = n+1;
        k = 2*pi*n/Lambda_y;
        EXPR(end+1) = (a-1)*( u(z,1)^2 + k^2*u(z)^2 ) ...
            +(a-1)*( w(z,2)^2/k^2 + 2*w(z,1)^2 + k^2*w(z)^2 ) ...
            + Re*( a+D1PHI )*u(z)*w(z);
    end

    %% Solve the problem
    % The aim of the optimization is to minimize the upper bound U on the energy
    % dissipation. To specify the extra constraints on the optimization variables,
    % we use the command quinopt() with five inputs: EXPR and BC to specify the
    % integral inequalities, U as the objective, an options structure,  and
    % the additional constraints PHI_BC. We use the "outer" option to obtain an
    % easier problem without estimates for the discretization error.
    OPTS.method = 'outer';
    OPTS.N = 50;
    OPTS.YALMIP = sdpsettings('solver','mosek');
    OPTS.YALMIP.savesolveroutput = true;
    OPTS.solve = false;
    [SOL,CNSTR,DATA] = quinopt(EXPR,BC,U,OPTS);
    if savedata
        SDP = export(CNSTR,U,sdpsettings('solver','sedumi'));
        SDP.At = SDP.A;
        SDP = rmfield(SDP,'A');
        FNAME = sprintf('couette-Re%i.mat',Re);
        save(FNAME,'Re','SDP')
    end
end

%% END CODE