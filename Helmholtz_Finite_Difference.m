function P_numeric = Helmholtz_Finite_Difference(npts, L1, L2, c1, c2, f)
    % Initialize parameters
    w = 2 * pi * f;       % Angular frequency
    dx = (L1 + L2) / (npts - 1);  % Grid spacing
    x = linspace(0, L1 + L2, npts); % x-coordinates
    rho = 1000.0;
    L = L1+L2;
    % Initialize the system matrix (A) and right-hand side (B)
    A = sparse(npts, npts);
    B = sparse(npts, 1);
    
    % Apply Dirichlet boundary condition at the left boundary (x = 0)
    A(1, 1) = 1;  % Dirichlet condition (p = 1 at x = 0)
    B(1) = 1;

    for ii = 2:npts-1
    
        if x(ii) < L1 % In the first medium
            c = c1;
            K_1 = w / c1; % Wave number for the current medium
            K = K_1;
            k1 = 1/dx^2;
            k2 = -2/dx^2 + K^2;
            k3 = 1/dx^2;
            A(ii, ii-1:ii+1) = [k1, k2, k3];
        elseif x(ii) == L1 % at common node
    
            K_m = 1;
            k1 = 1/dx;
            k2 = -2/dx;
            k3 = 1/dx;
            A(ii, ii-1:ii+1) = K_m * [k1, k2, k3];
        elseif x(ii) > L1 && x(ii)<=L  % In the second medium
    
            K_2 = w / c2; % Wave number for the current medium
            K= K_2;
            k1 = 1/dx^2;
            k2 = -2/dx^2 + K^2;
            k3 = 1/dx^2;
            A(ii, ii-1:ii+1) = [k1, k2, k3];
        else
            A(ii, ii-1:ii+1) = NaN;
        end
    
    
    
    end
    
    % z = rho * c2;  % Acoustic impedance for the second medium (Robin condition)
    % impedance_coeff = -1 - 2 * 1i * w * rho / z;
    % A(end, end-1:end) = impedance_coeff * [-1/dx, 1/dx];
    % B(end) = 0;
    Z = rho * c2; % Impedance condition for fluid M2
    A(end, end-1) = -Z/(1i*rho*w*dx);
    A(end, end) = (1 + Z/(1i*rho*w*dx));

    % Solve the system of equations
    P_numeric = A \ B;  % Pressure solution at each grid point

end