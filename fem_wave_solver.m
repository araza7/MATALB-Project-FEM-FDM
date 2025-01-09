function P = fem_wave_solver(ne, L1, L2, c1, c2, rho, f)

    % Derived parametersa
    w = 2 * pi * f; % Angular frequency

    % Total length and discretization
    L = L1 + L2;
    dx = L / ne;
    x = 0:dx:L;
    n_l = 1:ne;
    n_r = 2:ne + 1;

    % Initialize system matrices
    K = sparse(ne + 1, ne + 1);
    M = sparse(ne + 1, ne + 1);
    C = sparse(ne + 1, ne + 1);

    % Assemble stiffness and mass matrices
    for ii = 1:ne
        Le = dx;
        ndof = [n_l(ii) n_r(ii)]; % Degrees of freedom for the element

        if x(ii) <= L1
            % Fluid M1
            [ke1, me1] = keme(Le, rho, c1);
            K(ndof, ndof) = K(ndof, ndof) + ke1 ;
            M(ndof, ndof) = M(ndof, ndof) + me1;
        elseif x(ii)==L1
            M(n_l(ii),n_l(ii-1))=Le/(6*c1^2*rho);
            M(n_l(ii),n_l(ii))=Le/(3*rho*(1/c1^2+1/c2^2));
            M(n_l(ii),n_l(ii+1))=Le/(6*c2^2*rho);
            K(n_l(ii),n_l(ii-1))=-1/(rho*Le);
            K(n_l(ii),n_l(ii))=2/(rho*Le);
            K(n_l(ii),n_l(ii+1))=-1/(rho*Le);
        elseif x(ii)> L1
            % Fluid M2
            [ke1, me2] = keme(Le, rho, c2);
            K(ndof, ndof) = K(ndof, ndof) + ke1 ;
            M(ndof, ndof) = M(ndof, ndof) + me2;
        end

    end

    % Impedance boundary condition at the right end (damping term)
    Z = rho * c2; % Impedance using c2
    C(end, end) = 1 / Z;

    % System matrix assembly with frequency-based scaling
    KK = K + 1i * w * C - w^2 * M;

    % Right-hand side
    %Drichlet conditions on left side
    F = zeros(ne+1, 1);
    KK(1,:) = 0;
    KK(1,1) = 1;
    F(1) = 1;

    % Solve the system for pressure
    P = (KK \ F); % Pressure solution
    
end
