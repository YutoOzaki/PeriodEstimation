function [P, z, b] = PIE(x, Pmax, lmd)
    %%
    T = numel(x);
    d = 1:Pmax;
    ephi = eulerPhi(d);
    g = sum(ephi);
    
    %% Binary mask
    w = ones(T, 1);
    w(isnan(x)) = 0;
    x(isnan(x)) = 0;
    
    %% Ramanujan periodic dictionary
    A = [];
    N = (0:(T - 1))';
    for i=1:numel(d)
        C = ramsum(N, d(i));
        P = arrayfun(@(j) circshift(C, j), 0:(ephi(i) - 1), 'UniformOutput', false)';
        P = cat(2, P{:});
        A = [A, P];
    end
    
    A = round(real(A));
    
    %% Matrix H
    iH = sparse(g, g);
    p = 0;
    for i=1:numel(ephi)
        for j=0:(ephi(i) - 1)
            p = p + 1;
            iH(p, p) = 1/d(i)^2;
        end
    end
    
    %% Block diagonal adjacency matrix
    G = sparse(g, g);
    p = 1;
    q = 1;
    for i=1:numel(d)
        k = ephi(i);
        G(p:p + k - 1, q:q + k - 1) = 1;
        p = p + k;
        q = q + k;
    end
    
    %% Laplacian matrix of G
    L = sparse(g, g);
    df = sum(G, 1);
    for i=1:g
        L(i, i) = df(i);
    end
    L = L - G;
    
    %% Algorithm
    rho = 1;
    z = x;
    y = zeros(g, 1);
    al = zeros(g, 1);
    be = zeros(g, 1);
    Ahat = A*iH;
    I = sparse(g, g);
    for i=1:g
        I(i, i) = 1;
    end 
    Q = inv(Ahat'*Ahat + 2*lmd(3).*L + rho.*I);
    
    O = ones(T, 1);
    C1 = (2*lmd(1)).*w.*x;
    C2 = O + 2*lmd(1).*w;
    
    eps_b = 1e-4;
    L_z = 100;
    L_z_old = 0;
    eps_z = 1e-4;
    
    while abs(1 - L_z/L_z_old) > eps_z
        L_b = 100;
        L_b_old = 0;
    
        while abs(1 - L_b/L_b_old) > eps_b
            be = Q*(Ahat'*z + rho.*al - y);
            v = be + y./rho;
            al = max(abs(v) - lmd(2)/rho, 0).*sign(v);
            y = y + (be - al);
            
            L_b_old = L_b;
            L_b = 0.5*sum((z - Ahat*be).^2) + lmd(2)*sum(abs(be)) + lmd(3)*be'*L*be;
        end
        
        z = (C1 + Ahat*be)./C2;
        
        L_z_old = L_z;
        L_z = 0.5*sum((z - Ahat*be).^2) + lmd(1)*sum((w.*(x - z)).^2) + ...
            lmd(2)*sum(abs(be)) + lmd(3)*be'*L*be;
        fprintf('%s: L_z = %3.4f\n', datetime, L_z);
    end
    
    b = iH*al;
    
    %% Period strength
    P = zeros(Pmax, 1);
    p = 1;
    for i=1:Pmax
        k = ephi(i);
        P(i) = sum(b(p:p + k - 1).^2);
        p = p + k;
    end
end