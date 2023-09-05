function C = ramsum(g, d)
    C = 0;
    
    for k=1:d
        if h_gcd(d, k) == 1
            C = C + exp((1j*2*pi*k/d).*g);
        end
    end
end

function gcd = h_gcd(a, b)
    r = 1;
    
    while r ~= 0
        r = mod(a, b);
        a = b;
        b = r;
    end
    
    gcd = a;
end

%{
Tenneti, S. V. & Vaidyanathan, P. P. (2015). Nested Periodic Matrices and
Dictionaries: New Signal Representations for Period Estimation. IEEE
TRANSACTIONS ON SIGNAL PROCESSING, 63(14), 3736-3750.

% equations (13) & (14)
for k=1:5
    C = ramsum(0:(k - 1), k);
    disp(real(C));
end

% equations (15) & (16)
P = 4;
q = 4;
phi_q = eulerPhi(q);
c_q = ramsum((0:(P - 1))', q);
C_q = [c_q, circshift(c_q, 1)];
%}