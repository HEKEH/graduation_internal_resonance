K = 100;
m_idx = 6;
n_idx = 4;

F_mati = 1e-5;
Fc = @(xc) F_mati * sin(2 * pi * xc);
Fb = @(xb) F_mati * m * sin(2 * pi * xb);