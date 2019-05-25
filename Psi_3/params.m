func_sym;
run('../Configs.m');
run('../Phi/phi_params.m');
C_idx = m_idx;

Gamma1_m = Gamma1(m_idx);
Gamma2_m = Gamma2(m_idx);
Gamma1_n = Gamma1(n_idx);
Gamma2_n = Gamma2(n_idx);
C_m = C(:, m_idx);
C_n = C(:, n_idx);
omega_m = omegas(m_idx);
omega_n = omegas(n_idx);

alpha_bk = sqrt(sqrt(P^2 + 16 * beta1 * omega_m^2) / 2 / beta1 + P / 2 / beta1);
beta_bk = sqrt(sqrt(P^2 + 16 * beta1 * omega_m^2) / 2 / beta1 - P / 2 / beta1);

PHI1 = @(xc)H*lambda^2*phi_c_p2(xc, m_idx)*integral(@(xc)phi_c(xc, m_idx),0, 1);
PHI2 = @(xc)H*lambda^2*phi_c_p2(xc, n_idx)*integral(@(xc)phi_c(xc, n_idx),0, 1);
PHI3 = @(xc)H*lambda^2*(phi_c_p2(xc, m_idx)*integral(@(xc)phi_c(xc, n_idx),0, 1)+phi_c_p2(xc, n_idx)*integral(@(xc)phi_c(xc, m_idx),0, 1));
PHI4 = @(xc)H*lambda^2*(1/2-H*tan(theta))*phi_c_p2(xc, m_idx)*phi_c(1, m_idx);
PHI5 = @(xc)H*lambda^2*(1/2-H*tan(theta))*phi_c_p2(xc, n_idx)*phi_c(1, n_idx);
PHI6 = @(xc)H*lambda^2*(1/2-H*tan(theta))*(phi_c_p2(xc, m_idx)*phi_c(1, n_idx)+phi_c_p2(xc, n_idx)*phi_c(1, m_idx));
I1 = H*lambda^2*integral(@(xc)phi_c_p(xc, m_idx).^2,0, 1);
I2 = H*lambda^2*integral(@(xc)phi_c_p(xc, n_idx).^2,0, 1);
I3 = 2*H*lambda^2*integral(@(xc)phi_c_p(xc, m_idx).*phi_c_p(xc, n_idx),0, 1);

J_m = integral(@(xc)phi_c(xc, m_idx) .* (PHI2(xc) - PHI5(xc) - I2), 0, 1);
J_n = integral(@(xc)phi_c(xc, n_idx) .* (PHI3(xc) - PHI6(xc) - I3), 0, 1);

Ups_ck1 = (H*lambda^2*C_m(1)*(C_m(1)*(2*sin(omega_m)+omega_m*(cos(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma1_m))+C_m(2)*(2-2*cos(omega_m)+omega_m*(sin(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma2_m))))/(6*omega_m);

Ups_ck2 = (H*lambda^2*C_m(2)*(C_m(1)*(2*sin(omega_m)+omega_m*(cos(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma1_m))+C_m(2)*(2-2*cos(omega_m)+omega_m*(sin(omega_m)*(-1+2*H*tan(theta))+(1+2*H*tan(theta))*Gamma2_m))))/(6*omega_m);

temp1 = (lambda^2*(2*sin(omega_m)+cos(omega_m)*omega_m*(-1+2*H*tan(theta))))/(8*omega_m^3-lambda^2*omega_m*(1+2*H*tan(theta)));

temp2 = (lambda^2*(2-2*cos(omega_m)+sin(omega_m)*omega_m*(-1+2*H*tan(theta))))/(8*omega_m^3-lambda^2*omega_m*(1+2*H*tan(theta)));

temp3 = (H*lambda^2*omega_m*((sin(2*omega_m)-2*omega_m)*C_m(1)^2+4*sin(omega_m)^2*C_m(1)*C_m(2)-(sin(2*omega_m)+2*omega_m)*C_m(2)^2))/(2*(-8*omega_m^2+lambda^2*(1+2*H*tan(theta))));

Ups_ck3 = temp1 * Ups_ck1 + temp2 * Ups_ck2 + temp3;

Gamma_k1 = (2*lambda^2*cos(omega_m)*sin(omega_m)-lambda^2*cos(2*omega_m)*omega_m+2*H*lambda^2*cos(2*omega_m)*omega_m*tan(theta))/(omega_m*(-lambda^2+8*omega_m^2-2*H*lambda^2*tan(theta)));

Gamma_k2 = (2*lambda^2*sin(omega_m)^2-lambda^2*sin(2*omega_m)*omega_m+2*H*lambda^2*sin(2*omega_m)*omega_m*tan(theta))/(omega_m*(-lambda^2+8*omega_m^2-2*H*lambda^2*tan(theta)));