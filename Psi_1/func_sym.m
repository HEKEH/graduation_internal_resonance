syms Ups_ck1 xc omega_m Ups_ck2 Ups_ck3 Ups_ck4 omega_n Ups_ck5 xb Ups_bk1 Ups_bk2 Ups_bk3 Ups_bk4 alpha_bk beta_bk k1 k2 k3 k4 Gamma_k1 Gamma_k2 theta
psi_ck_known_sym = Ups_ck1 * xc * cos(omega_m * xc) + Ups_ck2 * xc * sin(omega_m * xc)...
    + Ups_ck3 * cos(omega_n * xc) + Ups_ck4 * sin(omega_n * xc) + Ups_ck5;
psi_ck_known_sym_p = diff(psi_ck_known_sym, xc);
psi_ck_known_sym_p2 = diff(psi_ck_known_sym, xc, 2);
psi_ck_known_sym_int = int(psi_ck_known_sym, xc, 0, 1);

psi_bk_known_sym = Ups_bk1 * xb .* cos(alpha_bk * xb) + Ups_bk2 * xb .* sin(alpha_bk * xb)...
    + Ups_bk3 * xb .* cosh(beta_bk * xb) + Ups_bk4 * xb .* sinh(beta_bk * xb);
psi_bk_known_sym_p = diff(psi_bk_known_sym, xb);
psi_bk_known_sym_p2 = diff(psi_bk_known_sym, xb, 2);
psi_bk_known_sym_p3 = diff(psi_bk_known_sym, xb, 3);

psi_ck_uncertain_sym = [cos(omega_m * xc) + Gamma_k1, sin(omega_m * xc) + Gamma_k2, 0, 0, 0, 0];
psi_ck_uncertain_sym_p = diff(psi_ck_uncertain_sym, xc);
psi_ck_uncertain_sym_p2 = diff(psi_ck_uncertain_sym, xc, 2);
psi_ck_uncertain_sym_int = int(psi_ck_uncertain_sym, xc, 0, 1);
psi_bk_uncertain_sym = [0, 0, cos(alpha_bk * xb), sin(alpha_bk * xb), cosh(beta_bk * xb), sinh(beta_bk * xb)];
psi_bk_uncertain_sym_p = diff(psi_bk_uncertain_sym, xb);
psi_bk_uncertain_sym_p2 = diff(psi_bk_uncertain_sym, xb, 2);
psi_bk_uncertain_sym_p3 = diff(psi_bk_uncertain_sym, xb, 3);

matr_A = sym(zeros(6, 6));
vec_b = sym(zeros(6, 1));

matr_A(1, :) = subs(psi_ck_uncertain_sym, xc, 0);
vec_b(1) = -subs(psi_ck_known_sym, xc, 0);

matr_A(2, :) = subs(psi_bk_uncertain_sym, xb, 0);
vec_b(2) = -subs(psi_bk_known_sym, xb, 0);

matr_A(3, :) = subs(psi_bk_uncertain_sym_p, xb, 0);
vec_b(3) = -subs(psi_bk_known_sym_p, xb, 0);

matr_A(4, :) = subs(psi_bk_uncertain_sym_p2, xb, 1); 
vec_b(4) = -subs(psi_bk_known_sym_p2, xb, 1);

matr_A(5, :) = subs(psi_ck_uncertain_sym, xc, 1) - subs(psi_bk_uncertain_sym, xb, 1) * cos(theta)^2;
vec_b(5) = -subs(psi_ck_known_sym, xc, 1) + subs(psi_bk_known_sym, xb, 1) * cos(theta)^2;

matr_A(6, :) = -subs(psi_bk_uncertain_sym_p3, xb, 1) + k1 * subs(psi_ck_uncertain_sym_p, xc, 1) + k2 * psi_ck_uncertain_sym_int + k3 * subs(psi_bk_uncertain_sym, xb, 1) + k4 * subs(psi_bk_uncertain_sym_p, xb, 1);
vec_b(6) = subs(psi_bk_known_sym_p3, xb, 1) - k1 * subs(psi_ck_known_sym_p, xc, 1) - k2 * subs(psi_ck_known_sym_int, xc, 1) - k3 * subs(psi_bk_known_sym, xb, 1) - k4 * subs(psi_bk_known_sym_p, xb, 1);
C_k_sym = matr_A \ vec_b;