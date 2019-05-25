get_am_bar = @(a_m, a_n, gamma1) sqrt(a_m.^2+(cos(gamma1).*a_m.*a_n.^2.*((cos(theta).^3.*integral(@(xb)phi_b(xb,m_idx).*psi_bks{1}(xb), 0, 1)+m.*integral(@(xc)phi_c(xc,m_idx).*psi_cks{1}(xc), 0, 1))))+(1./4.*a_n.^4.*((cos(theta).^3.*integral(@(xb)psi_bks{1}(xb).^2, 0, 1)+m.*integral(@(xc)psi_cks{1}(xc).^2, 0, 1)))));
get_an_bar = @(a_m, a_n, gamma1) sqrt(a_n.^2+(cos(gamma1).*a_m.*a_n.^2.*((cos(theta).^3.*integral(@(xb)phi_b(xb,n_idx).*psi_bks{2}(xb), 0, 1)+m.*integral(@(xc)phi_c(xc,n_idx).*psi_cks{2}(xc), 0, 1))))+(1./4.*a_m.^2.*a_n.^2.*((cos(theta).^3.*integral(@(xb)psi_bks{2}(xb).^2, 0, 1)+m.*integral(@(xc)psi_cks{2}(xc).^2, 0, 1)))));

a_bar_cell = cell(2, size(a_cell, 2));
for idx = 1: size(a_cell, 2)
    a_bar_cell{1, idx} = get_am_bar(a_cell{1, idx}, a_cell{2, idx}, gamma1_angle{idx});
    a_bar_cell{2, idx} = get_an_bar(a_cell{1, idx}, a_cell{2, idx}, gamma1_angle{idx});
end