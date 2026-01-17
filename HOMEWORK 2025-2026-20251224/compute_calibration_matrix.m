function K = compute_calibration_matrix(v_vert_n, v_axis_n, v_trans_n, img)
% COMPUTE_CALIBRATION_MATRIX - Compute camera calibration matrix K
%   Uses orthogonality constraints on vanishing points.

[rows, cols, ~] = size(img);
cx = cols / 2;
cy = rows / 2;
scale = max(cx, cy);
T_n2p = [scale 0 cx; 0 scale cy; 0 0 1];

v1 = v_vert_n(:);
v2 = v_axis_n(:);
v3 = v_trans_n(:);

% Build constraint matrix for omega (IAC)
build_constraint = @(vi, vj) [vi(1)*vj(1), vi(2)*vj(2), ...
    vi(1)*vj(3)+vi(3)*vj(1), vi(2)*vj(3)+vi(3)*vj(2)];

A = zeros(4, 4);
b = zeros(4, 1);

A(1,:) = build_constraint(v1, v2); b(1) = -v1(3)*v2(3);
A(2,:) = build_constraint(v1, v3); b(2) = -v1(3)*v3(3);
A(3,:) = build_constraint(v2, v3); b(3) = -v2(3)*v3(3);
A(4,:) = [1, -1, 0, 0]; b(4) = 0; % Square pixel assumption

x = A \ b;
w1 = x(1); w3 = x(2); w4 = x(3); w5 = x(4); w6 = 1;

omega_n = [w1, 0, w4; 0, w3, w5; w4, w5, w6];

% Verify positive definiteness
eig_omega = eig(omega_n);
if any(eig_omega <= 0)
    A_alt = A(1:3, :);
    b_alt = b(1:3);
    x_alt = pinv(A_alt) * b_alt;
    w1 = x_alt(1); w3 = x_alt(2); w4 = x_alt(3); w5 = x_alt(4);
    omega_n = [w1, 0, w4; 0, w3, w5; w4, w5, w6];
end

% Extract K from omega via Cholesky
try
    R = chol(omega_n);
    K_norm = inv(R);
    K_norm = K_norm / K_norm(3,3);
catch
    [U, D, ~] = svd(omega_n);
    D(D < 0) = abs(D(D < 0));
    R = U * sqrt(D);
    K_norm = inv(R);
    K_norm = K_norm / K_norm(3,3);
end

% Denormalize
K = T_n2p * K_norm;
K = K / K(3,3);
end
