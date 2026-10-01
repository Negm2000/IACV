function K = compute_calibration_matrix(v_vert_n, v_axis_n, v_trans_n, img)
% COMPUTE_CALIBRATION_MATRIX - Compute camera calibration matrix K
%
% Theory:
%   - With 3 orthogonal vanishing points: 3 constraints (vi^T * ω * vj = 0)
%   - Assumes zero skew (s=0) and square pixels (fx = fy)
%   - This gives 4 constraints for 4 unknowns in IAC

[rows, cols, ~] = size(img);
cx = cols / 2;
cy = rows / 2;
scale = max(cx, cy);
T_n2p = [scale 0 cx; 0 scale cy; 0 0 1];

v1 = v_vert_n(:);
v2 = v_axis_n(:);
v3 = v_trans_n(:);

% IAC parametrization (zero skew: ω12 = 0)
%   ω = [w1,  0, w4]
%       [ 0, w3, w5]
%       [w4, w5, w6]
%
% Orthogonality: vi^T * ω * vj = 0
% => w1*vi1*vj1 + w3*vi2*vj2 + w4*(vi1*vj3+vi3*vj1) + w5*(vi2*vj3+vi3*vj2) + w6*vi3*vj3 = 0

build_row = @(vi, vj) [vi(1)*vj(1), vi(2)*vj(2), vi(1)*vj(3)+vi(3)*vj(1), vi(2)*vj(3)+vi(3)*vj(2)];
build_rhs = @(vi, vj) -vi(3)*vj(3);

A = zeros(4, 4);
b = zeros(4, 1);

% Three orthogonality constraints
A(1,:) = build_row(v1, v2);  b(1) = build_rhs(v1, v2);
A(2,:) = build_row(v1, v3);  b(2) = build_rhs(v1, v3);
A(3,:) = build_row(v2, v3);  b(3) = build_rhs(v2, v3);

% Square pixels assumption: fx = fy => w1 = w3
A(4,:) = [1, -1, 0, 0];  b(4) = 0;

% Solve
x = A \ b;
w1 = x(1); w3 = x(2); w4 = x(3); w5 = x(4); w6 = 1;

omega_n = [w1, 0, w4; 0, w3, w5; w4, w5, w6];

% Extract K via Cholesky: ω = K^(-T) * K^(-1)
try
    L = chol(omega_n, 'lower');
    K_norm = inv(L');
    K_norm = K_norm / K_norm(3,3);
catch
    % Handle non-positive-definite case
    [U, D, ~] = svd(omega_n);
    D_pos = abs(D);
    omega_fixed = U * D_pos * U';
    L = chol(omega_fixed, 'lower');
    K_norm = inv(L');
    K_norm = K_norm / K_norm(3,3);
end

% Denormalize to pixel coordinates
K = T_n2p * K_norm;
K = K / K(3,3);

% Ensure positive focal lengths
if K(1,1) < 0
    K = -K;
    K = K / K(3,3);
end
end
