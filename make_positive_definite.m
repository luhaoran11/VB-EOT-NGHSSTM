function M = make_positive_definite(M)
    % MAKE_POSITIVE_DEFINITE 确保矩阵是正定的
    min_eig = min(eig(M));
    if min_eig < 1e-6
        % 如果最小特征值太小，添加小量使其正定
        M = M + (abs(min_eig) + 1e-6) * eye(size(M));
    end
end