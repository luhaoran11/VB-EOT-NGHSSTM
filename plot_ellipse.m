  function plot_ellipse(M, center, scale, varargin)
    % PLOT_ELLIPSE 绘制扩展椭圆
    % 输入:
    %   M: 2x2协方差或形状矩阵
    %   center: [x; y] 椭圆中心
    %   scale: 缩放因子
    %   varargin: 绘图属性
    
    % 确保形状矩阵正定
    M = make_positive_definite(M);
    
    % 特征分解
    [V, D] = eig(M);
    eigen_values = diag(D);
    
    % 参数化角度
    theta = linspace(0, 2*pi, 200);
    
    % 计算椭圆坐标
    ellipse_points = V * sqrt(scale * D) * [cos(theta); sin(theta)] + center;
    
    % 绘制椭圆
    plot(ellipse_points(1, :), ellipse_points(2, :), varargin{:});
    
    % 添加主轴方向（50%长度）
    max_eig = max(eigen_values);
    for i = 1:2
        if eigen_values(i) > 1e-6 % 忽略极小特征值
            axis_direction = V(:, i) * sqrt(eigen_values(i)) * sqrt(scale) * 0.5;
            start_point = center - axis_direction;
            end_point = center + axis_direction;
            plot([start_point(1), end_point(1)], [start_point(2), end_point(2)], ...
                 'k--', 'LineWidth', 1);
        end
    end
end