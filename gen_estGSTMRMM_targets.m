function [estGSTM] = gen_estGSTMRMM_targets(ett_measurements,ett_ground_truth)
K=150;
nz=2;
nx=4;
F = [eye(2) eye(2);zeros(2) eye(2)];
H =[eye(2) zeros(2)];
Q0 =0.1*[1/3*eye(2) 1/2*eye(2);1/2*eye(2) eye(2)];
R0=2*eye(2);%2
xkk= [100;100;0;0];     % [25;25;12;12]
Pkk=10*diag([10 10 1 1]) ;    %0.5*diag([10 10 1 1])
extentkk = [20 0; 0 20];

IWvkk=60;
IWVkk=extentkk*(IWvkk-2*nz-2);
tao=50;
for t=1:K 
    xk1k=F*xkk;
    Pk1k=F*Pkk*F'+Q0;
    IWvk1k=exp(-1/tao)*IWvkk;
    IWVk1k=(exp(-1/tao)*IWvkk-nz-1)*IWVkk/(IWvk1k-nz-1);
    %避免不正定
    IWVk1k = (IWVk1k + IWVk1k')/2;
    min_eig = min(eig(IWVk1k));
    if min_eig <= 0
        IWVk1k = IWVk1k + (abs(min_eig) + 1e-8) * eye(size(IWVk1k));
    end
    %初始化
    q0=0.5;
    g0=0.5;
    E_extentkk=IWVk1k/(IWvk1k-2*nz-2);
    E_i_extentkk=IWvk1k*inv(IWVk1k);
    E_log_taok=psi(q0)-psi(1);
    E_log_1_taok=psi(q0)-psi(1);
    E_log_paik=psi(g0)-psi(1);
    E_log_1_paik=psi(g0)-psi(1);
    E_thitak=0.9;
    E_log_thitak=0;
    E_pk=1;
    E_log_pk=0;
    E_tk=0.8;
    E_gamak=0.8;
    c=0.25;

    nk=10;
    ZK=ett_measurements{t};
    zk_bar1=mean(ZK');
    zk_bar=zk_bar1';
    Zk_bar=(ZK-repmat(zk_bar,1,nk))*(ZK-repmat(zk_bar,1,nk))';
    N=10;
    w=5;
    v=5;
     for i=1:10
        %估计xkk，Pkk
        xkk_i=xk1k;
        Pk1k_hat=Pk1k/(E_tk+(1-E_tk)*E_thitak);
        R_hat=R0/( E_gamak+(1- E_gamak)*E_pk);
        Sk=H*Pk1k_hat*H'+(c*E_extentkk+R_hat)/nk;
        Kk=Pk1k_hat*H'*inv(Sk);
        xkk=xk1k+Kk*(zk_bar-H*xk1k);
        Pkk=Pk1k_hat-Kk*Sk*Kk';
        Pkk = (Pkk + Pkk')/2;
        min_eig = min(eig(Pkk));
        if min_eig <= 0
            Pkk = Pkk + (abs(min_eig) + 1e-8) * eye(size(Pkk));
        end
        %估计扩展状态
        Nk=(zk_bar-H*xk1k)*(zk_bar-H*xk1k)';
        IWvkk=IWvk1k+nk;
        IWVkk=IWVk1k+inv(Sk)*Nk+Zk_bar;
        E_extentkk=IWVkk/(IWvkk-2*nz-2);
         %估计tk和yk
        Ak=Pkk+(xkk-xk1k)*(xkk-xk1k)';
        Bk=H*Pkk*H'+(zk_bar-H*xkk)*(zk_bar-H*xkk)';
        ptk1=exp( E_log_taok-0.5*trace(Ak*inv(Pk1k)));
        ptk0=exp(E_log_1_taok-0.5*nx* E_log_thitak-0.5*trace(Ak*E_thitak*inv(Pk1k)));
        if  ptk1<=1e-98
            ptk1= ptk1+1e-98;
        end
        if  ptk1>=0.98
            ptk1= 0.98;
        end
        if  ptk0<=1e-98
            ptk0= ptk0+1e-98;
        end
        if  ptk0>=0.98
            ptk0= 0.98;
        end
        E_tk= ptk1/(ptk1+ptk0);
        if E_tk<=1e-98
            E_tk= E_tk+1e-98;
        end
        if  E_tk>=0.98
            E_tk= 0.98;
        end
       
        pgamak1=exp(E_log_paik-0.5*trace(Bk*inv(R0)));
        pgamak0=exp(E_log_1_paik-0.5*nz*E_log_pk-0.5*trace(Bk*E_pk*inv(R0)));
        if  pgamak1<=1e-98
            pgamak1= pgamak1+1e-98;
        end
        if  pgamak1>=0.98
            pgamak1= 0.98;
        end
         if  pgamak0<=1e-98
            pgamak0= pgamak0+1e-98;
        end
        if  pgamak0>=0.98
            pgamak0= 0.98;
        end
        E_gamak=pgamak1/(pgamak1+pgamak0);
        if  E_gamak<=1e-98
            E_gamak= E_gamak+1e-98;
        end
        if  E_gamak>=0.98
           E_gamak= 0.98;
        end
         %估计thitak和pk
         arfa1=0.5*nx*(1-E_tk)+0.5*w;
         arfa2=0.5*(1-E_tk)*trace(Ak*inv(Pk1k))+0.5*w;
         beita1=0.5*nz*(1- E_gamak)+0.5*v;
         beita2=0.5*(1-E_gamak)*trace(Bk*inv(R0))+0.5*v;
         E_thitak=arfa1/arfa2;
         E_pk=beita1/beita2;
         E_log_thitak=psi(arfa1)-log(arfa2);
         E_log_pk=psi(beita1)-log(beita2);
         %估计taok和pai
         gkk=q0+E_tk;
         hkk=2-q0-E_tk;
         ekk=g0+E_gamak;
         fkk=2-g0-E_gamak;
         E_log_taok=psi(gkk)-psi(2);
         E_log_1_taok=psi(hkk)-psi(2);
         E_log_paik=psi(ekk)-psi(2);
         E_log_1_paik=psi(fkk)-psi(2);
     end
     estextents(:, :, t)=E_extentkk;
     eststates(:, t)= xkk;
     estGSTM.states = eststates;
     estGSTM.extents = estextents;
end
%  hold on;
%     grid on;
%     box on;
%     axis equal;
%     xlabel('X位置');
%     ylabel('Y位置');
%     title('所有时刻的估计和真实扩展');
%     
%     % 创建颜色映射，不同时间步用不同颜色
%     colors = lines(K);  % 使用distinct colors
%     
%     % 循环绘制每个时间步的数据
%     for k = 1:K
%         % 提取当前时间步数据
%         Y_k = ett_measurements{k};
%         center = ett_ground_truth.states(1:2, k);
%         extent = ett_ground_truth.extents(:,:,k);
%         estcenter = eststates(1:2, k);
%         estextent = estextents(:,:,k);
%         
%         % 绘制测量点（带透明度）
%         scatter(Y_k(1, :), Y_k(2, :), 40, colors(1, :), 'filled', ...
%                'MarkerFaceAlpha', 0.6, ...
%                'MarkerEdgeColor', 'none');
%         
%         % 绘制扩展椭圆
%         plot_ellipse(extent, center, 1, ...
%             'Color', colors(2, :), ...
%             'LineWidth', 1.5);
%         plot_ellipse(estextent, estcenter, 1, ...
%             'Color', colors(4, :), ...
%             'LineWidth', 1.5);
%         
%         % 绘制目标中心位置
%         plot(estcenter(1), center(2), 'o', ...
%             'MarkerFaceColor', colors(2, :), ...
%             'MarkerEdgeColor', 'k', ...
%             'MarkerSize', 6);
%         plot(estcenter(1), estcenter(2), 'o', ...
%             'MarkerFaceColor', colors(4, :), ...
%             'MarkerEdgeColor', 'k', ...
%             'MarkerSize', 6);
%     end
%     
% %     % 添加颜色条显示时间步
% %     colormap(colors);
% %     caxis([1, K]);
% %     cbar = colorbar;
% %     cbar.Label.String = '时间步';
% %     cbar.Ticks = unique(round(linspace(1, K, min(10, K)))); % 最多10个刻度
%     
%     % 添加图例
%     h = zeros(2, 1);
%     h(1) = plot(NaN, NaN, 'ko', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
%     h(2) = plot(NaN, NaN, 'k-', 'LineWidth', 1.5);
%     legend(h, {'目标中心位置', '真实扩展'}, 'Location', 'best');
%     
%     % 增强视觉效果
%     set(gca, 'FontSize', 12);
%     set(gcf, 'Color', 'w');
%     
%     % 添加数据统计信息
%     num_measurements = sum(cellfun(@(y) size(y, 2), ett_measurements));
%     text(0.02, 0.98, sprintf('总测量点数: %d', num_measurements), ...
%          'Units', 'normalized', 'FontSize', 12, ...
%          'BackgroundColor', [1, 1, 1, 0.7]); 
%      function plot_ellipse(M, center, scale, varargin)
%     % PLOT_ELLIPSE 绘制扩展椭圆
%     % 输入:
%     %   M: 2x2协方差或形状矩阵
%     %   center: [x; y] 椭圆中心
%     %   scale: 缩放因子
%     %   varargin: 绘图属性
%     
%     % 确保形状矩阵正定
%     M = make_positive_definite(M);
%     
%     % 特征分解
%     [V, D] = eig(M);
%     eigen_values = diag(D);
%     
%     % 参数化角度
%     theta = linspace(0, 2*pi, 200);
%     
%     % 计算椭圆坐标
%     ellipse_points = V * sqrt(scale * D) * [cos(theta); sin(theta)] + center;
%     
%     % 绘制椭圆
%     plot(ellipse_points(1, :), ellipse_points(2, :), varargin{:});
%     
%     % 添加主轴方向（50%长度）
%     max_eig = max(eigen_values);
%     for i = 1:2
%         if eigen_values(i) > 1e-6 % 忽略极小特征值
%             axis_direction = V(:, i) * sqrt(eigen_values(i)) * sqrt(scale) * 0.5;
%             start_point = center - axis_direction;
%             end_point = center + axis_direction;
%             plot([start_point(1), end_point(1)], [start_point(2), end_point(2)], ...
%                  'k--', 'LineWidth', 1);
%         end
%     end
% end
% 
% function M = make_positive_definite(M)
%     % MAKE_POSITIVE_DEFINITE 确保矩阵是正定的
%     min_eig = min(eig(M));
%     if min_eig < 1e-6
%         % 如果最小特征值太小，添加小量使其正定
%         M = M + (abs(min_eig) + 1e-6) * eye(size(M));
%     end
% end
% function plot_handle = drawEllipse(M, center, n)
% % drawEllipse:
% %   Draws an ellipse defined by a positive-definite matrix M, centered at 'center',
% %   scaled by factor n.
% %
% % Inputs:
% %   M: 2x2 covariance or shape matrix defining the ellipse
% %   center: [x; y] center of the ellipse
% %   n: scaling factor (e.g., 1 for 1-sigma ellipse)
% %
% % Output:
% %   plot_handle: handle to the generated plot object
% 
%     [V,D] = eig(M); % Eigen decomposition to get principal axes
%     t = linspace(0, 2*pi, 100); % Parametric angle
%     ellipse = V * sqrt(D) * [n*cos(t); n*sin(t)] + center; 
%     plot_handle = plot(ellipse(1,:), ellipse(2,:));
% end


end
