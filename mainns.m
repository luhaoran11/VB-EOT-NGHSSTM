profile on;
clear all;
close all;
clc; 
nxp=1
K=150
for expt = 1:nxp
 v_val=5;
 t_val=0.8;
 iter_param=10;
 [est_measurements,est_ground_truth] = gen_data_targets()   ;
 [estVB_EOT_NGHSSTM] = gen_estVB_EOT_NGHSSTM_targets(est_measurements,est_ground_truth,v_val,t_val, iter_param);
 [estGSTMRMM] = gen_estGSTMRMM_targets(est_measurements,est_ground_truth);
 [estVB_EOT_SN] = gen_estVB_EOT_SN_targets(est_measurements,est_ground_truth);
 [estVB_EOT_G] = gen_estVB_EOT_G_targets(est_measurements,est_ground_truth);
  for t=1:K
      mse_VB_EOT_NGHSSTM(1,t,expt)=(est_ground_truth.states(1,t)-estVB_EOT_NGHSSTM.states(1,t))^2+(est_ground_truth.states(2,t)-estVB_EOT_NGHSSTM.states(2,t))^2;
      wst_VB_EOT_NGHSSTM(1,t,expt)  = calc_wasserstein_sq(est_ground_truth.states([1,2],t), est_ground_truth.extents([1,2],[1,2],t), estVB_EOT_NGHSSTM.states([1,2],t), estVB_EOT_NGHSSTM.extents([1,2],[1,2],t));
     
      mse_GSTMRMM(1,t,expt)=(est_ground_truth.states(1,t)-estGSTMRMM.states(1,t))^2+(est_ground_truth.states(2,t)-estGSTMRMM.states(2,t))^2;
      wst_GSTMRMM(1,t,expt)  = calc_wasserstein_sq(est_ground_truth.states([1,2],t), est_ground_truth.extents([1,2],[1,2],t), estGSTMRMM.states([1,2],t), estGSTMRMM.extents([1,2],[1,2],t));
      
      mse_VB_EOT_SN(1,t,expt)=(est_ground_truth.states(1,t)-estVB_EOT_SN.states(1,t))^2+(est_ground_truth.states(2,t)-estVB_EOT_SN.states(2,t))^2;
      wst_VB_EOT_SN(1,t,expt)  = calc_wasserstein_sq(est_ground_truth.states([1,2],t), est_ground_truth.extents([1,2],[1,2],t), estVB_EOT_SN.states([1,2],t), estVB_EOT_SN.extents([1,2],[1,2],t));
        
      mse_VB_EOT_G(1,t,expt)=(est_ground_truth.states(1,t)-estVB_EOT_G.states(1,t))^2+(est_ground_truth.states(2,t)-estVB_EOT_G.states(2,t))^2;
      wst_VB_EOT_G(1,t,expt)  = calc_wasserstein_sq(est_ground_truth.states([1,2],t), est_ground_truth.extents([1,2],[1,2],t), estVB_EOT_G.states([1,2],t), estVB_EOT_G.extents([1,2],[1,2],t));
       
  end
end
rmse_VB_EOT_NGHSSTM=sqrt(mean(mse_VB_EOT_NGHSSTM,3));
rmse_GSTMRMM=sqrt(mean(mse_GSTMRMM,3));
rmse_VB_EOT_SN=sqrt(mean(mse_VB_EOT_SN,3));
rmse_VB_EOT_G=sqrt(mean(mse_VB_EOT_G,3));

%%%%%%%RMSE smooth
srmse_VB_EOT_NGHSSTM=smooth(rmse_VB_EOT_NGHSSTM,10);
srmse_GSTMRMM=smooth(rmse_GSTMRMM,10);
srmse_VB_EOT_SN=smooth(rmse_VB_EOT_SN,10);
srmse_VB_EOT_G=smooth(rmse_VB_EOT_G,10);

% %wst
avg_wasserstein_NGHSSTM = mean(wst_VB_EOT_NGHSSTM, 3);
avg_wasserstein_GSTM = mean(wst_GSTMRMM, 3);
avg_wasserstein_SN = mean(wst_VB_EOT_SN, 3);
avg_wasserstein_G = mean(wst_VB_EOT_G, 3);
% 平滑处理
smooth_wasserstein_NGHSSTM = smooth(avg_wasserstein_NGHSSTM, 10);
smooth_wasserstein_GSTM = smooth(avg_wasserstein_GSTM, 10);
smooth_wasserstein_SN = smooth(avg_wasserstein_SN, 10);
smooth_wasserstein_G = smooth(avg_wasserstein_G, 10);
%rmse
figure;
j = 1:t;
plot(j,srmse_VB_EOT_NGHSSTM,'-sg','MarkerFaceColor','g')
hold on;
plot(j,srmse_GSTMRMM,'-sb','MarkerFaceColor','b')
hold on;
 plot(j,srmse_VB_EOT_SN,'-sr','MarkerFaceColor','r')
  hold on;
  plot(j,srmse_VB_EOT_G,'-sy','MarkerFaceColor','y')
  hold on;
xlabel('Time step');
ylabel('RMSE (m)');
legend('VB-EOT-NGHSSTM','VB-EOT-GSTM','VB-EOT-SN','VB-EOT-G');
axis tight;
fprintf('NGHSSTM pos: %d\n',sum(srmse_VB_EOT_NGHSSTM)/K);
fprintf('GSTM: %d\n',sum(srmse_GSTMRMM)/K);
  fprintf('SN: %d\n',sum(srmse_VB_EOT_SN)/K);
  fprintf('G: %d\n',sum(srmse_VB_EOT_G)/K);


% %wst绘图
figure;
j = 1:K;
plot(j, smooth_wasserstein_NGHSSTM, '-sg', 'MarkerFaceColor','g');%MarkerSize','1.5', 'MarkerFaceColor','g'
hold on;
plot(j, smooth_wasserstein_GSTM, '-sb', 'MarkerFaceColor','b');
 plot(j, smooth_wasserstein_SN, '-sr', 'MarkerFaceColor','r');
 plot(j, smooth_wasserstein_G, '-sy', 'MarkerFaceColor','y');
xlabel('Time step');
ylabel('Gaussion Wasserstein Distance (m)');
legend('VB-EOT-NGHSSTM','GSTMRMM','VB-EOT-SN','VB-EOT-G');
axis tight;
fprintf('NGHSSTM: %d\n',sum(smooth_wasserstein_NGHSSTM)/K);
fprintf('GSTM: %d\n',sum(smooth_wasserstein_GSTM)/K);
fprintf('SN: %d\n',sum(smooth_wasserstein_SN)/K);
fprintf('G: %d\n',sum(smooth_wasserstein_G)/K);



snapshot_frames = [20, 70, 120]; 

figure('Position', [100, 100, 1200, 800], 'Color', 'w');


subplot(2, 3, [1, 2, 3]); 
hold on; grid on; box on; axis equal;
xlabel('X (m)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Y (m)', 'Interpreter', 'latex', 'FontSize', 12);
title('Overall Tracking Trajectory', 'Interpreter', 'latex', 'FontSize', 14);


plot(est_ground_truth.states(1,:), est_ground_truth.states(2,:), 'c-', 'LineWidth', 2.5);
plot(estVB_EOT_NGHSSTM.states(1,:), estVB_EOT_NGHSSTM.states(2,:), 'g-', 'LineWidth', 1.5);
plot(estGSTMRMM.states(1,:), estGSTMRMM.states(2,:), 'b-', 'LineWidth', 1.5);
plot(estVB_EOT_SN.states(1,:), estVB_EOT_SN.states(2,:), 'r-', 'LineWidth', 1.5);
plot(estVB_EOT_G.states(1,:), estVB_EOT_G.states(2,:), 'y-', 'LineWidth', 1.5);


plot_interval = 20; 
for k = 1:plot_interval:K
    plot_ellipse(est_ground_truth.extents(:,:,k), est_ground_truth.states(1:2, k), 1, 'Color', 'c', 'LineWidth', 1);
    plot_ellipse(estVB_EOT_NGHSSTM.extents(:,:,k), estVB_EOT_NGHSSTM.states(1:2, k), 1, 'Color', 'g', 'LineWidth', 1);
end


zoom_width = 80; % 虚线框宽度
zoom_height = 80; % 虚线框高度
for i = 1:length(snapshot_frames)
    k_snap = snapshot_frames(i);
    center_snap = est_ground_truth.states(1:2, k_snap);
    % 画虚线框
    rectangle('Position', [center_snap(1)-zoom_width/2, center_snap(2)-zoom_height/2, zoom_width, zoom_height], ...
              'EdgeColor', 'k', 'LineStyle', '-.', 'LineWidth', 1.5);
    % 标上序号
    text(center_snap(1)+zoom_width/2 + 5, center_snap(2), sprintf('Snap %d', i), ...
         'FontSize', 12, 'FontWeight', 'bold', 'Interpreter', 'latex');
end

% 主图图例
h_main = gobjects(5, 1);
h_main(1) = plot(NaN, NaN, 'c-', 'LineWidth', 2.5);
h_main(2) = plot(NaN, NaN, 'g-', 'LineWidth', 2);
h_main(3) = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);
h_main(4) = plot(NaN, NaN, 'r-', 'LineWidth', 1.5);
h_main(5) = plot(NaN, NaN, 'y-', 'LineWidth', 1.5);
legend(h_main, {'True Trajectory', 'VB-EOT-NGHSSTM (Proposed)', 'GSTMRMM', 'VB-EOT-SN', 'VB-EOT-G'}, ...
       'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 11);
set(gca, 'FontSize', 12);


for i = 1:length(snapshot_frames)
    k = snapshot_frames(i);
    subplot(2, 3, 3 + i); % 放在第 4, 5, 6 个位置
    hold on; grid on; box on; axis equal;
    
    Y_k = est_measurements{k};
    center = est_ground_truth.states(1:2, k);
    extent = est_ground_truth.extents(:,:,k);
    
    % 1. 绘制测量点
    scatter(Y_k(1, :), Y_k(2, :), 60, 'c', 'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
    
    % 2. 绘制所有算法的扩展椭圆
    plot_ellipse(extent, center, 1, 'Color', 'k', 'LineWidth', 2); % 真实轮廓用黑色，更清晰
    plot_ellipse(estGSTMRMM.extents(:,:,k), estGSTMRMM.states(1:2, k), 1, 'Color', 'b', 'LineWidth', 1.5); 
    plot_ellipse(estVB_EOT_SN.extents(:,:,k), estVB_EOT_SN.states(1:2, k), 1, 'Color', 'r', 'LineWidth', 1.5);
    plot_ellipse(estVB_EOT_G.extents(:,:,k), estVB_EOT_G.states(1:2, k), 1, 'Color', 'y', 'LineWidth', 1.5); 
    plot_ellipse(estVB_EOT_NGHSSTM.extents(:,:,k), estVB_EOT_NGHSSTM.states(1:2, k), 1, 'Color', 'g', 'LineWidth', 2.5, 'LineStyle', '--');
    
    % 3. 绘制所有算法的中心位置
    plot(center(1), center(2), 'kp', 'MarkerFaceColor','k', 'MarkerSize', 10); % 真实中心
    plot(estGSTMRMM.states(1, k), estGSTMRMM.states(2, k), 'bo', 'MarkerFaceColor','b', 'MarkerSize', 5);
    plot(estVB_EOT_SN.states(1, k), estVB_EOT_SN.states(2, k), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 5);
    plot(estVB_EOT_G.states(1, k), estVB_EOT_G.states(2, k), 'yo', 'MarkerFaceColor','y', 'MarkerSize', 5);
    plot(estVB_EOT_NGHSSTM.states(1, k), estVB_EOT_NGHSSTM.states(2, k), 'go', 'MarkerFaceColor','g', 'MarkerSize', 7); % 你的中心稍微大一点

    % 4. 动态限制坐标轴实现局部放大
    buffer = 40; % 边距缓冲，根据点云散布大小可微调
    xlim([center(1)-buffer, center(1)+buffer]);
    ylim([center(2)-buffer, center(2)+buffer]);
    
    % 设置子图标题和标签
    title(sprintf('Snapshot %d ($t = %d$)', i, k), 'Interpreter', 'latex', 'FontSize', 13);
    xlabel('X (m)', 'Interpreter', 'latex', 'FontSize', 12);
    if i == 1
        ylabel('Y (m)', 'Interpreter', 'latex', 'FontSize', 12);
    end
    set(gca, 'FontSize', 11);
end
function w2sq = calc_wasserstein_sq(true_mean, true_cov, est_mean, est_cov)
    % 确保协方差矩阵对称正定
    true_cov = (true_cov + true_cov')/2 + 1e-6*eye(size(true_cov));
    est_cov = (est_cov + est_cov')/2 + 1e-6*eye(size(est_cov));
    
    mean_diff = norm(true_mean - est_mean)^2;
    trace_term = trace(true_cov + est_cov - 2*sqrtm(sqrtm(true_cov)*est_cov*sqrtm(true_cov)));
    
    w2sq = mean_diff + trace_term;
end


