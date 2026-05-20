%生成真实值和测量值

function [ett_measurements,ett_ground_truth] = gen_data_targets()
profile on;
clear all;
close all;
clc;
K = 150;

% Generate model parameters, ground truth and measurements 
x_dim = 4;
z_dim = 2;
model.x_dim = x_dim;                           % dimension of target state
model.z_dim = z_dim;                           % dimension of measurements
                                       %？ dimension of target extent目标范围尺寸
model.K = K;                                   % number of scans扫描次数
model.c = 0.25;
model.delta = 500;
model.Ak = eye(z_dim)/sqrt(model.delta);
% motion model parameters for Gaussian distribution
Ts = 1;
model.Ts = Ts;
model.F = [eye(2) Ts*eye(2);zeros(2) eye(2)];  % transition model（目标做匀速运动）
                                % process noise covariance         
model.Q0 = 0.1*[1/3*eye(2) 1/2*eye(2);1/2*eye(2) eye(2)];%用克罗内克积求过程噪声协方差

% measurement model
model.H = [1 0 0 0;0 1 0 0];%提取位置 前两维是位置
% sigma_v = 1;
% model.D= eye(z_dim)*sigma_v^2; 
% model.R0= model.D*model.D';                     % observation noise covariance
model.R0= 2*eye(2);
x=[0;0;10;10];

extent = [20 0; 0 20];       % True extent (elliptical shape matrix)

states = zeros(4, K);  % Ground truth states [x, y, vx, vy, theta, omega]
extents = zeros(2, 2, K); % Ground truth extents
Y_k = cell(1, K);

nz=10;
c=0.25;
nu=5;


% Vkk=extent/model.delta


%%%%KF (Kalman filter)   KFTNCM

for t=1:K
        
        %%%%%%%Simulate non-stationary noise
        %%%%%%%Gaussian noise
        if t<=50
            p1=0.8;
            p2=0.8;
        end
        %%%%%%%Slightly heavy-tailed noise
        if (t>50)&&(t<=100)
            p1=0.8;
            p2=0.8;
        end
        %%%%%%%Moderately heavy-tailed noise
        if (t>100)&&(t<=150)
            p1=0.8;
            p2=0.8;
        end
        
        
        %%%%Simulate true state and measurement
        test1=rand;
        test2=rand;
        if test1<=p1
            Q=model.Q0;
            wk=utchol(Q)*randn(x_dim,1);%%高斯噪声
        else
            Q=10*model.Q0;
            wk=utchol(Q)*randn(x_dim,1);
           % wk=skew_noise(diag([1.0;1.0;1.0;1.0]),model.Q0,nu);
            %             if rand<1
            %             wk=skew_noise(diag([1.0;1.0;1.0;1.0]),Q0,5);
            %             end
        end
        x=model.F*x+wk;
        
        extent=wishrnd(model.Ak*extent*model.Ak, model.delta);
  %extent=extent;
        extent = (extent + extent')/2; % 强制对称
        [V,D] = eig(extent);
        d = diag(D);
        if any(d <= 0)
            % 替换非正特征值
            d(d <= 0) = 1e-10;
            % 调整对角元素保持均值不变
            diag_adjustment = trace(D) - sum(d);
            d = d + diag_adjustment/(p*10);
            extent = V*diag(d)*V';
        end
       

        
        if test2<=p2 %%
            R=model.R0;
            [V, D] = eig(extent);
            d = diag(D);
            d = max(d, 1e-10);  % 特征值截断
            inv_sqrt_extent = V * diag(1./sqrt(d)) * V';
            %避免不正定
            M = c * extent + R;
            M = (M + M')/2;
            [V_m, D_m] = eig(M);
            d_m = diag(D_m);
            d_m = max(d_m, 1e-10);
            sqrt_M = V_m * diag(sqrt(d_m)) * V_m';
            B = sqrt_M * inv_sqrt_extent;
        
          
            Y_k{t} = mvnrnd([0, 0], B*extent*B', nz)' +  model.H * x;
        else
           
            n=size(model.R0,1);         
            Lambda=zeros(n);            
            r1=20;            
           
            for i=1:n
                Lambda(i,i)=gamrnd(nu/2,2/nu);
            end
                    
            inv_Lambda=inv(Lambda);
            R=10*model.R0;%1
            [V, D] = eig(extent);
            d = diag(D);
            d = max(d, 1e-10);  % 特征值截断
            inv_sqrt_extent = V * diag(1./sqrt(d)) * V';
            %避免不正定
            M = c * extent + R;
            M = (M + M')/2;
            [V_m, D_m] = eig(M);
            d_m = diag(D_m);
            d_m = max(d_m, 1e-10);
            sqrt_M = V_m * diag(sqrt(d_m)) * V_m';
            B = sqrt_M * inv_sqrt_extent;
            R0=B*extent*B';
            R0 = (R0 + R0')/2;
            min_eig = min(eig(R0));
            if min_eig <= 0
                R0 = R0 + (abs(min_eig) + 1e-8) * eye(size(R0));
            end
            
            u=abs(utchol(r1*inv_Lambda)*randn(n,1)); 
%             %P= ((extent*Lambda) + (extent*Lambda)')/2;
%             min_eig = min(eig( P));
%             if min_eig <= 0
%                 P =  P + (abs(min_eig) + 1e-8) * eye(size( P));
%             end
            Y_k{t} =diag([1.0;1.0])*u+utchol(R0*inv_Lambda)*randn(n,nz)+  model.H * x;
            %diag([1.0;1.0])*
        end
 
        
%         mu = (abs(utchol(d.*inv_Lambda))*randn(z_dim,1))'
      
    
        extents(:, :, t)=extent;
        states(:, t)= x;
        % Plotting the current scenario:
    ett_measurements = Y_k;
    ett_ground_truth.states = states;
    ett_ground_truth.extents = extents;
   
end

function X = iwishrnd(Psi, nu)
% 生成服从逆威沙特分布的随机矩阵
% 输入：
%   Psi - 尺度矩阵 (p×p 对称正定)
%   nu - 自由度 (nu > p-1)
% 输出：
%   X - 服从逆威沙特分布的随机矩阵

p = size(Psi, 1);
validateattributes(Psi, {'double'}, {'size', [p p], 'symmetric'});
if nu <= p - 1
    error('自由度nu必须大于p-1 (当前p=%d, nu=%d)', p, nu);
end

% 通过威沙特分布生成逆矩阵
% 首先生成威沙特分布随机矩阵: W ~ Wishart(inv(Psi), nu)
W = wishrnd(inv(Psi), nu);

% 取逆得到逆威沙特分布
X = inv(W);

% 数值稳定性处理
[~, p] = chol(X, 'lower');
if p ~= 0
    % 添加正则化项确保正定性
    X = X + eye(size(X)) * 1e-6 * norm(X, 'fro');
end



