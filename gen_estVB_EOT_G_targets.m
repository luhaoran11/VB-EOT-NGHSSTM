function [estVB_EOT_G] = gen_estVB_EOT_G_targets(ett_measurements,ett_ground_truth)
K=150;
nz=2;
nx=4;
F = [eye(2) eye(2);zeros(2) eye(2)];
H =[eye(2) zeros(2)];
Q0 =0.1*[1/3*eye(2) 1/2*eye(2);1/2*eye(2) eye(2)];
R0=2*eye(2);
xkk= [100;100;0;0] ;     % [25;25;12;12]
Pkk=10*diag([10 10 1 1]) ;    %0.5*diag([10 10 1 1])
extentkk = [20 0; 0 20];
delta =500;
IWvkk=60;%60
IWVkk=extentkk*(IWvkk-2*nz-2);
A = eye(2)/sqrt(delta);
dk=[0  0 ]';
kasaikk=1;
sitakk2=0.01;
k=0.01;
c=0.25;
for t=1:K 

    xk1k=F*xkk;
    Pk1k=F*Pkk*F'+Q0;
    lambda=IWvkk-2*nz-2;
    IWvk1k=2*delta*(lambda+1)*(lambda-1)*(lambda-2)/(lambda*lambda*(lambda+delta))+2*nz+4;
    IWVk1k=delta*(IWvk1k-2*nz-2)*A*IWVkk*(A)'/(lambda-1);
    
    %避免不正定
    IWVk1k = (IWVk1k + IWVk1k')/2;
    min_eig = min(eig(IWVk1k));
    if min_eig <= 0
        IWVk1k = IWVk1k + (abs(min_eig) + 1e-8) * eye(size(IWVk1k));
    end
    kasaik1k=kasaikk;
    sitak1k2=sitakk2+k^2;
    %初始化
    xkk=xk1k;
    Pkk=Pk1k;
    IWvkk=IWvk1k;
    IWVkk=IWVk1k;
    nk=10;
    ZK=ett_measurements{t};
    zk_bar1=mean(ZK');
    zk_bar=zk_bar1';
    Zk_bar=(ZK-repmat(zk_bar,1,nk))*(ZK-repmat(zk_bar,1,nk))';
    N=10;
    E_extentkk=IWVk1k/(IWvk1k-2*nz-2);

    for i=1:10;
        
         [V, D] = eig(E_extentkk);
    d = diag(D);
    d = max(d, 1e-10);  % 特征值截断
    inv_sqrt_E_extentkk = V * diag(1./sqrt(d)) * V';
    M = c * E_extentkk + R0;
    M = (M + M')/2;
    [V_m, D_m] = eig(M);
    d_m = diag(D_m);
    d_m = max(d_m, 1e-10);
    sqrt_M = V_m * diag(sqrt(d_m)) * V_m';
    B = sqrt_M * inv_sqrt_E_extentkk;
    
        deltak=nk*sitak1k2*dk'*inv(B*IWVkk*B')*dk;
        kasaikk=(kasaik1k+deltak*(zk_bar-H*xkk)'*(zk_bar-H*xkk))/(deltak+1);
        sitakk2= sitak1k2/(deltak+1);   
        taok= normcdf(-kasaikk/sqrt(sitakk2))/(1- normcdf(-kasaikk/sqrt(sitakk2)));
        VK=0;
        E_pk=kasaikk+sqrt(sitakk2)*taok;
        IWvkk=IWvk1k+nk;
        Dk=(zk_bar-H*xkk-dk*E_pk)*(zk_bar-H*xkk-dk*E_pk)'+H*Pkk*H';
        IWVkk=IWVk1k+inv(B)*(Zk_bar+nk*Dk)*inv(B');
        E_extentkk=IWVkk/(IWvkk-2*nz-2);
        Sk=H*Pk1k*H'+(B*E_extentkk*B')/nk;
        Kk=Pk1k*H'*inv(Sk);
        Pkk=Pk1k-Kk*H*Pk1k;
        xkk=xk1k+Kk*(zk_bar-H*xk1k-dk*E_pk);
    end
     estextents(:, :, t)=E_extentkk;
     eststates(:, t)= xkk;
     estVB_EOT_G.states = eststates;
     estVB_EOT_G.extents = estextents;
end
end