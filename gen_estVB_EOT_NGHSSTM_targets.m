function [est] = gen_estVB_EOT_NGHSSTM_targets(ett_measurements,ett_ground_truth, v_param ,t_param, iter_param)
K=150;
nz=2;
nx=4;
F = [eye(2) eye(2);zeros(2) eye(2)];
H =[eye(2) zeros(2)];
Q0 =0.1*[1/3*eye(2) 1/2*eye(2);1/2*eye(2) eye(2)];
R0=2*eye(2);
%xkk=ett_ground_truth.x00;  
xkk= [100;100;0;0]   ;            % [25;25;12;12]
Pkk=10*diag([10 10 1 1]) ;    %0.5*diag([10 10 1 1])
extentkk = [20 0; 0 20];
delta = 500;%500
A = eye(2)/sqrt(delta);
IWvkk=60;
IWVkk=extentkk*(IWvkk-2*nz-2);
for t=1:K 
    %预测
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
   
    %初始化
    sigama1=1;
    sigama2=1;
    q0=0.5;
    %hkk=0.5
    %ekk=0.5
    g0=0.5;
    E_extentkk=IWVk1k/(IWvk1k-2*nz-2);
    E_i_extentkk=IWvk1k*inv(IWVk1k);
    E_log_taok=psi(q0)-psi(1);
    E_log_1_taok=psi(q0)-psi(1);
    E_log_paik=psi(g0)-psi(1);
    E_log_1_paik=psi(g0)-psi(1);
    E_thitak=1;
    E_i_thitak=1/E_thitak;
    E_log_thitak=0;
    E_pk=1;
    E_i_pk=1/E_pk;
    E_log_pk=0;
    bk_bar=0.01*ones(nx,1);%0.001
    dk_bar=0.01*ones(nz,1);
    E_bk=bk_bar;
    P_bk=sigama1*eye(nx);
    E_dk=dk_bar;
    P_dk=sigama2*eye(nz);
    E_bkbk=E_bk*E_bk'+P_bk;
    E_dkdk=E_dk*E_dk'+P_dk;
    E_tk=t_param;
    E_gamak=t_param;
    c=0.25;
    nk=10;
    ZK=ett_measurements{t};
    zk_bar1=mean(ZK');
    zk_bar=zk_bar1';
    Zk_bar=(ZK-repmat(zk_bar,1,nk))*(ZK-repmat(zk_bar,1,nk))';
    N=10;
    w=v_param;
    v=v_param;
    for i=1:iter_param
        
        [V, D] = eig(E_extentkk);
        d = diag(D);
        d = max(d, 1e-10);  % 特征值截断
        inv_sqrt_E_extentkk = V * diag(1./sqrt(d)) * V';
        %避免不正定
        M = c * E_extentkk + R0;
        M = (M + M')/2;
        [V_m, D_m] = eig(M);
        d_m = diag(D_m);
        d_m = max(d_m, 1e-10);
        sqrt_M = V_m * diag(sqrt(d_m)) * V_m';
        B = sqrt_M * inv_sqrt_E_extentkk;
    
        %估计xkk，Pkk
        xkk_i=xk1k;
        Pk1k_hat=Pk1k/(E_tk+(1-E_tk)*E_i_thitak);
        R_hat=B*inv(E_i_extentkk)*B'/(nk*( E_gamak+(1- E_gamak)*E_i_pk));
        mk=((1-E_tk)*E_i_thitak*(xk1k-E_thitak*E_bk)+E_tk*xk1k)/(E_tk+(1-E_tk)*E_i_thitak);
        yitak=((1- E_gamak)* E_pk*E_i_pk*E_dk)/(E_gamak+(1- E_gamak)*E_i_pk);
        Kk=Pk1k_hat*H'*inv(H*Pk1k_hat*H'+R_hat);
        xkk=mk+Kk*(zk_bar-H*mk-yitak);
        Pkk=Pk1k_hat-Kk*H*Pk1k_hat;
        
        %估计bk和dk
        Lbk=(Pk1k* E_thitak)/(1-E_tk);
        Kbk=sigama1*E_thitak*inv((E_thitak^2*sigama1*eye(nx)+ Lbk));
        E_bk=bk_bar+Kbk*(xkk-xk1k-E_thitak*bk_bar);
        P_bk=sigama1*(eye(nx)-E_thitak*Kbk);
        
        
        Ldk=(E_pk* B*inv(E_i_extentkk)*B')/(nk*(1-E_gamak));
        Kdk=sigama2*E_pk*inv((E_pk^2*sigama2*eye(nz)+ Ldk));
        E_dk=dk_bar+Kdk*(zk_bar-H*xkk-E_pk*dk_bar);
        P_dk=sigama2*(eye(nz)-E_pk*Kdk);
       
        
        E_bkbk=E_bk*E_bk'+P_bk;
        E_dkdk=E_dk*E_dk'+P_dk;
        %估计扩展状态
      
        Ak=H*Pkk*H'+(zk_bar-H*xkk)*(zk_bar-H*xkk)';
        Bk=Ak-E_pk*E_dk*(zk_bar-H*xkk)'-E_pk*(zk_bar-H*xkk)*E_dk'+E_pk*E_pk*E_dkdk;
        IWvkk=IWvk1k+nk;
        IWVkk=IWVk1k+inv(B')*(E_gamak*Ak+(1-E_gamak)*Bk*E_i_pk)*nk*inv(B)+inv(B')*(E_gamak+(1-E_gamak)*E_i_pk)*Zk_bar*inv(B);
        E_extentkk=IWVkk/(IWvkk-2*nz-2);
        E_i_extentkk=IWvkk*inv(IWVkk);
        %估计tk和yk
        Ck=Pkk+(xkk-xk1k)*(xkk-xk1k)';
        Dk=Ck-E_thitak* E_bk*(xkk-xk1k)'-E_thitak*(xkk-xk1k)* E_bk'+E_thitak*E_thitak*E_bkbk;
        ptk1=exp( E_log_taok-0.5*trace(Ck*inv(Pk1k)));
        ptk0=exp(E_log_1_taok-0.5*nx* E_log_thitak-0.5*trace(Dk*E_i_thitak*inv(Pk1k)));
        if ptk1<=1e-100
            ptk1= 1e-100;
        end
        if ptk0<=1e-100
            ptk0= 1e-100;
        end
        E_tk= ptk1/(ptk1+ptk0);
        if E_tk<=1e-98
            E_tk= 1e-98;
        end
        if  E_tk>=0.9999
            E_tk= 0.9999;
        end
       
        
        
        pgamak1=exp(E_log_paik-0.5*trace(((Ak*nk)+Zk_bar)*inv(B*B')*E_i_extentkk));
        pgamak0=exp(E_log_1_paik-0.5*nz*nk*E_log_pk-0.5*trace(((Bk*nk)+Zk_bar)* E_i_pk*inv(B*B')*E_i_extentkk));
        if pgamak1<=1e-100
            pgamak1= 1e-100;
        end
        if pgamak0<=1e-100
            pgamak0= 1e-100;
        end
        E_gamak=pgamak1/(pgamak1+pgamak0);
        if  E_gamak<=1e-98
            E_gamak= 1e-98;
        end
        if  E_gamak>=0.99999
           E_gamak=0.999999;
        end
        
        
        %估计thitak和pk
        Y1=nx*(1-E_tk);
        Y2=(1-E_tk)*trace(Ck*inv(Pk1k));
        Y3=(1-E_tk)*trace(E_bkbk*inv(Pk1k));

        E_thitak=(-(Y1+w+2)+sqrt((Y1+w+2)^2+4*Y3*(Y2+w)))/(2*Y3);
        E_i_thitak=1/E_thitak;
        E_log_thitak=log(E_thitak);
        
        Y4=nz*nk*(1-E_gamak);
        Y5=(1-E_gamak)*trace(inv(B') * E_i_extentkk * inv(B) * ((Ak*nk) + Zk_bar));
        Y6=(1-E_gamak)*trace(inv(B') * E_i_extentkk * inv(B) * E_dkdk * nk);

        E_pk=(-(Y4+v+2)+sqrt((Y4+v+2)^2+4*Y6*(Y5+v)))/(2*Y6);
        E_i_pk=1/E_pk;
        E_log_pk=log(E_pk);
        
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
    est.states = eststates;
    est.extents = estextents;
end


end
