close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implementation reproduces the fig.9 in paper.
% 
% 2020.05.30 li jiayuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
rho_d = 10^(20/10);
k = 20;
num = 100;
tau2 = 0.01;
var2_range = [0:0.05:1];
M = 500;

for idx = 1:21
    var2 = var2_range(idx);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. fixed apmlitude distribution %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    at4 = 10^(-4/10);   
    bt4 = 10^(4/10);
    ar4 = 10^(-4/10);
    br4 = 10^(4/10);
    alpha_amp_b4 = 10^0; 
    var2_amp_b4 = var2; 
    pd2_4 = makedist('Normal','mu',alpha_amp_b4,'sigma',sqrt(var2_amp_b4));
    pd_amp_t4 = truncate(pd2_4,at4,bt4);
    pd_amp_r4 = truncate(pd2_4,ar4,br4);

    % amplitude-error-related parameter after truncate
    if idx == 1
        % 单独计算当var2=0
        alpha_amp_t4 = alpha_amp_b4;
        alpha_amp_r4 = alpha_amp_b4;
        var2_amp_t4 = var2_amp_b4;
        var2_amp_r4 = var2_amp_b4;
    else
        at_hat4 = (at4-alpha_amp_b4)/sqrt(var2_amp_b4);
        bt_hat4 = (bt4-alpha_amp_b4)/sqrt(var2_amp_b4);
        ar_hat4 = (ar4-alpha_amp_b4)/sqrt(var2_amp_b4);
        br_hat4 = (br4-alpha_amp_b4)/sqrt(var2_amp_b4);
        Zt4 = Phi_func(bt_hat4) - Phi_func(at_hat4);
        Zr4 = Phi_func(br_hat4) - Phi_func(ar_hat4);
        alpha_amp_t4 = alpha_amp_b4 + (phi_func(at_hat4)-phi_func(bt_hat4))*sqrt(var2_amp_b4)/Zt4;  
        var2_amp_t4 = var2_amp_b4*( 1 ...
                                +( at_hat4*phi_func(at_hat4)-bt_hat4*phi_func(bt_hat4) )/Zt4 ...
                                -( (phi_func(at_hat4)-phi_func(bt_hat4))/Zt4 )^2 );
        alpha_amp_r4 = alpha_amp_b4 + (phi_func(ar_hat4)-phi_func(br_hat4))*sqrt(var2_amp_b4)/Zr4;  
        var2_amp_r4 = var2_amp_b4*( 1 ...
                                +( ar_hat4*phi_func(ar_hat4)-br_hat4*phi_func(br_hat4) )/Zr4 ...
                                -( (phi_func(ar_hat4)-phi_func(br_hat4))/Zr4 )^2 );   
    end
    At4 = alpha_amp_t4^2 + var2_amp_t4; 
    Ar4 = alpha_amp_r4^2 + var2_amp_r4; 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. fixed phase distribution %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    theta_b4 = 0;
    var2_phase4 = var2;
    theta1_4 = -50/180*3.14;
    theta2_4 = 50/180*3.14;
    pd1_4 = makedist('Normal','mu',theta_b4,'sigma',sqrt(var2_phase4));
    pd_phase4 = truncate(pd1_4,theta1_4,theta2_4);
    
    % phase-error-related parameter after truncate
    if idx == 1
        % 单独计算当var2_phase=0
        gt4 = exp(1i*theta_b4);
        gr4 = gt4;
    else
        gt4 = exp( -var2_phase4/2+1i*theta_b4 ) ...
                *( erfz( (theta2_4-theta_b4)/sqrt(2*var2_phase4)-1i*sqrt(var2_phase4/2) ) ...
                  -erfz( (theta1_4-theta_b4)/sqrt(2*var2_phase4)-1i*sqrt(var2_phase4/2) ) ) ...
                /( erf((theta2_4-theta_b4)/sqrt(2*var2_phase4)) - erf((theta1_4-theta_b4)/sqrt(2*var2_phase4)) );
        gr4 = gt4;
    end
    AI4 = alpha_amp_t4^2*alpha_amp_r4^2 ...
                        /(alpha_amp_t4^2+var2_amp_t4) ...
                        /(alpha_amp_r4^2+var2_amp_r4) ...
                        *abs(gt4)^2*abs(gr4)^2;
    BI4 = (1-tau2)*AI4*At4*Ar4/((1-tau2)*Ar4+tau2);
    BI4_hat = (1-tau2)*AI4;

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  2.   SINR Calculation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A. Analytical Result 
    SINR_mrt1 = rho_d*(((M-1)*AI4*At4+2*At4))*(k^2+rho_d*k*(k-1)*(rho_d*At4^2+2*At4))/(rho_d*(k-1)*At4+k)^3;
    SINR_zf1 = rho_d*(M-k)*AI4*At4*(k^2+rho_d*k*(k-1)*At4*(1-AI4)*(rho_d*At4*(1-AI4)+2))/(rho_d*(k-1)*At4*(1-AI4)+k)^3;
    ratio1(idx) = SINR_zf1/SINR_mrt1;
    
    SINR_mrt2 = rho_d*At4*( ((1-tau2)*Ar4*((M-1)*AI4+2)+tau2)/((1-tau2)*Ar4+tau2) ) ...
                *( (k^2+rho_d*k*(k-1)*(rho_d*At4^2+2*At4))/((rho_d*(k-1)*At4+k)^3) );
    SINR_zf2 = rho_d*(M-k)*BI4*(k^2+rho_d*k*(k-1)*(At4-BI4)*(rho_d*(At4-BI4)+2))/(rho_d*(k-1)*(At4-BI4)+k)^3;
    ratio2(idx) = SINR_zf2/SINR_mrt2;
    
    c1 = 1/(1-AI4);
    ratio3(idx) = c1;
    
    c2 = 1/(1-BI4_hat);
    ratio4(idx) = c2;
end
  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         4. Plotting         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
suptitle('Output SINR comparison of MRT and ZF');
subplot(1,2,1)
semilogy(var2_range,ratio1,'r--','linewidth',1.5); hold on;
semilogy(var2_range,ratio3,'b--','linewidth',1.5); hold on;
% plot(snr_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
legend('interpreter','latex',{'SINR_{ZF/MRT} (Analytic)','C_1^~ '})
ylim([10^0,10^2]);
xlabel('\rho_A^2 = \rho_D^2 (dB)')
ylabel('SINR Comparison Between ZF and MRT')
title('Reciprocity Error Only')

subplot(1,2,2)
semilogy(var2_range,ratio2,'r--','linewidth',1.5); hold on;
semilogy(var2_range,ratio4,'b--','linewidth',1.5); hold on;
% plot(snr_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
legend('interpreter','latex',{'SINR_{ZF/MRT}^{err} (Analytic)','C_1'})
ylim([10^0,10^2]);
xlabel('\rho_A^2 = \rho_D^2 (dB)')
ylabel('SINR Comparison Between ZF and MRT')
title('Compound Effect')

toc;




% Defined Functions
function result = phi_func(num)
    result = sqrt(1/(2*pi))*exp(-num^2/2);
end

function result = Phi_func(num)
    result = ( 1+erf(num/sqrt(2)) )/2;
end

