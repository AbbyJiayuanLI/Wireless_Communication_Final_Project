close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the following simulation, perfect 
% reciprocity, normal level reciprocity error, and high 
% level reciprocity error is labeled seperately by 1, 2, 3.
% 
% Note that the running result of this file is not the figs in
% the paper, but it contains all the needed codes. If you want 
% to see the original figs, please manually uncomment the code,
% or refer to m files named as reproduce_paper_figx.
% 
% 2020.05.30 li jiayuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
k = 20;
num = 100;
% M_range = [50:50:500];
% rho_d = 10;
snr_range = [-10:10:40];
M = 500;

% for idx = 1:10
%     M = M_range(idx);
for idx = 1:6
    rho_d = 10^(snr_range(idx)/10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. fixed apmlitude distribution %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha_amp_b1 = 10^0; 
    var2_amp_b1 = 0; 

    at2 = 10^(-1/10);   
    bt2 = 10^(1/10);
    ar2 = 10^(-1/10);
    br2 = 10^(1/10);
    alpha_amp_b2 = 10^0; 
    var2_amp_b2 = 0.5; 
    pd2_2 = makedist('Normal','mu',alpha_amp_b2,'sigma',sqrt(var2_amp_b2));
    pd_amp_t2 = truncate(pd2_2,at2,bt2);
    pd_amp_r2 = truncate(pd2_2,ar2,br2);

    at3 = 10^(-4/10);   
    bt3 = 10^(4/10);
    ar3 = 10^(-4/10);
    br3 = 10^(4/10);
    alpha_amp_b3 = 10^0; 
    var2_amp_b3 = 1; 
    pd2_3 = makedist('Normal','mu',alpha_amp_b3,'sigma',sqrt(var2_amp_b3));
    pd_amp_t3 = truncate(pd2_3,at3,bt3);
    pd_amp_r3 = truncate(pd2_3,ar3,br3);

    % amplitude-error-related parameter after truncate
    alpha_amp_t1 = alpha_amp_b1;
    alpha_amp_r1 = alpha_amp_b1;
    var2_amp_t1 = var2_amp_b1;
    var2_amp_r1 = var2_amp_b1;
    At1 = alpha_amp_t1^2 + var2_amp_t1; 

    at_hat2 = (at2-alpha_amp_b2)/sqrt(var2_amp_b2);
    bt_hat2 = (bt2-alpha_amp_b2)/sqrt(var2_amp_b2);
    ar_hat2 = (ar2-alpha_amp_b2)/sqrt(var2_amp_b2);
    br_hat2 = (br2-alpha_amp_b2)/sqrt(var2_amp_b2);
    Zt2 = Phi_func(bt_hat2) - Phi_func(at_hat2);
    Zr2 = Phi_func(br_hat2) - Phi_func(ar_hat2);
    alpha_amp_t2 = alpha_amp_b2 + (phi_func(at_hat2)-phi_func(bt_hat2))*sqrt(var2_amp_b2)/Zt2;  
    var2_amp_t2 = var2_amp_b2*( 1 ...
                            +( at_hat2*phi_func(at_hat2)-bt_hat2*phi_func(bt_hat2) )/Zt2 ...
                            -( (phi_func(at_hat2)-phi_func(bt_hat2))/Zt2 )^2 );
    alpha_amp_r2 = alpha_amp_b2 + (phi_func(ar_hat2)-phi_func(br_hat2))*sqrt(var2_amp_b2)/Zr2;  
    var2_amp_r2 = var2_amp_b2*( 1 ...
                            +( ar_hat2*phi_func(ar_hat2)-br_hat2*phi_func(br_hat2) )/Zr2 ...
                            -( (phi_func(ar_hat2)-phi_func(br_hat2))/Zr2 )^2 );   
    At2 = alpha_amp_t2^2 + var2_amp_t2; 

    at_hat3 = (at3-alpha_amp_b3)/sqrt(var2_amp_b3);
    bt_hat3 = (bt3-alpha_amp_b3)/sqrt(var2_amp_b3);
    ar_hat3 = (ar3-alpha_amp_b3)/sqrt(var2_amp_b3);
    br_hat3 = (br3-alpha_amp_b3)/sqrt(var2_amp_b3);
    Zt3 = Phi_func(bt_hat3) - Phi_func(at_hat3);
    Zr3 = Phi_func(br_hat3) - Phi_func(ar_hat3);
    alpha_amp_t3 = alpha_amp_b3 + (phi_func(at_hat3)-phi_func(bt_hat3))*sqrt(var2_amp_b3)/Zt3;  
    var2_amp_t3 = var2_amp_b3*( 1 ...
                            +( at_hat3*phi_func(at_hat3)-bt_hat3*phi_func(bt_hat3) )/Zt3 ...
                            -( (phi_func(at_hat3)-phi_func(bt_hat3))/Zt3 )^2 );
    alpha_amp_r3 = alpha_amp_b3 + (phi_func(ar_hat3)-phi_func(br_hat3))*sqrt(var2_amp_b3)/Zr3;  
    var2_amp_r3 = var2_amp_b3*( 1 ...
                            +( ar_hat3*phi_func(ar_hat3)-br_hat3*phi_func(br_hat3) )/Zr3 ...
                            -( (phi_func(ar_hat3)-phi_func(br_hat3))/Zr3 )^2 );   
    At3 = alpha_amp_t3^2 + var2_amp_t3; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. fixed phase distribution %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta_b1 = 0;
    var2_phase1 = 0;

    theta_b2 = 0;
    var2_phase2 = 0.5;
    theta1_2 = -20/180*3.14;
    theta2_2 = 20/180*3.14;
    pd1_2 = makedist('Normal','mu',theta_b2,'sigma',sqrt(var2_phase2));
    pd_phase2 = truncate(pd1_2,theta1_2,theta2_2);

    theta_b3 = 0;
    var2_phase3 = 1;
    theta1_3 = -50/180*3.14;
    theta2_3 = 50/180*3.14;
    pd1_3 = makedist('Normal','mu',theta_b3,'sigma',sqrt(var2_phase3));
    pd_phase3 = truncate(pd1_3,theta1_3,theta2_3);

    % phase-error-related parameter after truncate
    gt1 = exp(1i*theta_b1);
    gr1 = gt1;
    AI1 = alpha_amp_t1^2*alpha_amp_r1^2 ...
                    /(alpha_amp_t1^2+var2_amp_t1) ...
                    /(alpha_amp_r1^2+var2_amp_r1) ...
                    *abs(gt1)^2*abs(gr1)^2;

    gt2 = exp( -var2_phase2/2+1i*theta_b2 ) ...
            *( erfz( (theta2_2-theta_b2)/sqrt(2*var2_phase2)-1i*sqrt(var2_phase2/2) ) ...
              -erfz( (theta1_2-theta_b2)/sqrt(2*var2_phase2)-1i*sqrt(var2_phase2/2) ) ) ...
            /( erf((theta2_2-theta_b2)/sqrt(2*var2_phase2)) - erf((theta1_2-theta_b2)/sqrt(2*var2_phase2)) );
    gr2 = gt2;
    AI2 = alpha_amp_t2^2*alpha_amp_r2^2 ...
                        /(alpha_amp_t2^2+var2_amp_t2) ...
                        /(alpha_amp_r2^2+var2_amp_r2) ...
                        *abs(gt2)^2*abs(gr2)^2;

    gt3 = exp( -var2_phase3/2+1i*theta_b3 ) ...
            *( erfz( (theta2_3-theta_b3)/sqrt(2*var2_phase3)-1i*sqrt(var2_phase3/2) ) ...
              -erfz( (theta1_3-theta_b3)/sqrt(2*var2_phase3)-1i*sqrt(var2_phase3/2) ) ) ...
            /( erf((theta2_3-theta_b3)/sqrt(2*var2_phase3)) - erf((theta1_3-theta_b3)/sqrt(2*var2_phase3)) );
    gr3 = gt3;
    AI3 = alpha_amp_t3^2*alpha_amp_r3^2 ...
                        /(alpha_amp_t3^2+var2_amp_t3) ...
                        /(alpha_amp_r3^2+var2_amp_r3) ...
                        *abs(gt3)^2*abs(gr3)^2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  2.   SINR Calculation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A. Analytical Result
    SINR_mrt1 = rho_d*(((M-1)*AI1*At1+2*At1))*(k^2+rho_d*k*(k-1)*(rho_d*At1^2+2*At1))/(rho_d*(k-1)*At1+k)^3;
    SINR_zf1 = rho_d*(M-k)*AI1*At1*(k^2+rho_d*k*(k-1)*At1*(1-AI1)*(rho_d*At1*(1-AI1)+2))/(rho_d*(k-1)*At1*(1-AI1)+k)^3;
    analytic_SINR_mrt1(idx) = 10*log10(SINR_mrt1);
    analytic_SINR_zf1(idx) = 10*log10(SINR_zf1);

    SINR_mrt2 = rho_d*(((M-1)*AI2*At2+2*At2))*(k^2+rho_d*k*(k-1)*(rho_d*At2^2+2*At2))/(rho_d*(k-1)*At2+k)^3;
    SINR_zf2 = rho_d*(M-k)*AI2*At2*(k^2+rho_d*k*(k-1)*At2*(1-AI2)*(rho_d*At2*(1-AI2)+2))/(rho_d*(k-1)*At2*(1-AI2)+k)^3;
    analytic_SINR_mrt2(idx) = 10*log10(SINR_mrt2);
    analytic_SINR_zf2(idx) = 10*log10(SINR_zf2);

    SINR_mrt3 = rho_d*(((M-1)*AI3*At3+2*At3))*(k^2+rho_d*k*(k-1)*(rho_d*At3^2+2*At3))/(rho_d*(k-1)*At3+k)^3;
    SINR_zf3 = rho_d*(M-k)*AI3*At3*(k^2+rho_d*k*(k-1)*At3*(1-AI3)*(rho_d*At3*(1-AI3)+2))/(rho_d*(k-1)*At3*(1-AI3)+k)^3;
    analytic_SINR_mrt3(idx) = 10*log10(SINR_mrt3);
    analytic_SINR_zf3(idx) = 10*log10(SINR_zf3);

    % B. Simulation Result
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         4. Plotting         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% plot(M_range,analytic_SINR_mrt1,'r--*','linewidth',1.5); hold on;
% plot(M_range,analytic_SINR_mrt2,'r--*','linewidth',1.5); hold on;
% plot(M_range,analytic_SINR_mrt3,'r--*','linewidth',1.5); hold on;
% % plot(M_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
% plot(M_range,analytic_SINR_zf1,'b--*','linewidth',1.5); hold on;
% plot(M_range,analytic_SINR_zf2,'b--*','linewidth',1.5); hold on;
% plot(M_range,analytic_SINR_zf3,'b--*','linewidth',1.5); hold on;
% % plot(M_range,simu_SINR_zf,'b--o','linewidth',1.5); hold on;
% % legend('QPSK','16QAM','QPSK with permutation','16QAM with permutation')
% xlabel('Amplitude Error covairance')
% ylabel('SINR(dB)')
% title('Output SINR over M');


% upper bound calculation
ub_mrt1 = M*AI1/k;
ub_mrt2 = M*AI2/k;
ub_mrt3 = M*AI3/k;

ub_zf1 = M*(AI1/(1-AI1))/k;
ub_zf2 = M*(AI2/(1-AI2))/k;
ub_zf3 = M*(AI3/(1-AI3))/k;

figure;
plot(snr_range,analytic_SINR_mrt1,'r--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_mrt2,'r--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_mrt3,'r--*','linewidth',1.5); hold on;
% plot(snr_range,ub_mrt1*ones(1,6),'r:','linewidth',1.5); hold on;
% plot(snr_range,ub_mrt2*ones(1,6),'r:','linewidth',1.5); hold on;
% plot(snr_range,ub_mrt3*ones(1,6),'r:','linewidth',1.5); hold on;
% plot(snr_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_zf1,'b--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_zf2,'b--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_zf3,'b--*','linewidth',1.5); hold on;
% plot(snr_range,ub_zf1*ones(1,6),'b:','linewidth',1.5); hold on;
% plot(snr_range,ub_zf2*ones(1,6),'b:','linewidth',1.5); hold on;
% plot(snr_range,ub_zf3*ones(1,6),'b:','linewidth',1.5); hold on;
% plot(snr_range,simu_SINR_zf,'b--o','linewidth',1.5); hold on;
legend('SINR_ MRT','SINR_ ZF')
ylim([0,25]);
xlabel('SNR \rho_d (dB)')
ylabel('SINR(dB)')
title('Output SINR over SNR');

toc;




% Defined Functions
function result = phi_func(num)
    result = sqrt(1/(2*pi))*exp(-num^2/2);
end

function result = Phi_func(num)
    result = ( 1+erf(num/sqrt(2)) )/2;
end









