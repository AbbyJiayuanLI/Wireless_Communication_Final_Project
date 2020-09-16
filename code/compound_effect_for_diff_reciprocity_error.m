close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This implementation reproduces the fig.8 in paper.
% 
% 2020.05.30 li jiayuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
k = 20;
num = 100;
tau2 = 0.1;
snr_range = [-10:5:40];
M = 500;

for idx = 1:11
    rho_d = 10^(snr_range(idx)/10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. fixed apmlitude distribution %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alpha_amp_b1 = 10^0; 
    var2_amp_b1 = 0; 

    at3 = 10^(-1/10);   
    bt3 = 10^(1/10);
    ar3 = 10^(-1/10);
    br3 = 10^(1/10);
    alpha_amp_b3 = 10^0; 
    var2_amp_b3 = 0.5; 
    pd2_3 = makedist('Normal','mu',alpha_amp_b3,'sigma',sqrt(var2_amp_b3));
    pd_amp_t3 = truncate(pd2_3,at3,bt3);
    pd_amp_r3 = truncate(pd2_3,ar3,br3);

    at4 = 10^(-4/10);   
    bt4 = 10^(4/10);
    ar4 = 10^(-4/10);
    br4 = 10^(4/10);
    alpha_amp_b4 = 10^0; 
    var2_amp_b4 = 1; 
    pd2_4 = makedist('Normal','mu',alpha_amp_b4,'sigma',sqrt(var2_amp_b4));
    pd_amp_t4 = truncate(pd2_4,at4,bt4);
    pd_amp_r4 = truncate(pd2_4,ar4,br4);

    % amplitude-error-related parameter after truncate
    alpha_amp_t1 = alpha_amp_b1;
    alpha_amp_r1 = alpha_amp_b1;
    var2_amp_t1 = var2_amp_b1;
    var2_amp_r1 = var2_amp_b1;
    At1 = alpha_amp_t1^2 + var2_amp_t1; 
    Ar1 = alpha_amp_r1^2 + var2_amp_r1; 
    
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
    Ar3 = alpha_amp_r3^2 + var2_amp_r3; 

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
    At4 = alpha_amp_t4^2 + var2_amp_t4; 
    Ar4 = alpha_amp_r4^2 + var2_amp_r4; 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. fixed phase distribution %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta_b1 = 0;
    var2_phase1 = 0;

    theta_b3 = 0;
    var2_phase3 = 0.5;
    theta1_3 = -20/180*3.14;
    theta2_3 = 20/180*3.14;
    pd1_3 = makedist('Normal','mu',theta_b3,'sigma',sqrt(var2_phase3));
    pd_phase3 = truncate(pd1_3,theta1_3,theta2_3);

    theta_b4 = 0;
    var2_phase4 = 1;
    theta1_4 = -50/180*3.14;
    theta2_4 = 50/180*3.14;
    pd1_4 = makedist('Normal','mu',theta_b4,'sigma',sqrt(var2_phase4));
    pd_phase4 = truncate(pd1_4,theta1_4,theta2_4);

    % phase-error-related parameter after truncate
    gt1 = exp(1i*theta_b1);
    gr1 = gt1;
    AI1 = alpha_amp_t1^2*alpha_amp_r1^2 ...
                    /(alpha_amp_t1^2+var2_amp_t1) ...
                    /(alpha_amp_r1^2+var2_amp_r1) ...
                    *abs(gt1)^2*abs(gr1)^2;
    BI1 = AI1*At1*Ar1/Ar1;
    
    BI2 = (1-tau2)*AI1*At1*Ar1/((1-tau2)*Ar1+tau2);
    
    gt3 = exp( -var2_phase3/2+1i*theta_b3 ) ...
            *( erfz( (theta2_3-theta_b3)/sqrt(2*var2_phase3)-1i*sqrt(var2_phase3/2) ) ...
              -erfz( (theta1_3-theta_b3)/sqrt(2*var2_phase3)-1i*sqrt(var2_phase3/2) ) ) ...
            /( erf((theta2_3-theta_b3)/sqrt(2*var2_phase3)) - erf((theta1_3-theta_b3)/sqrt(2*var2_phase3)) );
    gr3 = gt3;
    AI3 = alpha_amp_t3^2*alpha_amp_r3^2 ...
                        /(alpha_amp_t3^2+var2_amp_t3) ...
                        /(alpha_amp_r3^2+var2_amp_r3) ...
                        *abs(gt3)^2*abs(gr3)^2;
    BI3 = (1-tau2)*AI3*At3*Ar3/((1-tau2)*Ar3+tau2);
    
    gt4 = exp( -var2_phase4/2+1i*theta_b4 ) ...
            *( erfz( (theta2_4-theta_b4)/sqrt(2*var2_phase4)-1i*sqrt(var2_phase4/2) ) ...
              -erfz( (theta1_4-theta_b4)/sqrt(2*var2_phase4)-1i*sqrt(var2_phase4/2) ) ) ...
            /( erf((theta2_4-theta_b4)/sqrt(2*var2_phase4)) - erf((theta1_4-theta_b4)/sqrt(2*var2_phase4)) );
    gr4 = gt4;
    AI4 = alpha_amp_t4^2*alpha_amp_r4^2 ...
                        /(alpha_amp_t4^2+var2_amp_t4) ...
                        /(alpha_amp_r4^2+var2_amp_r4) ...
                        *abs(gt4)^2*abs(gr4)^2;
    BI4 = (1-tau2)*AI4*At4*Ar4/((1-tau2)*Ar4+tau2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  2.   SINR Calculation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A. Analytical Result
    SINR_mrt1 = rho_d*(((M-1)*AI1*At1+2*At1))*(k^2+rho_d*k*(k-1)*(rho_d*At1^2+2*At1))/(rho_d*(k-1)*At1+k)^3;
    SINR_zf1 = rho_d*(M-k)*AI1*At1*(k^2+rho_d*k*(k-1)*At1*(1-AI1)*(rho_d*At1*(1-AI1)+2))/(rho_d*(k-1)*At1*(1-AI1)+k)^3;
    analytic_SINR_mrt1(idx) = 10*log10(SINR_mrt1);
    analytic_SINR_zf1(idx) = 10*log10(SINR_zf1);
    
    SINR_mrt2 = rho_d*At1*( ((1-tau2)*Ar1*((M-1)*AI1+2)+tau2)/((1-tau2)*Ar1+tau2) ) ...
                *( (k^2+rho_d*k*(k-1)*(rho_d*At1^2+2*At1))/((rho_d*(k-1)*At1+k)^3) );
    SINR_zf2 = rho_d*(M-k)*BI2*(k^2+rho_d*k*(k-1)*(At1-BI2)*(rho_d*(At1-BI2)+2))/(rho_d*(k-1)*(At1-BI2)+k)^3;
    analytic_SINR_mrt2(idx) = 10*log10(SINR_mrt2);
    analytic_SINR_zf2(idx) = 10*log10(SINR_zf2);
    
    SINR_mrt3 = rho_d*At3*( ((1-tau2)*Ar3*((M-1)*AI3+2)+tau2)/((1-tau2)*Ar3+tau2) ) ...
                *( (k^2+rho_d*k*(k-1)*(rho_d*At3^2+2*At3))/((rho_d*(k-1)*At3+k)^3) );
    SINR_zf3 = rho_d*(M-k)*BI3*(k^2+rho_d*k*(k-1)*(At3-BI3)*(rho_d*(At3-BI3)+2))/(rho_d*(k-1)*(At3-BI3)+k)^3;
    analytic_SINR_mrt3(idx) = 10*log10(SINR_mrt3);
    analytic_SINR_zf3(idx) = 10*log10(SINR_zf3);

    SINR_mrt4 = rho_d*At4*( ((1-tau2)*Ar4*((M-1)*AI4+2)+tau2)/((1-tau2)*Ar4+tau2) ) ...
                *( (k^2+rho_d*k*(k-1)*(rho_d*At4^2+2*At4))/((rho_d*(k-1)*At4+k)^3) );
    SINR_zf4 = rho_d*(M-k)*BI4*(k^2+rho_d*k*(k-1)*(At4-BI4)*(rho_d*(At4-BI4)+2))/(rho_d*(k-1)*(At4-BI4)+k)^3;
    analytic_SINR_mrt4(idx) = 10*log10(SINR_mrt4);
    analytic_SINR_zf4(idx) = 10*log10(SINR_zf4);

    % B. Simulation Result
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
ub_mrt1 = M*((rho_d*BI1)/(rho_d+1/At1))/k;
ub_mrt2 = M*((rho_d*BI2)/(rho_d+1/At1))/k;
ub_mrt3 = M*((rho_d*BI3)/(rho_d+1/At3))/k;
ub_mrt4 = M*((rho_d*BI4)/(rho_d+1/At4))/k;

ub_zf1 = (M-k)*((rho_d*BI1)/(rho_d*(1-BI1)+1/At1))/k;
ub_zf2 = (M-k)*((rho_d*BI2)/(rho_d*(1-BI2)+1/At1))/k;
ub_zf3 = (M-k)*((rho_d*BI3)/(rho_d*(1-BI3)+1/At3))/k;
ub_zf4 = (M-k)*((rho_d*BI4)/(rho_d*(1-BI4)+1/At4))/k;

figure;
subplot(1,2,1)
plot(snr_range,analytic_SINR_mrt1,'r--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_mrt2,'r--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_mrt3,'r--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_mrt4,'r--*','linewidth',1.5); hold on;
plot(snr_range,ub_mrt1*ones(1,11),'r:','linewidth',1.5); hold on;
plot(snr_range,ub_mrt2*ones(1,11),'r:','linewidth',1.5); hold on;
plot(snr_range,ub_mrt3*ones(1,11),'r:','linewidth',1.5); hold on;
plot(snr_range,ub_mrt4*ones(1,11),'r:','linewidth',1.5); hold on;
% plot(snr_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
% legend('QPSK','16QAM','QPSK with permutation','16QAM with permutation')
ylim([2,16]);
xlabel('\rho_d (dB)')
ylabel('SINR(dB)')
title('MRT Output SINR over SNR');

subplot(1,2,2)
plot(snr_range,analytic_SINR_zf1,'b--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_zf2,'b--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_zf3,'b--*','linewidth',1.5); hold on;
plot(snr_range,analytic_SINR_zf4,'b--*','linewidth',1.5); hold on;
plot(snr_range,ub_zf1*ones(1,11),'b:','linewidth',1.5); hold on;
plot(snr_range,ub_zf2*ones(1,11),'b:','linewidth',1.5); hold on;
plot(snr_range,ub_zf3*ones(1,11),'b:','linewidth',1.5); hold on;
plot(snr_range,ub_zf4*ones(1,11),'b:','linewidth',1.5); hold on;
% plot(snr_range,simu_SINR_zf,'b--o','linewidth',1.5); hold on;
% legend('QPSK','16QAM','QPSK with permutation','16QAM with permutation')
ylim([0,25]);
xlabel('\rho_d (dB)')
ylabel('SINR(dB)')
title('ZF Output SINR over SNR');

toc;




% Defined Functions
function result = phi_func(num)
    result = sqrt(1/(2*pi))*exp(-num^2/2);
end

function result = Phi_func(num)
    result = ( 1+erf(num/sqrt(2)) )/2;
end









