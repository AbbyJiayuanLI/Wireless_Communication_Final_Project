close all;
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For too many comparisons in the original plot, here we only
% implementent basic realization. 
% 
% Note that the running result of this file is not the figs in
% the paper, but it contains all the needed codes. If you want 
% to see the original figs, please manually uncomment the code,
% or refer to m files named as reproduce_paper_figx.
% 
%  2020.05.30 li jiayuan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
M = 500;
k = 20;
rho_d = 10;
var2_phase_b_range = [0:0.05:0.5];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. fixed apmlitude distribution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
at = 10^(-1/10);   
bt = 10^(1/10);
ar = 10^(-1/10);
br = 10^(1/10);
alpha_amp_b = 10^0; 
var2_amp_b = 0.5; 
pd2 = makedist('Normal','mu',alpha_amp_b,'sigma',sqrt(var2_amp_b));
pd_amp_t = truncate(pd2,at,bt);
pd_amp_r = truncate(pd2,ar,br);

% amplitude-error-related parameter after truncate
at_hat = (at-alpha_amp_b)/sqrt(var2_amp_b);
bt_hat = (bt-alpha_amp_b)/sqrt(var2_amp_b);
ar_hat = (ar-alpha_amp_b)/sqrt(var2_amp_b);
br_hat = (br-alpha_amp_b)/sqrt(var2_amp_b);
Zt = Phi_func(bt_hat) - Phi_func(at_hat);
Zr = Phi_func(br_hat) - Phi_func(ar_hat);
alpha_amp_t = alpha_amp_b + (phi_func(at_hat)-phi_func(bt_hat))*sqrt(var2_amp_b)/Zt;  
var2_amp_t = var2_amp_b*( 1 ...
                        +( at_hat*phi_func(at_hat)-bt_hat*phi_func(bt_hat) )/Zt ...
                        -( (phi_func(at_hat)-phi_func(bt_hat))/Zt )^2 );
alpha_amp_r = alpha_amp_b + (phi_func(ar_hat)-phi_func(br_hat))*sqrt(var2_amp_b)/Zr;  
var2_amp_r = var2_amp_b*( 1 ...
                        +( ar_hat*phi_func(ar_hat)-br_hat*phi_func(br_hat) )/Zr ...
                        -( (phi_func(ar_hat)-phi_func(br_hat))/Zr )^2 );   
At = alpha_amp_t^2 + var2_amp_t; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2.   SINR Calculation      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for idx =  1 : length(var2_phase_b_range)
    for idx2 = [0:2]                       % different phase mean
    % for idx3 = [1:4]                     % different bound for phase distribution
    
    
    % A. Phase Distribution
        % theta_b = 0;                                  %%%%%% uncomment when use idx3 %%%%%
        theta_b = idx2/18*3.14;
        var2_phase = var2_phase_b_range(idx);
        theta1 = -40/180*3.14;
        theta2 = 40/180*3.14;
        % theta1 = -idx3/18*3.14;                       %%%%%% uncomment when use idx3 %%%%%
        % theta2 = idx3/18*3.14;
        pd1 = makedist('Normal','mu',theta_b,'sigma',sqrt(var2_phase));
        pd_phase = truncate(pd1,theta1,theta2);

        
    % B. Phase-error-related parameter after truncate
        if idx == 1
            % µ¥¶À¼ÆËãµ±var2_phase=0
            gt = exp(1i*theta_b);
            gr = gt;
        else
            gt = exp( -var2_phase/2+1i*theta_b ) ...
                    *( erfz( (theta2-theta_b)/sqrt(2*var2_phase)-1i*sqrt(var2_phase/2) ) ...
                      -erfz( (theta1-theta_b)/sqrt(2*var2_phase)-1i*sqrt(var2_phase/2) ) ) ...
                    /( erf((theta2-theta_b)/sqrt(2*var2_phase)) - erf((theta1-theta_b)/sqrt(2*var2_phase)) );
            gr = gt;
        end
        AI = alpha_amp_t^2*alpha_amp_r^2 ...
                /(alpha_amp_t^2+var2_amp_t) ...
                /(alpha_amp_r^2+var2_amp_r) ...
                *abs(gt)^2*abs(gr)^2;
        
            
    % C. Analytical Result
        SINR_mrt = rho_d*(((M-1)*AI*At+2*At))*(k^2+rho_d*k*(k-1)*(rho_d*At^2+2*At))/(rho_d*(k-1)*At+k)^3;
        SINR_zf = rho_d*(M-k)*AI*At*(k^2+rho_d*k*(k-1)*At*(1-AI)*(rho_d*At*(1-AI)+2))/(rho_d*(k-1)*At*(1-AI)+k)^3;
        analytic_SINR_mrt(idx2+1,idx) = 10*log10(SINR_mrt);
        analytic_SINR_zf(idx2+1,idx) = 10*log10(SINR_zf);
        % analytic_SINR_mrt(idx3,idx) = 10*log10(SINR_mrt);     %%%%%% uncomment when use idx3 %%%%%
        % analytic_SINR_zf(idx3,idx) = 10*log10(SINR_zf);
    

    % D. Simulation Result
        % transmit signal
        s1_BPSK = round(unifrnd(0,1,[k,1]))*2-1;
        s1_BPSK(find(s1_BPSK>=0))=1/sqrt(2);
        s1_BPSK(find(s1_BPSK<0))=-1/sqrt(2);

        % channel matrix
        if idx == 1
            phase = ones(1,M)*theta_b;
        else 
            phase = random(pd_phase,[1,M]);
        end
        amp_t = random(pd_amp_t,[1,M]);
        amp_r = random(pd_amp_r,[1,M]);
        Hbt = diag(amp_t.*exp(1i*phase));
        Hbr = diag(amp_r.*exp(1i*phase));
        H = normrnd(0,0.5,[M,k])+1i*normrnd(0,0.5,[M,k]);
        Hu_hat = Hbr * H; % perfect estimation with tau=0
        Hd = H.'*Hbt;
        Hd_hat = Hu_hat.';

        % receive signal
        n = normrnd(0,0.5,[k,1])+1i*normrnd(0,0.5,[k,1]);
        
        % MRT
        w_mrt = conj(Hbr)*conj(H); 
        lambda_mrt = sqrt(1/trace(w_mrt*w_mrt'));
        y_mrt = sqrt(rho_d)*lambda_mrt*Hd*w_mrt*s1_BPSK+n;
        Pyk_mrt = abs(sqrt(rho_d)*lambda_mrt*H(:,1).'*Hbt*w_mrt(:,1)*s1_BPSK(1))^2;
        PI_mrt = 0;
        for idx4 = 2:k
            PI_mrt=PI_mrt+sqrt(rho_d)*lambda_mrt*H(:,1).'*Hbt*w_mrt(:,idx4)*s1_BPSK(idx4);
        end
        PI_mrt = abs(PI_mrt)^2;
        Pn_mrt = abs(n(1))^2;
        simu_SINR_mrt(idx2+1,idx) = 10*log10(Pyk_mrt/(PI_mrt+Pn_mrt));
        % simu_SINR_mrt(idx3,idx) = 10*log10(Pyk_mrt/(PI_mrt+Pn_mrt));      %%%%%% uncomment when use idx3 %%%%%

        % ZF
        w_zf = Hd_hat'/(Hd_hat*Hd_hat'); 
        lambda_zf = sqrt(1/trace(w_zf*w_zf'));
        y_zf = sqrt(rho_d)*lambda_zf*Hd*w_zf*s1_BPSK+n;
        Pyk_zf = abs(sqrt(rho_d)*lambda_zf*H(:,1).'*Hbt*w_zf(:,1)*s1_BPSK(1))^2;
        PI_zf = 0;
        for idx4 = 2:k
            PI_zf=PI_zf+sqrt(rho_d)*lambda_zf*H(:,1).'*Hbt*w_zf(:,idx4)*s1_BPSK(idx4);
        end
        PI_zf = abs(PI_zf)^2;
        Pn_zf = abs(n(1))^2;
        simu_SINR_zf(idx2+1,idx) = 10*log10(Pyk_zf/(PI_zf+Pn_zf));
        % simu_SINR_zf(idx3,idx) = 10*log10(Pyk_zf/(PI_zf+Pn_zf));          %%%%%% uncomment when use idx3 %%%%%

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         3. Plotting         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A. when use idx2 
figure;
plot(var2_phase_b_range,analytic_SINR_mrt,'r--*','linewidth',1.5); hold on;
% plot(var2_phase_b_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
% legend('QPSK','16QAM','QPSK with permutation','16QAM with permutation')
xlabel('Amplitude Error covairance')
ylabel('SINR(dB)')
title('Output SINR with MRT Precoding');

figure;
plot(var2_phase_b_range,analytic_SINR_zf,'r--*','linewidth',1.5); hold on;
% plot(var2_phase_b_range,simu_SINR_zf,'b--o','linewidth',1.5); hold on;
% legend({'SINR(Analytic)','SINR(simulated)'})
xlabel('Amplitude Error covairance')
ylabel('SINR(dB)')
title('Output SINR with ZF Precoding');

% % B. when use idx3
% figure;
% plot(var2_phase_b_range,analytic_SINR_mrt,'r--*','linewidth',1.5); hold on;
% % plot(var2_phase_b_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
% % legend('QPSK','16QAM','QPSK with permutation','16QAM with permutation')
% xlabel('Amplitude Error covairance')
% ylabel('SINR(dB)')
% title('Output SINR with MRT Precoding');
% 
% figure;
% plot(var2_phase_b_range,analytic_SINR_zf,'r--*','linewidth',1.5); hold on;
% % plot(var2_phase_b_range,simu_SINR_zf,'b--o','linewidth',1.5); hold on;
% legend({'SINR(Analytic)','SINR(simulated)'})
% xlabel('Amplitude Error covairance')
% ylabel('SINR(dB)')
% title('Output SINR with ZF Precoding');

toc;




% Defined Functions
function result = phi_func(num)
    result = sqrt(1/(2*pi))*exp(-num^2/2);
end

function result = Phi_func(num)
    result = ( 1+erf(num/sqrt(2)) )/2;
end






