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
num = 1;  % number of repeated trials


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. fixed phase distribution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_b = 0;
var2_phase = 0.5;
theta1 = -20/180*3.14;
theta2 = 20/180*3.14;
pd1 = makedist('Normal','mu',theta_b,'sigma',sqrt(var2_phase));
pd_phase = truncate(pd1,theta1,theta2);

% phase-error-related parameter after truncate
gt = exp( -var2_phase/2+1i*theta_b ) ...
        *( erfz( (theta2-theta_b)/sqrt(2*var2_phase)-1i*sqrt(var2_phase/2) ) ...
          -erfz( (theta1-theta_b)/sqrt(2*var2_phase)-1i*sqrt(var2_phase/2) ) ) ...
        /( erf((theta2-theta_b)/sqrt(2*var2_phase)) - erf((theta1-theta_b)/sqrt(2*var2_phase)) );
gr = gt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2.   SINR Calculation      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at = 10^(-1/10);   
% bt = 10^(1/10);
ar = 10^(-1/10);
br = 10^(1/10);
var2_amp_b_range = [0:0.1:0.5]; 
alpha_amp_b = 10^0; 


for idx =  1 : length(var2_amp_b_range)
    for idx2 = [1:4]  
        
        % A. Amplitude Distribution
            at = 10^(-idx2/10);   
            bt = 10^(idx2/10);
            % ar = 10^(-idx2/10);
            % br = 10^(idx2/10);
            var2_amp_b = var2_amp_b_range(idx); 
            pd2 = makedist('Normal','mu',alpha_amp_b,'sigma',sqrt(var2_amp_b));
            pd_amp_t = truncate(pd2,at,bt);
            pd_amp_r = truncate(pd2,ar,br);


        % B. Amplitude-error-related parameter after truncate
            if idx == 1
                % 单独计算当var2_amp_b=0
                alpha_amp_t = alpha_amp_b;
                alpha_amp_r = alpha_amp_b;
                var2_amp_t = var2_amp_b;
                var2_amp_r = var2_amp_b;
            else
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

    %%%%%%%%% 或者直接从truncate函数里拿均值方差 %%%%%%%%%% 
    %         alpha_amp_t = mean(pd_amp_t);
    %         alpha_amp_r = mean(pd_amp_r);
    %         var2_amp_t = var(pd_amp_t); 
    %         var2_amp_r = var(pd_amp_r);

            end
            At = alpha_amp_t^2 + var2_amp_t; 
            AI = alpha_amp_t^2*alpha_amp_r^2 ...
                    /(alpha_amp_t^2+var2_amp_t) ...
                    /(alpha_amp_r^2+var2_amp_r) ...
                    *abs(gt)^2*abs(gr)^2;


        % C. Analytical Result
            SINR_mrt = rho_d*(((M-1)*AI*At+2*At))*(k^2+rho_d*k*(k-1)*(rho_d*At^2+2*At))/(rho_d*(k-1)*At+k)^3;
            SINR_zf = rho_d*(M-k)*AI*At*(k^2+rho_d*k*(k-1)*At*(1-AI)*(rho_d*At*(1-AI)+2))/(rho_d*(k-1)*At*(1-AI)+k)^3;
            analytic_SINR_mrt(idx2,idx) = 10*log10(SINR_mrt);
            analytic_SINR_zf(idx2,idx) = 10*log10(SINR_zf);


        % D. Simulation Result
        sinr1 = 0;
        sinr2 = 0;
        for idx3 = [1:num]
            % transmit signal
            s1_BPSK = round(unifrnd(0,1,[k,1]))*2-1;
            s1_BPSK(find(s1_BPSK>=0))=1/sqrt(2);
            s1_BPSK(find(s1_BPSK<0))=-1/sqrt(2);

            % channel matrix
            phase = random(pd_phase,[1,M]);
            if idx == 1
                amp_t = ones(1,M);
                amp_r = ones(1,M);
            else
                amp_t = random(pd_amp_t,[1,M]);
                amp_r = random(pd_amp_r,[1,M]);
            end
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
            sinr1 = sinr1+10*log10(Pyk_mrt/(PI_mrt+Pn_mrt));
%             simu_SINR_mrt(idx2,idx) = 10*log10(Pyk_mrt/(PI_mrt+Pn_mrt));
            
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
            sinr2 = sinr2+10*log10(Pyk_zf/(PI_zf+Pn_zf));
%             simu_SINR_zf(idx2,idx) = 10*log10(Pyk_zf/(PI_zf+Pn_zf));
        end
        simu_SINR_mrt(idx2,idx) = sinr1/num;
        simu_SINR_zf(idx2,idx) = sinr2/num;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         3. Plotting         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,1,1);
plot(var2_amp_b_range,analytic_SINR_mrt,'r--*','linewidth',1.5); hold on;
% plot(var2_amp_b_range,simu_SINR_mrt,'b--o','linewidth',1.5); hold on;
% legend('QPSK','16QAM','QPSK with permutation','16QAM with permutation')
xlabel('Amplitude Error covairance')
ylabel('SINR(dB)')
% title('Empirical Symbol Error Rate');

subplot(2,1,2);
plot(var2_amp_b_range,analytic_SINR_zf,'r--*','linewidth',1.5); hold on;
% plot(var2_amp_b_range,simu_SINR_zf,'b--o','linewidth',1.5); hold on;
legend({'SINR(Analytic)','SINR(simulated)'})
xlabel('Amplitude Error covairance')
ylabel('SINR(dB)')
% title('Empirical Symbol Error Rate');
toc;




% Defined Functions
function result = phi_func(num)
    result = sqrt(1/(2*pi))*exp(-num^2/2);
end

function result = Phi_func(num)
    result = ( 1+erf(num/sqrt(2)) )/2;
end






