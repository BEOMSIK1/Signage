clc, clear
%% SVD Analog precoder
%% Parameters
fft_size = 64;
mod_type = 4;                     %1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
cp_size = fft_size / 4;
data_size = fft_size*mod_type;
tx_ant = 64;
rx_ant = 1;
N_u = 16;
N_rf = 2*N_u;
N_s = N_u;
snr = -10:5:20;
path = 7;
scatter = 10;
iter = 300;

%% SCM
model = SCM();
model.n_path = path;
model.n_mray = scatter;
model.ant(rx_ant,tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;
%% test
model.asd = 3;
model.zsd = 3;
model.asa = 3;
model.zsa = 3;

% model.fc = 30*10^9;
% model.fs = 0.25*10^9;

model.tx_ant(3) = 0.5;
model.rx_ant(3) = 0.5;
model.los = 0;


%% Initialize

 BER = zeros(1, length(snr) );
 SR = zeros(1, length(snr) );

%% Iter
tic
for i = 1:iter
    %% channel
    for d= 1:N_u
        temp = model.FD_channel(fft_size + cp_size);
        h(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
    end
    h_(:,:,:) = h(:,1,:,:);
    H = fft(h_, fft_size, 1);
    %% Data
    data = randi([0 1], N_s, data_size);
    sym = base_mod(data, mod_type);
    %% iter (snr)
    for snr_index = 1:length(snr)
       %% Precoding
        n = 1;
        F_hat = zeros(N_tx,N_rf);
        %% RF
        R = zeros(N_tx,N_tx);         % mean channel (cov)
        for k = 1:fft_size
            H_(:,:) = H(k,:,:);
            R = R + (H_'*H_);
        end
        R = R/fft_size;
        [~,~,V] = svd(R);
        F = V(:,1:N_rf);
        
        psn1 = find(angle(F)>=0 & angle(F)<(pi/2));
        psn2 = find(angle(F)>=(pi/2) & angle(F)<=(pi));
        psn3 = find(angle(F)>=(-pi) & angle(F)<(-pi/2));
        psn4 = find(angle(F)>=(-pi/2) & angle(F)<(0));
        
        f_mean1 = mean(angle(F(psn1)));
        f_mean2 = mean(angle(F(psn2)));
        f_mean3 = mean(angle(F(psn3)));
        f_mean4 = mean(angle(F(psn4)));
       
%                 f_mean1 = pi/4;
%         f_mean2 = 3*pi/4;
%         f_mean3 = -3*pi/4;
%         f_mean4 = -pi/4;
        
        B = 2;
        for ant = 1:N_tx
             sw = zeros(2^B,N_rf);
             sp = zeros(1,2^B);
            for rf = 1:N_rf
                if angle(F(ant,rf))>=0 & angle(F(ant,rf))<(pi/2)
                    sw(1,rf) = 1; 
                    sp(1) = 1;
                elseif angle(F(ant,rf))>=(pi/2) & angle(F(ant,rf))<=(pi)
                    sw(2,rf) = 1; 
                    sp(2) = 1;
                elseif angle(F(ant,rf))>=(-pi) & angle(F(ant,rf))<(-pi/2)
                    sw(3,rf) = 1; 
                    sp(3) = 1;
                elseif angle(F(ant,rf))>=(-pi/2) & angle(F(ant,rf))<(0)
                    sw(4,rf) = 1; 
                    sp(4) = 1;
                end
            end
            
            pss = diag([exp(1j*f_mean1),exp(1j*f_mean2),exp(1j*f_mean3),exp(1j*f_mean4)]);

            F_hat(ant,:) = sp * pss * sw;
            
        end
           
        %% BB
        for k = 1:fft_size   
            H_(:,:) = H(k,:,:);
            He_ = H_*F_hat;
            G = (He_'*inv(He_*He_'));
            NF(n) = trace(F_hat*G*G'*F_hat');
            pc_sym(:,k) = G*(sym(:,k)/sqrt(NF(n)));
            n = n+1;
        end
        ofdm_sym = ifft(pc_sym,fft_size,2)*sqrt(fft_size);
        cp_sym_ = [ofdm_sym(:,fft_size-cp_size+1:end) ofdm_sym];
        cp_sym = F_hat*cp_sym_;
        
        
        %% Pass Channel
        for r = 1 : N_u
            for t = 1 : N_tx
                receive(t,:) = conv(cp_sym(t,:),h_(:,r,t).');
            end
            hx(r,:) = sum(receive,1);
        end
        %% Receive
        [y, No] = awgn_noise( hx, snr(snr_index) );
        cp_remove = y(:,cp_size+1:fft_size+cp_size);
        y_hat = cp_remove.* sqrt(NF);
        ofdm_sym_rx = fft(y_hat,fft_size,2)./sqrt(fft_size);
        rx_data = base_demod(ofdm_sym_rx, mod_type);
        
        s = hx(:,cp_size+1:fft_size+cp_size);
        S = fft(s,fft_size,2)/sqrt(fft_size);
        
        SR(snr_index) = SR(snr_index) + mean( sum( log2( 1 + abs(S).^2 / (No*rx_ant) )) );
        num_error(i,snr_index) = biterr(data,rx_data);

%         BER(snr_index) = BER(snr_index) + sum( sum( data ~= rx_data ) )/ (data_size * N_s);
    end
end
toc
% BER = BER / iter;
SR = SR / iter;
BER = (sum(num_error,1)/(data_size*N_s))/iter;
%% Plot
figure(1)
plot(snr, SR, 'r-o');
title('Sum Rate Performance')
legend('test')
ylabel('Average Spectral Efficiency (bps/Hz)')
xlabel('SNR (dB)')
grid on
hold on

figure(2)

semilogy(snr, BER, 'r-o');
title('BER Performance')
legend('test')
ylabel('BER')
xlabel('SNR (dB)')
grid on
hold on
