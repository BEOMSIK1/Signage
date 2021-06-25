
%% SVD Analog precoder
%% Parameters
fft_size = 64;
mod_type = 4;                     %1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
cp_size = fft_size / 4;
data_size = fft_size*mod_type;
tx_ant = 8;
rx_ant = 1;
N_u = 8;
N_rf = N_u:2*N_u;
N_s = N_u;
snr = 20;
path = 7;
scatter = 10;
iter = 500;

%% RF parameters
P_t = 0.500;  % Watt
P_rf = 0.300;
P_ps = 0.010; %        1bit : 2.5mW  2bit: 10mW  inf: 78mW
B = 2;        %         resolution (bit)  1  2  5(inf)
P_bb = 0.200;
P_sw = 0.0025;
P_sp = 0.020;
P_cb = 0.020;
%% SCM
model = SCM();
model.n_path = path;
model.n_mray = scatter;
model.ant(rx_ant,tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;
%% test
model.asd = 15;
model.zsd = 5;
model.asa = 15;
model.zsa = 5;

%model.fc = 28*10^9;
%model.fs = 0.25*10^9;

model.tx_ant(3) = 0.5;
model.rx_ant(3) = 0.5;
model.los = 0;


%% Initialize

BER = zeros(1, length(N_rf) );
SR1 = zeros(1, length(N_rf) );

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
    for rf_index = 1:length(N_rf)
        %% Precoding
        n = 1;
        %% RF
        R = zeros(N_tx,N_tx);         % mean channel (cov)
        for k = 1:fft_size
            H_(:,:) = H(k,:,:);
            R = R + (H_'*H_);
        end
        R = R/fft_size;
        [~,~,V] = svd(R);
        F = V(:,1:N_rf(rf_index));
        F_hat = zeros(N_tx,N_rf(rf_index));
        
        for nt = 1:N_tx
            for nrf = 1:N_rf(rf_index)
                if B == 2
                    if angle(F(nt,nrf))>=0 && angle(F(nt,nrf))<(pi/2)
                        
                        F_hat(nt,nrf) = exp(1j*pi/4);
                        
                    elseif angle(F(nt,nrf))>=(pi/2) && angle(F(nt,nrf))<=(pi)
                        
                        F_hat(nt,nrf) = exp(1j*3*pi/4);
                        
                    elseif angle(F(nt,nrf))>=(-pi) && angle(F(nt,nrf))<(-pi/2)
                        
                        F_hat(nt,nrf) = exp(1j*-3*pi/4);
                        
                    elseif angle(F(nt,nrf))>=(-pi/2) && angle(F(nt,nrf))<(0)
                        
                        F_hat(nt,nrf) = exp(1j*-pi/4);
                        
                    end
                    
                elseif B == 1
                    
                    if angle(F(nt,nrf))>0 && angle(F(nt,nrf))<=(pi)
                        F_hat(nt,nrf) = exp(1j*pi/2);
                        
                        
                    elseif angle(F(nt,nrf))>=(-pi) && angle(F(nt,nrf))<=(0)
                        F_hat(nt,nrf) = exp(1j*-pi/2);
                        
                        
                    end
                    
                elseif B == 5
                    F_hat(nt,nrf) = exp(1j*angle(F(nt,nrf)));  % inf
                    
                end
            end
        end
        
        %% BB
        pc_sym = zeros(N_rf(rf_index),fft_size);
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
        [y, No] = awgn_noise( hx, snr);
        cp_remove = y(:,cp_size+1:fft_size+cp_size);
        y_hat = cp_remove.* sqrt(NF);
        ofdm_sym_rx = fft(y_hat,fft_size,2)./sqrt(fft_size);
        rx_data = base_demod(ofdm_sym_rx, mod_type);
        
        s = hx(:,cp_size+1:fft_size+cp_size);
        S = fft(s,fft_size,2)/sqrt(fft_size);
        
        SR1(rf_index) = SR1(rf_index) + mean( sum( log2( 1 + abs(S).^2 / (No*rx_ant) )) )/iter;
        %         num_error(i,rf_index) = biterr(data,rx_data);
        %
        ee_(i,rf_index) = SR1(rf_index)/(P_t+N_rf(rf_index)*P_rf+N_rf(rf_index)*N_tx*P_ps+N_rf(rf_index)*P_sp+N_tx*P_cb+P_bb);
        
        
        
        %         BER(snr_index) = BER(snr_index) + sum( sum( data ~= rx_data ) )/ (data_size * N_s);
    end
end
toc
ee = sum(ee_)/iter;
%% Plot
figure(1)
plot(N_rf,ee, '-o')
title('Energy Efficiency Performance')
legend('Opt')
ylabel('Energy Efficiency (bps/Hz/W)')
xlabel('RF chain')
grid on
hold on

figure(2)
plot(N_rf, SR1, '-o');
title('Sum Rate Performance')
legend('Opt')
ylabel('Average Spectral Efficiency (bps/Hz)')
xlabel('RF chain')
grid on
hold on
