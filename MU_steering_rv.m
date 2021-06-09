clc, clear all;
%% Parameters
fft_len = 64;
mod_type = 4;                     %1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
cp_len = fft_len / 4;
data_len = fft_len*mod_type;
rx_node = 4;
tx_ant = 16;
rx_ant= 1;
N_s = rx_node;
N_d = rx_node;
N_rf = rx_node;
path = 7;
scatter = 10;
snr = -10:3:20;
iter = 100;
%% SCM
model = SCM();
model.n_path = path;
model.n_mray = scatter;
model.ant(rx_ant,tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;
%% Result(preallocation)
result_ber = zeros(1,length(snr));
result_cap = zeros(1,length(snr));
%% Channel initializaton
H = zeros(path, fft_len+cp_len, N_rx * N_d, N_tx);
sel_angle = zeros(4, N_d);
t_He = zeros(path, N_d, N_rf);
He = zeros(fft_len, N_d, N_s);
%% Hybrid beamforming OFDM based
tic
for i = 1:iter
    for j = 1:length(snr)
        for d = 1:N_d
            [temp,rx_angle] = model.FD_channel(fft_len + cp_len);
            H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
            path_gain(:,:) = abs(temp(:,1,1,:));
            
            [~, idx] = max( path_gain );
            idx_ = mode(idx);
            sel_angle(:,d) = rx_angle(:,idx_);
        end
        Wt = steer_precoding(model.fc, model.tx_ant, sel_angle(1:2,:));
        Wr = steer_precoding(model.fc, model.rx_ant, sel_angle(3:4,:),2);
        
        for k = 1:path
            tmp(:,:) = H(k,1,:,:);
            t_He(k,:,:)= Wr.' * tmp * Wt;
        end
        bit = randi([0 1], N_s, data_len);
        Dsym = base_mod(bit, mod_type);
        He = fft(t_He, fft_len, 1);
        [Dsym, ~, Wd] = ZF_precoding(Dsym, He);
        factor = zeros(1,fft_len);
        for k = 1:fft_len
            t_Wd(:,:) = Wd(k,:,:);
            temp = Wt * t_Wd;
            factor(k) = sqrt( trace( temp * temp' ) );
        end
        Dsym = Dsym ./ factor;
        Isym = ifft(Dsym, fft_len, 2) * sqrt(fft_len);
        tx_ofdm = [ Isym(:, fft_len - cp_len + 1 : end) Isym ];
        tx_ofdm = Wt * tx_ofdm;
        [rx_ofdm, No] = awgn_noise( model.FD_fading( tx_ofdm, H ), snr(j) );
        rx_ofdm = Wr.' * rx_ofdm;
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        rx_Dsym = rx_Dsym .* factor;
        rx_bit = base_demod(rx_Dsym, mod_type);
        result_ber(j) = result_ber(j) + sum( sum( bit ~= rx_bit, 2 ), 1 ) / (data_len * N_s);
        rx_ofdm = Wr.' * model.FD_fading( tx_ofdm, H );
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        result_cap(j) = result_cap(j) + mean( sum( log2( 1 + abs(rx_Dsym).^2 / ( No * rx_ant ) ) ) );
        
    end
end
toc
figure(1)
result_ber = result_ber / iter;
result_cap = result_cap / iter;
plot(snr, result_cap, 'r-*');
title('Sum Rate Performance')
legend('ZF-test')
ylabel('Average Spectral Efficiency (bps/Hz)')
xlabel('SNR (dB)')
grid on
hold on

figure(2)
semilogy(snr, result_ber, 'r-*');
title('BER Performance')
legend('ZF-test')
ylabel('BER')
xlabel('SNR (dB)')
grid on
hold on
% axis([snr(1) max(snr) 10^-5 1])