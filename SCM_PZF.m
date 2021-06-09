clc, clear
%% ZF test MU-MIMO
fft_size = 64;
mod_type = 4;                     %1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
cp_size = fft_size / 4;
data_size = fft_size*mod_type;
UE = 4;
tx_ant = 16;
rx_ant = 1;
N_u = UE;
N_s = UE;
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
        %% Precoding(ZF)
        n = 1;
        R = zeros(N_tx,N_tx);         % mean channel (cov)
        for k = 1:fft_size
            H_(:,:) = H(k,:,:);
            R = R + (H_'*H_);
        end
        R = (1/fft_size).*R;
        A = exp(j.*angle(R'));
        
        
        for k = 1:fft_size
            H_(:,:) = H(k,:,:);
            
            He_ = H_*A;
            G = He_'*inv(He_*He_');

            gamma(n) = trace(A*G*G'*A');
            pc_sym(:,k) = G*(sym(:,k)/sqrt(gamma(n)));
            n = n+1;
        end
        ofdm_sym = ifft(pc_sym,fft_size,2)*sqrt(fft_size);
        cp_sym = [ofdm_sym(:,fft_size-cp_size+1:end) ofdm_sym];
        hpc_sym = A*cp_sym;
        
        
        %% Pass Channel
        
        for r = 1 : N_u
            for t = 1 : N_tx
                receive(t,:) = conv(hpc_sym(t,:),h_(:,r,t).');
            end
            hx(r,:) = sum(receive,1);
        end
        %% Receive
        [y, No] = awgn_noise( hx, snr(snr_index) );
        cp_remove = y(:,cp_size+1:fft_size+cp_size);
        y_hat = cp_remove.* sqrt(gamma);
        ofdm_sym_rx = fft(y_hat,fft_size,2)/sqrt(fft_size);
        rx_sym = base_demod(ofdm_sym_rx, mod_type);
        
        s = hx(:,cp_size+1:fft_size+cp_size);
        S = fft(s,fft_size,2)/sqrt(fft_size);
        
        SR(snr_index) = SR(snr_index) + mean( sum( log2( 1 + abs(S).^2 / No )) );
        BER(snr_index) = BER(snr_index) + sum( sum( data ~= rx_sym ) )/ (data_size * N_s);
    end
end
toc
BER = BER / iter;
SR = SR / iter;
%% Plot
figure(1)
plot(snr, SR, 'g-o');
title('Sum Rate Performance')
legend('ZF')
ylabel('Average Spectral Efficiency (bps/Hz)')
xlabel('SNR (dB)')
grid on
hold on

figure(2)

semilogy(snr, BER, 'g-o');
title('BER Performance')
legend('ZF')
ylabel('BER')
xlabel('SNR (dB)')
grid on
hold on