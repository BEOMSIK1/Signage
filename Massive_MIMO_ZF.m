function [num_error,Sum_rate] = Massive_MIMO_ZF(P,power_index,Nt,K,SNR,modulation_order,data_size)

%% Channel genrate
H = (randn(K,Nt)+j*randn(K,Nt))/sqrt(2);
%% Data
X=randi([0 1], [K, data_size]);                                    % N_t개의 서로 다른 데이터 (송신측)
X_mod = base_mod(X, modulation_order);                              % modulation
%% Precoding (ZF)
tx_power=1*(10^(P(power_index)/10));                       % rho
gamma=trace(inv(H*H'))./tx_power;                          % normalize
precoded_matrix = (H'*inv(H*H'))./sqrt(gamma);             % 전처리 행렬
%% Tx-Rx
precoded_X = precoded_matrix*X_mod;              % 사용자 신호 전처리
Y_received = awgn_noise(H*precoded_X,SNR);       % H*전처리된 각 사용자 신호/sqrt(gamma) + noise
Y = Y_received.*sqrt(gamma);                     % 정규화 factor를 곱해주어 사용자 신호만 추출
Y_demod_ZF = base_demod(Y,modulation_order);
%% Error
num_error=biterr(X,Y_demod_ZF);
%% Sum rate
Sum_rate=K*log2(1+(tx_power*(Nt-K)/K));