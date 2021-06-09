clc, clear, close all;
%% Parameters
Nt = 16;                                     % 송신안테나 개수
Nr = 1;                                      % 수신안테나 개수 / 유저
K = 8;                                       % 유저의 수
modulation_order = 6;                        % 1:BPSK  2:QPSK  4: 16QAM  6: 64QAM  8: 256QAM
symbol_size = 128;
data_size = symbol_size*modulation_order;
SNR=0:3:30;                                  % SNR (dB)
P=0;                                         % transmitpower=1
iteration=1000;
%% Channel inversion
for SNR_index=1:length(SNR)
    tic
    for iter=1:iteration
        %% Channel genrate
        H=(randn(Nr*K,Nt)+j*randn(Nr*K,Nt))/sqrt(2);
        cnd(iter,SNR_index)=cond(H);
        
        %% TEST
        A = exp(j.*angle(H'));
        He = (H*A);
        %% Data
        X=randi([0 1], [Nr*K, data_size]);                                    % N_t개의 서로 다른 데이터 (송신측)
        X_mod = base_mod(X, modulation_order)./sqrt(K);                              % modulation
         %% ZF
        [precoded_X,gamma]=MU_ZF(P,X_mod,He);
        %% TEST
        

        
        %% RX
        Y_received = awgn_noise(He*precoded_X,SNR(SNR_index));       % H*전처리된 각 사용자 신호/sqrt(gamma) + noise
        Y = Y_received.*sqrt(gamma).*sqrt(K);                                % 정규화 factor를 곱해주어 사용자 신호만 추출
        Y_demod = base_demod(Y,modulation_order);
        %% Error
        num_error_ZF(iter,SNR_index)=biterr(X,Y_demod);                      % x,y의 총 error 발생 횟수
    end
    toc
end
cnd5=find(cnd(:,8)>5);
cnd4=find(cnd(:,8)>=4&cnd(:,8)<5);
cnd3=find(cnd(:,8)>=3&cnd(:,8)<4);
cnd2=find(cnd(:,8)>=2&cnd(:,8)<3);
mean5=mean(num_error_ZF(cnd5,8));
mean4=mean(num_error_ZF(cnd4,8));
mean3=mean(num_error_ZF(cnd3,8));
mean2=mean(num_error_ZF(cnd2,8));
error_rate_ZF=(sum(num_error_ZF,1)/(data_size*Nr*K))/iteration;
%% graph
semilogy(SNR,error_rate_ZF,'-o')
title('BER Performance'), xlabel('SNR(dB)'),ylabel('BER'),legend('ZF')
grid on