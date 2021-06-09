function [precoded_X,gamma] = MU_ZF(P,X_mod,H)

%% Precoding (ZF)    
tx_power=1*(10^(P/10));                                    % rho 
gamma=trace(inv(H*H'))./tx_power;                          % normalize
precoded_matrix = (H'*inv(H*H'))./sqrt(gamma);             % 전처리 행렬
precoded_X = precoded_matrix*X_mod;                         % 사용자 신호 전처리
