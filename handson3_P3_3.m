close all;clear all;clc;
% Estimação com o conhecimento do sombreamento e do fading
% Parâmetros para geração do canal sintético
sPar.d0 = 5;                     % distância de referência d0
sPar.P0 = 0;                     % Potência medida na distância de referência d0 (em dBm)
sPar.nPoints = 50000;            % Número de amostras da rota de medição
sPar.totalLength = 100;          % Distância final da rota de medição
sPar.n = 4;                      % Expoente de perda de percurso
sPar.sigma = 6;                  % Desvio padrão do shadowing em dB
sPar.shadowingWindow = 200;      % Tamanho da janela de correlação do shadowing (colocar em função da distância de correlação)
sPar.m = 4;                      % Parâmetro de Nakagami
sPar.txPower = 0;                % Potência de transmissão em dBm
sPar.nCDF = 40;                  % Número de pontos da CDF normalizada
sPar.dW = 100;                   % Janela de estimação do sombreamento
sPar.chFileName  = 'Prx_sintetico';
% Distância entre pontos de medição
sPar.dMed = sPar.totalLength/sPar.nPoints;
%
% Chama função que gera o canal sintético
[vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm] = fGeraCanal(sPar);
%
% Mostra informações do canal sintético
disp('Canal sintético:')
disp(['   Média do sombreamento: ' num2str(mean(vtShadCorr)) ]);
disp(['   Std do sombreamento: ' num2str(std(vtShadCorr)) ]);
disp(['   Janela de correlação do sombreamento: ' num2str(sPar.shadowingWindow) ' amostras' ]);
disp(['   Expoente de path loss: ' num2str(sPar.n) ]);
disp(['   m de Nakagami: ' num2str(sPar.m) ]);
%
% Várias janelas de filtragem para testar a estimação
vtW = [10 50 150 200];
for iw = 1: length(vtW)
    % Configura valor da janela de filtragem
    sPar.dW = vtW(iw);
    % Chama função que estima o canal sintético
    sOut(iw) = fEstimaCanal(sPar);
    % Parser de variáveis
    vtDistEst = sOut(iw).vtDistEst;
    vtPathLossEst = sOut(iw).vtPathLossEst;
    dNEst = sOut(iw).dNEst;
    vtShadCorrEst = sOut(iw).vtShadCorrEst;
    dStdShadEst = sOut(iw).dStdShadEst;
    dStdMeanShadEst = sOut(iw).dStdMeanShadEst;
    vtDesPequeEst = sOut(iw).vtDesPequeEst;
    vtPrxEst = sOut(iw).vtPrxEst;
    vtXCcdfEst = sOut(iw).vtXCcdfEst;
    vtYCcdfEst = sOut(iw).vtYCcdfEst;
    vtDistLogEst = log10(vtDistEst);
    vtDistLog = log10(vtDist);
end
% Estimação cega via MLE
disp(' ')
disp('Estimação do Fading para várias janelas (estudo númerico sem conhecimento a priori do canal)');
disp('Resultados com fitdist do Matlab')
for iw = 1:length(vtW)%
    disp(['Janela W = ' num2str(vtW(iw))]);
    %
    sNaka(iw) = fitdist([sOut(iw).vtEnvNorm]','Nakagami');
    disp(['  Nakagami: m = ' num2str(sNaka(iw).mu) ', omega = ' num2str(sNaka(iw).omega)]);
    %
    sRice(iw) = fitdist([sOut(iw).vtEnvNorm]','Rician');
    K_rice = (sRice(iw).s)^2/(2*sRice(iw).sigma^2);
    disp(['  Rice: K = ' num2str(K_rice)]);
    %
    sRay(iw) = fitdist([sOut(iw).vtEnvNorm]','Rayleigh');
    disp(['  Rayleigh: sigma = ' num2str(sRay(iw).B)]);
    %
    sWei(iw) = fitdist([sOut(iw).vtEnvNorm]','Weibull');
    disp(['  Weibull: k = ' num2str(sWei(iw).B) ', lambda = ' num2str(sWei(iw).A)]);
    disp(' ')
end