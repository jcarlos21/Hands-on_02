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
    [sOut] = fEstimaCanal(sPar);
    % Parser de variáveis
    vtDistEst = sOut.vtDistEst;
    vtPathLossEst = sOut.vtPathLossEst;
    dNEst = sOut.dNEst;
    vtShadCorrEst = sOut.vtShadCorrEst;
    dStdShadEst = sOut.dStdShadEst;
    dStdMeanShadEst = sOut.dStdMeanShadEst;
    vtDesPequeEst = sOut.vtDesPequeEst;
    vtPrxEst = sOut.vtPrxEst;
    vtXCcdfEst = sOut.vtXCcdfEst;
    vtYCcdfEst = sOut.vtYCcdfEst;
    vtDistLogEst = log10(vtDistEst);
    vtDistLog = log10(vtDist);
    % MSE com Shadowing conhecido
    dMeiaJanela = round((sPar.dW-1)/2);
    vtMSEShad(iw) = immse(vtShadCorr(dMeiaJanela+1 : end-dMeiaJanela ), vtShadCorrEst);
    %
    % MSE com Fading conhecido
    vtMSEFad(iw) = immse(vtFading(dMeiaJanela+1 : end-dMeiaJanela ), vtDesPequeEst);
    %
    disp(['Estimação dos parâmetros de larga escala (W = ' num2str(sPar.dW) '):'])
    disp(['   Expoente de perda de percurso estimado n = ' num2str(dNEst)]);
    disp(['   Desvio padrão do sombreamento estimado = ' num2str(dStdShadEst)]);
    disp(['   Média do sombreamento estimado = ' num2str(dStdMeanShadEst)]);
    disp(['   MSE Shadowing = ' num2str(vtMSEShad(iw))]);
    disp('----');
    disp(' ');
end
% Display informação sobre o estudo das janelas
disp(['Estudo na melhor janela de filtragem']);
disp(['   Janelas utilizadas = ' num2str(vtW)]);
% Melhor janela com Shadowing conhecido
[valBestShad, posBestShad] = min(vtMSEShad);
disp(['   Melhor MSE relativo aos valores reais do Shadowing (melhor janela):'])
disp(['      Melhor janela W = ' num2str(vtW(posBestShad)) ': MSE Shadowing = ' num2str(valBestShad)]);
% Melhor janela com Fading conhecido
[valBestFad, posBestFad] = min(vtMSEFad);
disp(['   Melhor MSE relativo aos valores reais do Fading:'])
disp(['      Melhor janela W = ' num2str(vtW(posBestFad)) ': MSE Shadowing = ' num2str(valBestFad)]);
disp('----------------------------------------------------------------------------------');
disp(' ');