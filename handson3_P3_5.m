close all;clear all;clc;
addpath('./fitmethis')
%addpath('./CODES/HD_03_MATLAB/fitmethis')
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
% Teste KS
sDistNames = [{'Weibull'};{'Rician'};{'Rayleigh'};{'Nakagami'}];
% Test KS (h = 0, não rejeita H0 e os dados podem ser da distribuição especificada;
% k é o maior diferença engre a CDF dos dados e da distribuição especificada. Quanto menor k, melhor o fit)
for iw = 1:length(vtW)
    disp(' ')
    disp(['Janela W = ' num2str(vtW(iw)) ]);
    for ik = 1:length(sDistNames)
        data = [sOut(iw).vtEnvNorm]';
        pd = fitdist(data,sDistNames{ik});
        x = linspace(min(data),max(data),length(data));
        tCDF = [x' cdf(pd,x)'];
        [h,p,k,c] = kstest(data,'CDF',tCDF);
        sKtest(ik,iw).h = h;
        sKtest(ik,iw).p = p;
        sKtest(ik,iw).k = k;
        sKtest(ik,iw).c = c;
        disp(['   Distribuição ' sDistNames{ik} ': k = ' num2str(k) ', p-value = ' num2str(p)]);
        % Resultado do teste KS
        if (h == 0)
            disp('     h = 0 => Não rejeita a hipótese H0 com nível de significância $\alpha$ = 5% (p > 0.05).');
        elseif (h == 1)
            disp('     h = 1 => Rejeita a hipótese H0 com nível de significância $\alpha$ = 5% (p < 0.05).');
        end
    end
end
