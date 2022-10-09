close all;clear all;clc;
load("Prx_Real_2021_1.mat");
txPower = 0;
vtDist = dPath';
vtPrxdBm = dPrx';
%
% Transforma potência em mWatts
vtPtrxmW = 10.^(vtPrxdBm/10);
nSamples = length(vtPtrxmW);
%
vtW = [2 5 10];
for iw = 1: length(vtW)
    % Configura valor da janela de filtragem
    dW = vtW(iw);
    chFileName = ['DesvanecimentoPequenaEscala', num2str(dW)];
    % Vetores para canal estimado
    vtDesLarga = []; % Desvanecimento de larga escala
    vtDesPequeEst = []; % Desvanecimento de pequena escala
    %
    % Cálculo do desvanecimenro lento e rápido
    dMeiaJanela = round((dW-1)/2);  % Meia janela
    ij = 1;
    %
    for ik = dMeiaJanela + 1 : nSamples - dMeiaJanela
        % Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
        vtDesLarga(ij) = 10*log10(mean(vtPtrxmW(ik-dMeiaJanela:ik+dMeiaJanela)));
        % Desvanecimento de pequena escala [dB]
        vtDesPequeEst(ij) = vtPrxdBm(ik)-vtDesLarga(ij);
        ij = ij + 1;
    end
    %
    % Armazenamento das informações do vetor de desvanecimento em pequena escala
    dlmwrite([chFileName '.txt'], [vtDesPequeEst'], 'delimiter', '\t')
    % writematrix([chFileName '.txt'], [vtDesPequeEst'], 'delimiter', '\t')
    %
    % Cálculo da envoltória normalizada (para efeitos de cálculo do fading)
    indexes = dMeiaJanela+1 : nSamples-dMeiaJanela;
    vtPtrxmWNew = 10.^(vtPrxdBm(indexes)/10);
    desLarga_Lin = (10.^(vtDesLarga(1:length(indexes))./10));
    vtEnvNorm = sqrt(vtPtrxmWNew)./sqrt(desLarga_Lin);
    %
    % Ajuste no tamanho dos vetores devido a filtragem
    vtDistEst = vtDist( dMeiaJanela+1 : nSamples-dMeiaJanela );
    vtPrxdBm = vtPrxdBm( dMeiaJanela+1 : nSamples-dMeiaJanela );
    %
    % Cálculo reta de perda de percurso
    vtDistLog = log10(vtDist);
    vtDistLogEst = log10(vtDistEst);
    % Cálculo do coeficientes da reta que melhor se caracteriza a perda de percurso
    dCoefReta = polyfit(vtDistLogEst,vtPrxdBm,1); 
    % Expoente de perda de percurso estimado
    dNEst = -dCoefReta(1)/10;
    %
    disp(['Valor estimado dos parâmetros em larga escala (W = ' num2str(dW) '):']);
    disp(['Expoente estimado de perda de percurso n = ' num2str(dNEst)]);
    %
    % Perda de percurso estimada para os pontos de medição
    vtPathLossEst = polyval(dCoefReta,vtDistLogEst);  
    %
    % Sombreamento
    vtShadCorrEst = vtDesLarga - vtPathLossEst;
    % Calcula a variância do sombreamento estimado
    dStdShadEst = std(vtShadCorrEst);
    dStdMeanShadEst = mean(vtShadCorrEst);
    vtPathLossEst = - vtPathLossEst;
    vtPrxEst = txPower - vtPathLossEst + vtShadCorrEst + vtDesPequeEst;
    %
    disp(['Desvio padrão estimado do sombreamento = ' num2str(dStdShadEst)]);
    disp(['Média estimada do sombreamento = ' num2str(dStdMeanShadEst)]);
    %
    % Plot do desvanecimento
    figure;
    %
    % Potência recebida com canal completo
    plot(vtDistLogEst,vtPrxEst); hold all;
    % Potência recebida com path loss
    plot(vtDistLogEst,txPower-vtPathLossEst,'linewidth', 2)
    % Potência recebida com path loss e shadowing
    plot(vtDistLogEst,txPower-vtPathLossEst+vtShadCorrEst,'linewidth', 2)
    %
    title(['Potência recebida em RX x log da distância (W = ' num2str(dW) '):'])
    xlabel('log_{10}(d)');
    ylabel('Potência [dBm]');
    legend('PRX Canal Completo', 'PRX (apens com perda de percurso)', 'PRX (com perda de percurso e sombreamento)');
end