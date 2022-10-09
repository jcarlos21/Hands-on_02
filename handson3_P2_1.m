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
% Chama função que gera o canal sintético
[vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm] = fGeraCanal(sPar);
% Transforma potência em mWatts
vtPtrxmW = 10.^(vtPrxdBm/10);
nSamples = length(vtPtrxmW);
% Vetores para canal estimado
vtDesLarga = [];
vtDesPequeEst = [];
%
% Cálculo do desvanecimenro lento e rápido
dMeiaJanela = round((sPar.dW-1)/2);  % Meia janela
ij = 1;
for ik = dMeiaJanela + 1 : nSamples - dMeiaJanela
    % Desvanecimento de larga escala: perda de percurso + sombreamento [dB]
    vtDesLarga(ij) = 10*log10(mean(vtPtrxmW(ik-dMeiaJanela:ik+dMeiaJanela)));
    % Desvanecimento de pequena escala [dB]
    vtDesPequeEst(ij) = vtPrxdBm(ik)-vtDesLarga(ij);
    ij = ij + 1;
end
%
% Cálculo da envoltória normalizada (para efeitos de cálculo do fading)
indexes = dMeiaJanela+1 : nSamples-dMeiaJanela;
%vtPrxW = ((10.^(vtPrxdBm(indexes)./10))/1000);
vtPtrxmWNew = 10.^(vtPrxdBm(indexes)/10);
desLarga_Lin = (10.^(vtDesLarga(1:length(indexes))./10));
envNormal = sqrt(vtPtrxmWNew)./sqrt(desLarga_Lin);
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
disp(['Estimação dos parâmetros de larga escala (W = ' num2str(sPar.dW) '):'])
disp(['   Expoente de perda de percurso estimado n = ' num2str(dNEst)]);
% Perda de percurso estimada para os pontos de medição
vtPathLossEst = polyval(dCoefReta,vtDistLogEst);  
%
% Sombreamento
vtShadCorrEst = vtDesLarga - vtPathLossEst;
% Calcula a variância do sombreamento estimado
stdShad = std(vtShadCorrEst);
meanShad = mean(vtShadCorrEst);
disp(['   Desvio padrão do sombreamento estimado = ' num2str(stdShad)]);
disp(['   Média do sombreamento estimado = ' num2str(meanShad)]);
vtPathLossEst = - vtPathLossEst;
vtPrxEst = sPar.txPower - vtPathLossEst + vtShadCorrEst + vtDesPequeEst;
%
% Estimação da CDF do desvanecimento de pequena escala
% Cálculo dos pontos do eixo x da cdf (espacamento igual entre os pontos)
vtn = 1 : sPar.nCDF;
xCDF = 1.2.^(vtn-1) * 0.01;
%
% Cálculo da CDF
den = 0;
cdffn=zeros(1,sPar.nCDF);
for ik = 1:sPar.nCDF
    for ij = 1:length(envNormal)
        if envNormal(ij) <= xCDF(ik)
            den = den + 1;
        end
        cdffn(ik) = cdffn(ik) + den;
        den = 0;
    end
end
%
% Monta estrutura do histograma
xccdfEst = 20.*log10(xCDF);
yccdfEst = cdffn/(cdffn(end)); 
% Figuras do canal estimado
figure;
% Potência recebida com canal completo
plot(vtDistLogEst,vtPrxEst); hold all;
% Potência recebida com path loss
plot(vtDistLogEst,sPar.txPower-vtPathLossEst,'linewidth', 2)
% Potência recebida com path loss e shadowing
plot(vtDistLogEst,sPar.txPower-vtPathLossEst+vtShadCorrEst,'linewidth', 2)
%title('Canal estimado: Potência recebida no receptor vs. log da distância')
xlabel('log_{10}(d)');
ylabel('Potência [dBm]');
legend('Prx canal completo', 'Prx (somente perda de percurso)', 'Prx (perda de percurso + sombreamento)');
title('Prx original vs estimada');
%
% Figura do Path loss (original vs estimado) 
figure;
plot(vtDistLog,-vtPathLoss);hold on;plot(vtDistLogEst,-vtPathLossEst);
legend('Path Loss original','Path Loss estimado');
title('Perda de percurso original vs estimada');
%
% Figura do Sombreamento (original vs estimado) 
figure;
plot(vtDistLog,vtShadCorr);hold on;plot(vtDistLogEst,vtShadCorrEst);
legend('Shadowing original','Shadowing estimado');
title('Sombreamento original vs estimada');
%
% Figura do Fading (original vs estimado) 
figure;
plot(vtDistLog,vtFading);hold on;plot(vtDistLogEst,vtDesPequeEst);
legend('Fading original','Fading estimado');
title('Fading original vs estimada');
%
% Plot das CDFs normalizadas Nakagami (assumindo que sabemos que o canal é m-Nakagami)- para vários valores de m
figure;
plot( xccdfEst, yccdfEst, '--' );
legendaNaka = [{'CDF das amostras'}];
hold all;
vtm = [1 2 4 6];
xCDF = 10.^(xccdfEst/20);
tam_dist = length(gammainc(1*xCDF.^2,1)); % Tamanho da distribuição
for ik = 1:length(vtm)%
    im = vtm(ik);
    cdfnaka(ik,1:tam_dist) = gammainc(im*xCDF.^2,im);
    semilogy(20.*log10(xCDF),cdfnaka(ik,:));
    legendaNaka = [legendaNaka ; {['m = ' num2str(vtm(ik))]}];
end
legend(legendaNaka);
axis([-30 10 1e-5 1]);
title('Estudo do fading com o conhecimento da distribuição');
xlabel('x');
ylabel('F(x)');