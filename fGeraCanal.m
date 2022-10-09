envolt%%file fGeraCanal.m
function [vtDist, vtPathLoss, vtShadCorr, vtFading, vtPrxdBm] = fGeraCanal(sPar)
% Propósito: Gerar canal composta de path loss, shadowing e desvanecimento plano
%
% ENTRADAS: na estrutura sParams
%    - nPoints: Número de amostras
%    - totalLength: distância máxima da rota
%    - P0: potência de referência medida na distância d0
%    - d0: distância d0 de referência
%    - n: expoente de perda de percurso
%    - sigma: desvio padrão do sombreamento lognormal [dB]
%    - Tamanho da janela de correlação do sombreamento [amostras]
%    - m: parâmetro de Nakagami
%    - dMed: distância entre pontos de medição (totalLength/nPoints)
%    - txPower: potência de transmissão em dBm
%
% SAÍDAS:
%    - vtDist: pontos de medição [m]
%    - vtPathLoss: amostras da perda de percurso
%    - vtShadCorr: amostras do somrbeamento
%    - fading: amostras do desvanecimento de pequena escala
%    - vtPrx: potência recebida com o canal completo
%
% Parser dos parâmetros de entrada
nPoints = sPar.nPoints;
totalLength = sPar.totalLength;
P0 = sPar.P0;
d0 = sPar.d0;
n = sPar.n;
sigma = sPar.sigma;
shadowingWindow = sPar.shadowingWindow;
m = sPar.m;
dMed = sPar.dMed;
txPower = sPar.txPower;
%
% Distância do transmissor (além da distância de referência)
d = d0:dMed:totalLength;
nSamples = length(d);
%
% Geração da Perda de percurso (determinística)
vtPathLoss = P0 + 10*n*log10(d./d0);
%
% Geração do Sombreamento
nShadowSamples = floor(nSamples/shadowingWindow);
shadowing = sigma*randn(1,nShadowSamples);
% Amostras para a última janela
restShadowing = sigma*randn(1,1)*ones(1,mod(nSamples,shadowingWindow));
% Repetição do mesmo valor de sombreamento durante a janela de correlação
shadowing = ones(shadowingWindow,1)*shadowing;
% Amostras organizadas em um vetor
shadowing = [reshape(shadowing,1,nShadowSamples*shadowingWindow),restShadowing];
% Filtragem para evitar variação abrupta do sombreamento
jan = shadowingWindow/2;
iCont = 1;
for i = jan+1:nSamples-jan,
    vtShadCorr(iCont) = mean(shadowing(i-jan:i+jan)); %diminuir a variação brusca do sombreamento
    iCont = iCont+1;
end
% Ajuste do desvio padrão depois do filtro de correlação do sombreamento
vtShadCorr = vtShadCorr*std(shadowing)/std(vtShadCorr);
vtShadCorr = vtShadCorr - mean(vtShadCorr)+ mean(shadowing);
%
% Geração do desvanecimento de pequena escala
% Nakagami fading
nakagamiPdf = @(x)((2.*m.^m)./(gamma(m))).*x.^(2.*m-1).*exp(-(m.*x.^2));
nakagamiNormEnvelope = slicesample(1,nSamples,'pdf',nakagamiPdf);
nakagamiSamp=20.*log10(nakagamiNormEnvelope');
%
% Cálculo da Potência recebida
txPower = txPower*ones(1,nSamples);
% Ajuste do número de amostras devido a filtragem
txPower = txPower(jan+1:nSamples-jan);
vtPathLoss = vtPathLoss(jan+1:nSamples-jan);
vtFading = nakagamiSamp(jan+1:nSamples-jan);
vtDist = d(jan+1:nSamples-jan);
% Potência recebida
vtPrxdBm = txPower-vtPathLoss+vtShadCorr+vtFading;
% Salva variáveis do canal no Matlab
save([sPar.chFileName '.mat'],'vtDist', 'vtPathLoss', 'vtShadCorr', 'vtFading', 'vtPrxdBm');