% Parâmetros para geração do canal sintético
sPar.d0 = 5;                     % distância de referência d0
sPar.P0 = 0;                     % Potência medida na distância de referência d0
sPar.nPoints = 50000;            % Número de amostras da rota de medição
sPar.totalLength = 100;          % Distância final da rota de medição
sPar.n = 4;                      % Expoente de perda de percurso
sPar.sigma = 6;                  % Desvio padrão do shadowing em dB
sPar.shadowingWindow = 200;      % Tamanho da janela de correlação do shadowing (colocar em função da distância de correlação)
sPar.m = 4;                      % Parâmetro de Nakagami
sPar.txPower = 0;                % Potência de transmissão em dBm
sPar.nCDF = 40;                  % Número de pontos da CDF normalizada
sPar.chFileName  = 'Prx_sintetico';
% Distância entre pontos de medição
sPar.dMed = sPar.totalLength/sPar.nPoints;
% Vetor de distâncias do transmissor (além da distância de referência)
vtDist = sPar.d0:sPar.dMed:sPar.totalLength;
% Número de amostras geradas
nSamples = length(vtDist);
% Geração da Perda de percurso (determinística)
vtPathLoss = sPar.P0 + 10*sPar.n*log10(vtDist./sPar.d0);
% Geração do Sombreamento (V.A. Gaussiana com média zero e desvio padrão sigma)
nShadowSamples = floor(nSamples/sPar.shadowingWindow);
vtShadowing = sPar.sigma*randn(1,nShadowSamples);
% Amostras para a última janela
restShadowing = sPar.sigma*randn(1,1)*ones(1,mod(nSamples,sPar.shadowingWindow));
% Repetição do mesmo valor de sombreamento durante a janela de correlação
vtShadowing = ones(sPar.shadowingWindow,1)*vtShadowing;
% Amostras organizadas em um vetor
vtShadowing = [reshape(vtShadowing,1,nShadowSamples*sPar.shadowingWindow),restShadowing];
% Filtragem para evitar variação abrupta do sombreamento (filtro média móvel)
% O sombreamento tem menos "2*jan" amostras devido a filtragem
jan = sPar.shadowingWindow/2;
iCont = 1;
for i = jan+1:nSamples-jan,
    vtShadCorr(iCont) = mean(vtShadowing(i-jan:i+jan));
    iCont = iCont+1;
end
% Ajuste do desvio padrão depois do filtro de correlação do sombreamento
vtShadCorr = vtShadCorr*std(vtShadowing)/std(vtShadCorr);
vtShadCorr = vtShadCorr - mean(vtShadCorr)+ mean(vtShadowing);
%
% Geração do desvanecimento de pequena escala: Nakagami fading
% PDF da envolvtória normalizada
fpNakaPdf = @(x)((2.*sPar.m.^sPar.m)./(gamma(sPar.m))).*x.^(2.*sPar.m-1).*exp(-(sPar.m.*x.^2));
% Gerador de números aleatórios com distribuição Nakagami
vtNakagamiNormEnvelope = slicesample(1,nSamples,'pdf',fpNakaPdf);
% Fading em dB (Potência)
vtNakagamiSampdB = 20.*log10(vtNakagamiNormEnvelope');
%
% Cálculo da Potência recebida
vtTxPower = sPar.txPower*ones(1,nSamples);
% Ajuste do número de amostras devido ao filtro de correlação do
% sombreamento (tira 2*"Jan" amostras)
vtTxPower = vtTxPower(jan+1:nSamples-jan);
vtPathLoss = vtPathLoss(jan+1:nSamples-jan);
vtFading = vtNakagamiSampdB(jan+1:nSamples-jan);
vtDist = vtDist(jan+1:nSamples-jan);
% Potência recebida
vtPrx = vtTxPower-vtPathLoss+vtShadCorr+vtFading;
%
% Salvamento dos dados
%    Para excel:
dlmwrite([sPar.chFileName '.txt'], [vtDist',vtPrx'], 'delimiter', '\t');
%    Matlab
save([sPar.chFileName '.mat'],'vtDist', 'vtPathLoss', 'vtShadCorr', 'vtFading', 'vtPrx');
%
% Mostra informações do canal sintético
disp('Canal sintético:')
disp(['   Média do sombreamento: ' num2str(mean(vtShadCorr)) ]);
disp(['   Std do sombreamento: ' num2str(std(vtShadCorr)) ]);
disp(['   Janela de correlação do sombreamento: ' num2str(sPar.shadowingWindow) ' amostras' ]);
disp(['   Expoente de path loss: ' num2str(sPar.n) ]);
disp(['   m de Nakagami: ' num2str(sPar.m) ]);
%
% Plot do desvanecimento de larga escala (gráfico linear)
figure;
% Log da distância
log_distancia = log10(vtDist);
% Potência recebida com canal completo
plot(log_distancia,vtPrx); hold all;
% Potência recebida com path loss
plot(log_distancia,sPar.txPower-vtPathLoss,'linewidth', 2)
% Potência recebida com path loss e shadowing
plot(log_distancia,sPar.txPower-vtPathLoss+vtShadCorr,'linewidth', 2)
%title('Canal sintético: Potência recebida no receptor vs. log da distância')
xlabel('log_{10}(d)');
ylabel('Potência [dBm]');
legend('Prx canal completo', 'Prx (somente perda de percurso)', 'Prx (perda de percurso + sombreamento)');
xlim([0.7 1.6])
%
% Plot da geração do desvanecimento Nakagami
figure;
[f,x] = hist(vtNakagamiNormEnvelope,100);
bar(x,f/trapz(x,f)); % Histograma normalizado = PDF (das amostras)
hold all;
plot(x,fpNakaPdf(x),'r'); % PDF da distribuição
title('Canal sintético - desvanecimento de pequena escala');
legend('Histograma normalizado das amostras','PDF');