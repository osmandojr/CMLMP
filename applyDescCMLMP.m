%---------------------------------------------------------------
% CMLMP - Completed Mean Local Mapped Pattern 
%---
% # developed by Osmando Pereira Junior, 2019
% # mail: osmandoj@gmail.com
% # based on MLMP and CLMP descriptors
%---
% Input arguments:
% image - Gray scale image.
% radius - distance from neighborhood to the center pixel: radius fixed and
% equals to 2;
% numNeighbors - number of sample points in the neighborhood circle:
% numNeighbors fixed and equals to 8.
% beta - beta parameter for each CLMP componente: [beta_S, beta_M, beta_C]
% bin_S, bin_M, bin_C - number of bins of the local cell histogram; parameter b, that
% defines the histogram size for Signal, Magnitude and Center components.
%---------------------------------------------------------------

function [MLMP_S,MLMP_M,MLMP_C,t_desc] =...
            applyDescCMLMP_v02(image,radius,numNeighbors,...
                           beta_S,beta_M,beta_C,...
                           bin_S,bin_M,bin_C)

%---
% Declarar e inicializar variavel tempo
t_desc = zeros(1,3);
t1 = tic;
%---

%---
% Verificar versao do descritor a se considerar: invariante ou nao aa rotacao
% descMode = {'ri','classic'}
p = ones(1,numNeighbors);
for i = 1 : numNeighbors
    p(i) = 2^(i-1);
end
%---

%---
% Garantir que a imagem estah no formato double
image = double(image);
%---

%---
% Declarar e inicializar vizinhanca circular 
spoints = zeros(numNeighbors,2);
%
% Determinar passo angular para vizinhanca circular
a = 2*pi/numNeighbors;
%
% Determinar posicao de cada um dos N vizinhos
for i = 1 : numNeighbors
    spoints(i,1) = -radius*sin((i-1)*a);
    spoints(i,2) = radius*cos((i-1)*a);
end 
%---
    
%---
% Recuperar a dimensao da imagem a ser codificada
[ysize,xsize] = size(image);
%---

%---
% Determinacao do bloco limitante da vizinhanca para calculo da
% caracteristica MLMP
% Cada caracteristica MLMP eh determinada considerando-se um bloco de
% dimensao bsizey x bsizex, em que o pixel em questao ocupa a posicao
% central do bloco
miny = min(spoints(:,1));
maxy = max(spoints(:,1));
minx = min(spoints(:,2));
maxx = max(spoints(:,2));
%
% Block size
bsizey = ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex = ceil(max(maxx,0))-floor(min(minx,0))+1;
%
% Coordinates of origin (0,0) in the block
origy = 1-floor(min(miny,0));
origx = 1-floor(min(minx,0));
%---

%---
% Verificacao se imagem tem dimensao suficiente para aplicacao do descritor
% A imagem deve ser no minimo igual ao tamanho do bloco
if(xsize < bsizex || ysize < bsizey)
    error('Imagem muito pequena. Deve ter pelo menos (2*radius+1) x (2*radius+1).');
end
%---

%---
% Determinar janela da imagem para a qual serao determinadas as
% caracteristicas MLMP
% dx e dy eh a dimensao da imagem desconsiderando as bordas
dx = xsize - bsizex;
dy = ysize - bsizey;
%---

%---
% Recuperar nivel de cinza dos pixels centrais, pixels para os quais serao
% determinadas as caracteristicas MLMP
C = image(origy:origy+dy,origx:origx+dx);
%---

%---
% Declaracao e inicializacao da variavel que contera os niveis de cinza dos
% pixels da vizinhanca
N = zeros(size(C,1),size(C,2),numNeighbors);
%---
% Enlace sobre a vizinhanca para recuperacao do nivel de cinza dos N pixels
% da vizinhanca para todos os pixels centrais
for i = 1 : numNeighbors
    y = spoints(i,1)+origy;
    x = spoints(i,2)+origx;

    % Calculate floors, ceils and rounds for the x and y.
    fy = floor(y); 
    cy = ceil(y); 
    ry = round(y);
    fx = floor(x); 
    cx = ceil(x); 
    rx = round(x);

    %---
    % Verificar se eh necessario interpolacao
    if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
        % Interpolation is not needed, use original datatypes
        N(:,:,i) = image(ry:ry+dy,rx:rx+dx);
        %NArray(:,i) = reshape(N,prod(size(N)),1);
        %D = N >= C; 
    else
        % Interpolation needed, use double type images 
        ty = y - fy;
        tx = x - fx;

        % Calculate the interpolation weights.
        w1 = (1 - tx) * (1 - ty);
        w2 =      tx  * (1 - ty);
        w3 = (1 - tx) *      ty ;
        w4 =      tx  *      ty ;

        % Compute interpolated pixel values
        N(:,:,i) = w1*image(fy:fy+dy,fx:fx+dx) +...
                   w2*image(fy:fy+dy,cx:cx+dx) +...
                   w3*image(cy:cy+dy,fx:fx+dx) +...
                   w4*image(cy:cy+dy,cx:cx+dx);
        %---
    end
    % Fim verificao da necessidade de interpolacao
    %---
end
% Fim enlace sobre a vizinhanca para recuperacao do nivel de cinza dos N 
% pixels da vizinhanca para todos os pixels centrais
%---

%---
% Determinacao da media dos pixels da vizinhanca
% N_mean = (N(:,:,1)+N(:,:,2)+...+N(:,:,numNeighbours)/numNeighbours;
N_mean = mean(N,3);
%---

%---
% Determine the MLMP feature for each pixel
%---

%---
% Determine the MLMP S, M, and C features for each central pixel
%---
% 1. MLMP_S
%---
pert = zeros(size(N));
for i = 1 : numNeighbors
    pert(:,:,i) = (1./(1 + exp(-(N(:,:,i) - N_mean) / beta_S)))*p(i);
end
pert(pert<0) = 0;
pert = sum(pert,3)./sum(p);
MLMP_S = round(pert * (bin_S-1));
MLMP_S = hist(MLMP_S(:),bin_S);
t_desc(1) = toc(t1);
%
%---
% 2. MLMP_M
%---
% Determine the average magnitude
M = abs(N - N_mean);
M_mean = mean(image(:));
pert = zeros(size(M));
for i  = 1 : numNeighbors
    pert(:,:,i) = (1./(1 + exp(-(M(:,:,i) - M_mean)/beta_M)))*p(i);
end
pert(pert<0) = 0;
pert = sum(pert,3)./sum(p);
MLMP_M = round(pert * (bin_M-1));
MLMP_M = hist(MLMP_M(:),bin_M);
%
%---
% 3. MLMP_C
%---  
pert = 1./(1 + exp(-(C - mean(image(:))) / beta_C));
pert(pert<0) = 0;
MLMP_C = round(pert * (bin_C-1));
MLMP_C = hist(MLMP_C(:),bin_C);
%---

