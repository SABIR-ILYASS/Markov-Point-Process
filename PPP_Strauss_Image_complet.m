close all;
clear all;
clc;

% initialisation générateur aléatoire
rand('twister',sum(100*clock))

ImageOrig = im2double(imread('000000.bmp'));
dim = size(ImageOrig);
N_lig = dim(1);
N_col = dim(2);

figure(1)
imagesc(ImageOrig)
axis image;
title('Image originale')

ImageOrig = ImageOrig*255; % Valeur entière entre 0 et 255.

%% Taille de la zone
Cote = N_lig;
S = N_lig*N_col;

% Distance en lien avec la relation de voisinage
vD = 10;
% Parametre du processus
vGamma = 0.5;

% Parametre attache aux données
%Canal = 1;
Seuil = 70;
vData = 0.01;

% Paramètre d'intensité
vBeta = 0.02;

Nmax = 1000; % Nombre maximum de points
points = zeros(Nmax,2); % vecteur contenant les points d'une population
Data_term = zeros(Nmax,1); % vecteur contenant les valeurs du terme d'attache aux données

% Probabilités des mouvements de naissance et de mort d'un point
Pb = 0.5;
Pd = 1 - Pb;

NIte_max = 100000;

% Vecteur qui mémorise le nombre de points itération par itération
N_points = zeros(NIte_max,1);
Nx = 0;

% Vecteur qui mémorise le nombre de points en relation de voisinage itération par itération
N_vois = zeros(NIte_max,1);
N_v = 0;

% Vecteur qui mémorise la proportion de points bien placés
N_Data = zeros(NIte_max,1);

for i = 1 : NIte_max

  if(mod(i,100) == 0)
      i
  end
  
  R = -1;
  bDead = 0;
  N_v_temp = 0;
  Num_Point = 0;
  
  % Tirage d'un mouvement
  if(rand(1) < Pb) % Naissance d'un point
    % Tirage d'un point dans la zone
    points(Nx+1,:) = (Cote-2)*rand(1,2)+1;
    
    % Calcul du terme d'attache aux données
    Pix_hg = floor(points(Nx+1,:));
    bMal_Pos = 0;
    Canal = 1;
    Val_images1 = [ImageOrig(Pix_hg(1),Pix_hg(2),Canal) ImageOrig(Pix_hg(1)+1,Pix_hg(2),Canal) ImageOrig(Pix_hg(1),Pix_hg(2)+1,Canal) ImageOrig(Pix_hg(1)+1,Pix_hg(2)+1,Canal)];
    Canal = 2;
    Val_images2 = [ImageOrig(Pix_hg(1),Pix_hg(2),Canal) ImageOrig(Pix_hg(1)+1,Pix_hg(2),Canal) ImageOrig(Pix_hg(1),Pix_hg(2)+1,Canal) ImageOrig(Pix_hg(1)+1,Pix_hg(2)+1,Canal)]; 
    Data_term(Nx+1) = (mean(Val_images1)+mean(Val_images2))/2;
    if(Data_term(Nx+1) > Seuil)
        bMal_Pos = 1;
    end
    
    % Calcul du nombre de points voisins
    if(Nx > 0)
      dists = ((points(1:Nx,1) - points(Nx+1,1)).^2 + (points(1:Nx,2) - points(Nx+1,2)).^2).^0.5;
      Index = find(dists < vD);
      N_index = length(Index);
      if(N_index > 0)
        N_v_temp = length(Index);
      end
    end
    
    R = vBeta*(vGamma^N_v_temp)*Pd/Pb*S/(Nx+1);
    if(bMal_Pos)
        R = R*vData;
    end
    
  else % Mort d'un point
    bDead = 1;
    if(Nx > 0)
      % Tirage d'un point dans la population
      Num_Point = ceil(Nx*rand(1));
      N_v_temp = 0;
      if(Nx > 1)
        % Echange avec le dernier point
        point_tmp = points(Nx,:);
        Data_term_tmp = Data_term(Nx);
        points(Nx,:) = points(Num_Point,:);
        Data_term(Nx) = Data_term(Num_Point);
        points(Num_Point,:) = point_tmp;
        Data_term(Num_Point) = Data_term_tmp;
        % Calcul du nombre de points voisins
        dists = ((points(1:Nx-1,1) - points(Nx,1)).^2 + (points(1:Nx-1,2) - points(Nx,2)).^2).^0.5;
        Index = find(dists < vD);
        N_index = length(Index);
        if(N_index > 0)
            N_v_temp = length(Index);
        end
      end
      % Terme d'attache aux données
      bMal_Pos = 0;
      if(Data_term(Nx) > Seuil)
          bMal_Pos = 1;
      end
      
      R = 1/(vBeta*(vGamma^N_v_temp))*Pb/Pd*Nx/S;
      if(bMal_Pos)
          R = R/vData;
      end
      
    end
  end
  
  % Acceptation du mouvement
  alpha = min(1,R);
  v_unif = rand(1);
  if(v_unif < alpha) % mouvement accepte
    if(bDead)
      if(Nx > 0)
          points(Nx,:) = zeros(1,2);
          Data_term(Nx) = 0;
          Nx = Nx - 1;
          N_v = N_v - N_v_temp;
      end
    else
      Nx = Nx + 1;
      N_v = N_v + N_v_temp;
    end
  else
     points(Nx+1,:) = zeros(1,2);
     Data_term(Nx+1) = 0;
  end

  N_points(i) = Nx;
  N_vois(i) = N_v;
  if(Nx > 0)
    N_Data(i) = 100*length(find(Data_term(1:Nx) <= Seuil))/Nx;
  end
  
%   figure(1)
%   clf
%   plot(points(:,1),points(:,2),'+')
end

figure(2)
clf
plot(points(:,1),points(:,2),'+')

figure(3)
subplot(2,1,1)
plot(N_points)
title('Evolution du nombre de points')

Burn = 4000;
pas = 50;
CS_N_points = cumsum(N_points(Burn:pas:length(N_points)));
subplot(2,1,2)
plot([Burn:pas:length(N_points)]', CS_N_points./[1:length(CS_N_points)]')
title('Moyenne cumulée du nombre de points')

figure(4)
subplot(2,1,1)
plot(N_vois)
title('Evolution du nombre de points voisins')

CS_N_vois = cumsum(N_vois(Burn:pas:length(N_vois)));
subplot(2,1,2)
plot([Burn:pas:length(N_vois)]', CS_N_vois./[1:length(CS_N_vois)]')
title('Moyenne cumulée du nombre de points voisins')

CS2_N_points = cumsum(N_points(Burn:pas:length(N_points)).^2);
diff = CS2_N_points./[1:length(CS2_N_points)]' - (CS_N_points./[1:length(CS_N_points)]').^2;
Borne_sup = 1.96*diff.^(0.5)./[1:length(CS_N_points)].^0.5';

figure(5)
plot([Burn:pas:length(N_vois)]',Borne_sup)
title('Borne supérieure intervalle de confiance estimation du nombre de points moyens')

figure(6)
subplot(2,1,1)
plot(N_Data)
title('Evolution de la proportion de points bien placés')

CS_N_Data = cumsum(N_Data(Burn:pas:length(N_vois)));
subplot(2,1,2)
plot([Burn:pas:length(N_Data)]', CS_N_Data./[1:length(CS_N_Data)]')
title('Moyenne cumulée')

% Data_term2 = zeros(Nmax,1); % vecteur contenant les valeurs du terme d'attache aux données
%                             % Vérification

for i = 1 : Nx
    Pix_hg = floor(points(i,:));
    bMal_Pos = 0;
    val = 255;
    if(Data_term(i) > Seuil)
        bMal_Pos = 1;
        val = 0;
    end
%     Val_images = [ImageOrig(Pix_hg(1),Pix_hg(2),Canal) ImageOrig(Pix_hg(1)+1,Pix_hg(2),Canal) ImageOrig(Pix_hg(1),Pix_hg(2)+1,Canal) ImageOrig(Pix_hg(1)+1,Pix_hg(2)+1,Canal)]; 
%     Data_term2(i) = mean(Val_images);
    ImageOrig(Pix_hg(1),Pix_hg(2),:) = [val 0 0];
    ImageOrig(Pix_hg(1),Pix_hg(2)+1,:) = [val 0 0];
    ImageOrig(Pix_hg(1)+1,Pix_hg(2),:) = [val 0 0];
    ImageOrig(Pix_hg(1)+1,Pix_hg(2)+1,:) = [val 0 0];
end

figure(7)
imagesc(ImageOrig/255)
axis image;
title('Image originale')





