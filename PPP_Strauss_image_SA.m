clear

% Taille de la zone
Cote = 2;
S = Cote * Cote;

% Distance en lien avec la relation de voisinage
vD = 0.2;
% Parametre du processus
vGamma = 0.2; %0.1;

Nmax = 1000; % maximum number of points
points = zeros(Nmax,2);

Pb = 0.5;
Pd = 1 - Pb;

vBeta = 20;

NIte_max = 100000; % number of iterations

N_points = zeros(NIte_max,1);
Nx = 0;

N_vois = zeros(NIte_max,1);
N_v = 0;

for i = 1 : NIte_max

  R = -1;
  bDead = 0;
  N_v_temp = 0;
  Num_Point = 0;
  
  % Tirage d'un mouvement
  if(rand(1) < Pb) % Naissance d'un point
    % Tirage d'un point dans la zone
    points(Nx+1,:) = Cote * rand(1,2);
    
    % Calcul de N_v_temp

    % here...
    if(Nx > 0)
        dists = ((points(1:Nx,1) - points(Nx+1,1)).^2 + (points(1:Nx,2) - points(Nx+1,2)).^2).^0.5;
        Index = find(dists < vD);
        N_index = length(Index);
        if(N_index > 0)
            N_v_temp = length(Index);
        end
    end

    
    % Calcul de R
    % Py/Px = beta

    R = (Pd * vBeta * vGamma ^N_v_temp * S) / (Pb * (Nx + 1));
  
  else % Mort d'un point
    bDead = 1;
    if(Nx > 0)
      % Tirage d'un point dans la population
      Num_Point = ceil(Nx*rand(1));

      if(Nx > 1)
        % Echange avec le dernier point
        point_tmp = points(Nx,:);
        points(Nx,:) = points(Num_Point,:);
        points(Num_Point,:) = point_tmp;

        dists = (points(1:Nx - 1,1) - points(Nx, 1).^2 + (points(1: Nx - 1,2) - points(Nx,2)).^2);
        Index = find(dists < vD);
        N_index = length(Index);
        if (N_index > 0)
            N_v_temp = length(Index);
        end
      end
      
      % Calcul de N_v_temp
      % N_v_temp = 
      
      % Calcul de R
      R = (Pb * Nx) / (Pd * vBeta^N_v_temp * S);

    end
  end
  
  % Acceptation du mouvement
  alpha = min(1,R);
  v_unif = rand(1);
  if(v_unif < alpha) % mouvement accepte
    if(bDead)
      if(Nx > 0)
          points(Nx,:) = zeros(1,2);
          Nx = Nx - 1;
          N_v = N_v - N_v_temp;
      end
    else
      Nx = Nx + 1;
      N_v = N_v + N_v_temp;
    end
  else
     points(Nx+1,:) = zeros(1,2); 
  end

  N_points(i) = Nx;
  N_vois(i) = N_v;
  
%   figure(1)
%   clf
%   plot(points(:,1),points(:,2),'+')
end

figure(1)
clf
plot(points(:,1),points(:,2),'+')

figure(2)
subplot(2,1,1)
plot(N_points)
title('Evolution du nombre de points')

pas = 100; % pas ne change pas bq le résultat
burn = 150;
CS_N_points = cumsum(N_points(1:pas:length(N_points)));
subplot(2,1,2)
plot([1:pas:length(N_points)]', CS_N_points./[1:length(CS_N_points)]')
title('Moyenne cumulee du nombre de points')

figure(3)
subplot(2,1,1)
plot(N_vois)
title('Evolution du nombre de points voisins')

pas = 100;
CS_N_vois = cumsum(N_vois(1:pas:length(N_vois)));
subplot(2,1,2)
plot([1:pas:length(N_vois)]', CS_N_vois./[1:length(CS_N_vois)]')
title('Moyenne cumulee du nombre de points voisins')

% subplot(3,1,3)
% plot([1:length(CS_N_points)]',1.96*Ec_N_points./[1:length(CS_N_points)]' .^ 0.5);
% Index_precision = find((1.96*Ec_N_points./[1:length(CS_N_points)]').^0.5 < 0.1);

% figure(3);

ImageOrig = im2double(imread('000000.bmp'));
dim = size(ImageOrig);
N_lig = dim(1);
N_col = dim(2);
figure(1)
imagesc(ImageOrig)
axis image;
title('Image originale')
ImageOrig = ImageOrig*255; % Valeur entière entre 0 et 255.