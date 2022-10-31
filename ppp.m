clear

% Taille de la zone
Cote = 2;
S = Cote*Cote;

% Distance en lien avec la relation de voisinage
vD = 0.2;
% Parametre du processus
vGamma = 0.1;

Nmax = 1000;
points = zeros(Nmax,2);

Pb = 0.75; % probab of birth
Pd = 1 - Pb; % proba of death

vBeta = 100;

NIte_max = 100000;

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
    points(Nx+1,:) = Cote*rand(1,2);

  % Calcul de N_v_temp (numero di punti vicini)
  if(Nx > 0)
      dists = ((points(1:Nx,1) - points(Nx+1,1)).^2 + (points(1:Nx,2) - points(Nx+1,2)).^2).^0.5;
      Index = find(dists < vD);
      N_index = length(Index);
      if(N_index > 0)
          N_v_temp = length(Index);
      end
  end
    
  % Calcul de R
  R = Pd/Pb*vBeta*S/(Nx+1); %for the birth
  
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
      end
      
      % Calcul de N_v_temp
      if(Nx > 0)
        dists = ((points(1:Nx,1) - points(Nx+1,1)).^2 + (points(1:Nx,2) - points(Nx+1,2)).^2).^0.5;
        Index = find(dists < vD);
        N_index = length(Index);
        if(N_index > 0)
            N_v_temp = length(Index);
        end
    end
     
      
      % Calcul de R
      R = Pb*Nx/(Pd*vBeta*S); %for the death

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

pas = 100;
burn = 300;
CS_N_points = cumsum(N_points(1:pas:length(N_points)));

subplot(2,1,2)
plot([1:pas:length(N_points)]', CS_N_points./[1:length(CS_N_points)]');
title('Moyenne cumulee du nombre de points');

Ec_N_points = zeros(length(CS_N_points), 1);
for i = 1 : length(CS_N_points)
 if i*pas+burn < NIte_max
  Ec_N_points(i) =  std(N_points(burn:i*pas+burn));
 end
end

subplot(3,1,3)
plot([1:length(CS_N_points)]',1.96*Ec_N_points ./ [1:length(CS_N_points)]')
Index_precision = find((1.96*Ec_N_points./[1:length(CS_N_points)]') < 0.1);

figure(3)
subplot(2,1,1)
plot(N_vois)
title('Evolution du nombre de points voisins')

pas = 100;
burn = 300;
CS_N_vois = cumsum(N_vois(1:pas:length(N_vois)));
subplot(2,1,2)
plot([1:pas:length(N_vois)]', CS_N_vois./[1:length(CS_N_vois)]')
title('Moyenne cumulee du nombre de points voisins')
