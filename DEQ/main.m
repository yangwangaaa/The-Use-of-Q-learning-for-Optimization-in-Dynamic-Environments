close all;
clear;
clc;
format long g;
Pso.offline_err_FE=0;
for test=1 : 20
    Pso.test = test;
global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType; %#ok<*TLEV,NUSED>
init_parameters;
init_peaks;
Pso.swarm_number = 10;                                  % number of swarms
Pso.number_of_particle_in_each_swarm = 5;
Pso.best_swarm_number_of_quantum = 25;
Pso.non_best_swarm_number_of_quantum = 5;
% Pso.number_of_quantum_in_swarm = Pso.non_best_swarm_number_of_quantum * ones(1 , Pso.swarm_number);        % number of quantum in each swarm
Pso.number = Pso.swarm_number * Pso.number_of_particle_in_each_swarm;
% Pso.kind = zeros(Pso.number , 1);                       % 1 indicate neutral, 0 indicate quantum
% Pso.swarm_index = zeros(Pso.number , 1);                % swarm's index of each particle
% Pso.swarm_state = ones(1 , Pso.swarm_number)                 % 1 indicate PSO, 0 indicate QSO
% %%  initialize PSO
% for i = 1 : Pso.swarm_number
%    for j = 1 : Pso.N_num
%       Pso.kind((i-1)*Pso.number_of_particle_in_each_swarm + j) = 1;
%    end
%    for j = 1 : Pso.number_of_particle_in_each_swarm
%       Pso.swarm_index((i-1)*Pso.number_of_particle_in_each_swarm + j) = i;
%    end
% end
%%
Pso.dimension=5;
Pso.max_fe=100 * change_frequency;
Pso.FE = ones(1, Pso.max_fe);%
Pso.peak_number = number_of_peaks;
Pso.end = 0;
Pso.change_ferq = change_frequency;
Pso.fitness_function_bound=100; % [0 - 100]
Pso.max_itr=100000;
Pso.Pbest_value = inf(Pso.number , 1);
Pso.Pbest_position = zeros(Pso.number , Pso.dimension);
Pso.FE = 0;  %
Pso.fe = 0; %
Pso.velocity = rand(Pso.number , Pso.dimension);
Pso.change_accur = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pso.cloud = 0.5 * vlength * ones(Pso.swarm_number , Pso.dimension);                      % r_cloud = severity
Pso.excl = 0.5 * (100 / (Pso.peak_number ^ (1 / Pso.dimension)));   % r_excl = 0.5*d_boa 
Pso.conv = 40; %40;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pso.itr = 0;
Pso.X = Pso.fitness_function_bound * (rand(Pso.number , Pso.dimension));
Pso.particle_value = zeros(Pso.number , 1);
Pso = fitness(Pso.X , Pso);
Pso.particle_value =Pso.result;
Pso.Pbest_value = Pso.particle_value;
Pso.Pbest_position = Pso.X;
Pso.Gbest_value = inf(1,Pso.swarm_number);
Pso.Gbest_position = inf(Pso.swarm_number , Pso.dimension);
for i=1:Pso.swarm_number
   [val id] = min(Pso.Pbest_value(((i-1)*Pso.number_of_particle_in_each_swarm + 1) : (i*Pso.number_of_particle_in_each_swarm) ));
   id = id + ((i-1)*Pso.number_of_particle_in_each_swarm);
   Pso.Gbest_value(i) = val;
   Pso.Gbest_position(i,:) = Pso.Pbest_position(id,:);
end
Pso.Gbest_past_value  = Pso.Gbest_value;
Pso.is_converged = zeros(1 , Pso.swarm_number);
Pso.randomized = zeros(1 , Pso.swarm_number);
Pso.x = 0.729843788;
c1 = 2.05;
c2 = 2.05;
Pso.randomized = zeros(1 , Pso.swarm_number);
Pso.is_converged = zeros(1 , Pso.swarm_number);
Pso.itr=0;
Pso.success_rate = zeros(1 , Pso.swarm_number);
Pso.weight_min = 0.6;
Pso.weight_max = 1;
%%

% figure(1); hold on;
% [X,Y] = meshgrid(0:1:100);
% 
% Z = (peak(1,end) - peak(1,end-1)*(sqrt((X-peak(1,1)).^2 + (Y-peak(1,2)).^2)));
% i = 2;
% while i <= number_of_peaks
%   ZI =(peak(i,end) - peak(i,end-1)*(sqrt((X-peak(i,1)).^2 + (Y-peak(i,2)).^2)));
%   Z = max(Z,ZI);  
%   i = i+1;
% end
% 
% meshc(Z);
% 
% figure(2);
while(1==1)

%     Z = (peak(1,end) - peak(1,end-1)*(sqrt((X-peak(1,1)).^2 + (Y-peak(1,2)).^2)));
%     kk = 2;
%     while kk <= number_of_peaks
%       ZI =(peak(kk,end) - peak(kk,end-1)*(sqrt((X-peak(kk,1)).^2 + (Y-peak(kk,2)).^2)));
%       Z = max(Z,ZI);  
%       kk = kk + 1;
%     end
%    contour(X,Y,Z,20); hold on;    
%    
%    plot(Pso.Gbest_position(:,1),Pso.Gbest_position(:,2),'.','markersize',10,'markerfacecolor','g');
% %    
    Pso.itr = Pso.itr + 1;
%%

% %%%%%%%%%%%% Transform BEGIN
% if Pso.itr > 2
%     for ii=1 : Pso.swarm_number
%         clear tmp;
%         if Pso.swarm_state(ii == 1)
%             if ((Pso.Gbest_past_value(ii) - Pso.Gbest_value(ii)) < Pso.tarnsform_condition) && ((Pso.Gbest_past_value(ii) - Pso.Gbest_value(ii)) >= 0)
%                 Pso.swarm_state(ii) = 0;
%             end
%         end
%     end
% end
% %%%%%%%%%%%% Transform END

%%

%%%%%%%%%%%% Anti Convergence BEGIN
if(Pso.peak_number > Pso.swarm_number)% if p>m
    for ii=1 : Pso.swarm_number
       clear tmp;
       if Pso.is_converged(ii) == 0
          for jj=1 : Pso.number_of_particle_in_each_swarm
              for kk=1 : Pso.number_of_particle_in_each_swarm
                  tmp(ii , jj) = (sqrt( sum((Pso.X(((ii-1) * Pso.number_of_particle_in_each_swarm) + jj ,:)-Pso.X(((ii-1) * Pso.number_of_particle_in_each_swarm) + kk ,:)).^2))) > Pso.conv ; %#ok<SAGROW>
              end
          end
          if (sum(sum(tmp)) == 0)
              Pso.is_converged(ii) = 1;
%               Pso.swarm_state(ii) = 0;
          end
       end
    end
    if (sum(Pso.is_converged) == Pso.swarm_number) && (Pso.swarm_number < Pso.peak_number)
        [val id] = max(Pso.Gbest_value);
        Pso.randomized(id) = 1;
    end
end
%%%%%%%%%%%% Anti Convergence END
%%
%%%%%%%%%%%% Exclusion BEGIN
    for i=1:Pso.swarm_number
        for j=1:Pso.swarm_number
            if (i ~= j) && (sqrt( sum( (Pso.Gbest_position(i,:) - Pso.Gbest_position(j,:)).^2 ) ) <= Pso.excl)
                if (Pso.Gbest_value(i) < Pso.Gbest_value(j))
                    Pso.randomized(j) = 1;
                end
            end
        end
    end
%%%%%%%%%%%% Exclusion END
%%
%%%%%%%%%%%% Test for change 
Pso = fitness(Pso.Gbest_position , Pso);
temp = Pso.result;
for i=1 : Pso.swarm_number
    if (temp(i) ~= Pso.Gbest_value)
        Pso.change_accur = 1;
        break;
    end
end
if (Pso.change_accur == 1)
   Pso.cloud = 0.5 * vlength * ones(Pso.swarm_number , Pso.dimension);
   Pso.randomized = zeros(1 , Pso.swarm_number);
   for i=1 : Pso.swarm_number
           Pso = fitness(Pso.Pbest_position(( ((i-1)*Pso.number_of_particle_in_each_swarm + 1) : (i*Pso.number_of_particle_in_each_swarm) ) , :) , Pso);
           Pso.Pbest_value((((i-1)*Pso.number_of_particle_in_each_swarm + 1) : (i*Pso.number_of_particle_in_each_swarm))) = Pso.result;
   end
   for i=1:Pso.swarm_number
          [val id] = min(Pso.Pbest_value( ((i-1)*Pso.number_of_particle_in_each_swarm + 1) : (i*Pso.number_of_particle_in_each_swarm) ));
          id = id + ((i-1)*Pso.number_of_particle_in_each_swarm);
          Pso.Gbest_value(i) = val;
          Pso.Gbest_position(i,:) = Pso.Pbest_position(id,:);
   end
   Pso.change_accur = 0;
   tmp = rem(Pso.fe , Pso.change_ferq);
   temp = min(Pso.Pbest_value);
   for ii=1 : tmp
       Pso.FE(Pso.fe-ii+1) = temp;
   end
end
%%%%%%%%%%%% Change Accur END
%%
%%%%%%%%%%%% Execute Randomization BEGIN
for i=1 : Pso.swarm_number
    if (Pso.randomized(i) == 1)
        Pso.randomized(i) = 0;
        Pso.cloud (i ,:) = 0.5 * vlength * ones(1 , Pso.dimension);
        Pso.is_converged(i) = 0;
        Pso.X((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm) , :) = Pso.fitness_function_bound * (rand(Pso.number_of_particle_in_each_swarm , Pso.dimension));
        Pso = fitness(Pso.X((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm) , :) , Pso);
        Pso.particle_value((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm)) = Pso.result;
        Pso.Pbest_value((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm)) = Pso.particle_value((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm));
        Pso.Pbest_position((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm) , :) = Pso.X((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm) , :);
        [val id] = min(Pso.Pbest_value((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm)));
        Pso.Gbest_value(i) = val;
        Pso.Gbest_position(i,:) = Pso.Pbest_position((((i-1) * Pso.number_of_particle_in_each_swarm) + id), :);
    end
end
%%%%%%%%%%%% Execute Randomization END
%%
%%%%%%%%%%%% Update particles based on particles type
for i=1 : Pso.swarm_number
        for j=1:Pso.number_of_particle_in_each_swarm
               Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , :) = ...
                  Pso.x* ...
                  (Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , :) + ...
                  (c1*rand(1 , Pso.dimension).*(Pso.Pbest_position((i-1)*Pso.number_of_particle_in_each_swarm + j , :) - Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j , :))) + ...
                  (c2*rand(1 , Pso.dimension).*(Pso.Gbest_position(i,:) - Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j , :))));
              for k=1:Pso.dimension
                  if Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j , k) > 100
                     Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , k) = 0;
                     Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j , k) = 100;
                  elseif Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j , k) < 0
                     Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , k) = 0;
                     Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j , k) = 0;
                  end
              end
              for k=1:Pso.dimension
                  if Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , k) > 50
                     Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , k) = 50;
                  elseif Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , k) < -50
                     Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , k) = -50;
                  end
              end

              Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j,:) = Pso.X((i-1)*Pso.number_of_particle_in_each_swarm + j , :) + ...
                  Pso.velocity((i-1)*Pso.number_of_particle_in_each_swarm + j , :);
              Pso = fitness(Pso.X((((i-1) * Pso.number_of_particle_in_each_swarm) + j)  , :) , Pso);
              Pso.particle_value((((i-1) * Pso.number_of_particle_in_each_swarm) + j) ) = Pso.result;
              

        end
        for j=1 : Pso.number_of_particle_in_each_swarm
            if (Pso.particle_value(((i-1) * Pso.number_of_particle_in_each_swarm) + j) < Pso.Pbest_value(((i-1) * Pso.number_of_particle_in_each_swarm) + j))
                Pso.Pbest_value(((i-1) * Pso.number_of_particle_in_each_swarm) + j) = Pso.particle_value(((i-1) * Pso.number_of_particle_in_each_swarm) + j);
                Pso.Pbest_position((((i-1) * Pso.number_of_particle_in_each_swarm) + j) , :) = Pso.X((((i-1) * Pso.number_of_particle_in_each_swarm) + j) , :);
            end
        end
        [val id] = min(Pso.Pbest_value((((i-1) * Pso.number_of_particle_in_each_swarm) + 1) : (i * Pso.number_of_particle_in_each_swarm)));
        Pso.Gbest_value(i) = val;
        Pso.Gbest_position(i,:) = Pso.Pbest_position((((i-1) * Pso.number_of_particle_in_each_swarm) + id), :);
end
[notting best_swarm_id] = min(Pso.Gbest_value);
for i=1 : Pso.swarm_number
        if i  ~= best_swarm_id
           Pso.success_rate = 0;
           for j = 1 : Pso.non_best_swarm_number_of_quantum
               temp = Pso.Gbest_position(i,:) +  Pso.cloud(i , :) .* rands(1,Pso.dimension);   
               Pso = fitness(temp , Pso);
               temp_value = Pso.result;
               if temp_value <= Pso.Gbest_value(i)
                   Pso.Gbest_position(i,:) = temp;
                   Pso.Gbest_value(i) = temp_value;
                   Pso.success_rate = Pso.success_rate + 1;
               end
           end
           Pso.success_rate = Pso.success_rate / Pso.non_best_swarm_number_of_quantum;           
        elseif i == best_swarm_id
           Pso.success_rate = 0;
           for j = 1 : Pso.best_swarm_number_of_quantum
               temp = Pso.Gbest_position(i,:) +  Pso.cloud(i , :) .* rands(1,Pso.dimension);   
               Pso = fitness(temp , Pso);
               temp_value = Pso.result;
               if temp_value <= Pso.Gbest_value(i)
                   Pso.Gbest_position(i,:) = temp;
                   Pso.Gbest_value(i) = temp_value;
                   Pso.success_rate = Pso.success_rate + 1;                
               end
           end
           Pso.success_rate = Pso.success_rate / Pso.best_swarm_number_of_quantum;                     
        end
   Pso.weight = Pso.weight_min + (Pso.weight_max - Pso.weight_min) * Pso.success_rate;
   Pso.cloud(i , :) = Pso.cloud(i , :) .* Pso.weight;
end
    if (Pso.end == 1)
        break;
    end
%     drawnow;
%     hold off;        
end
Pso.offline_err_FE(test) = mean(Pso.FE); 
current_error_FE(test , :) = Pso.FE; %#ok<SAGROW>
end
Min = min(Pso.offline_err_FE) %#ok<*NOPTS>
Max = max(Pso.offline_err_FE)
Mean = mean(Pso.offline_err_FE)
Std = std(Pso.offline_err_FE)
XX = mean(current_error_FE);
semilogy(XX);
% 
% Offline_error = mean(Pso.FE) %#ok<NOPTS>
% %semilogy(Pso.FE);