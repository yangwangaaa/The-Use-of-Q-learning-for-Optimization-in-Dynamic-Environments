function Pso = fitness(x_in , Pso)
global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType; %#ok<NUSED>
    [l n]=size(x_in); %#ok<NASGU>
    Pso.result=zeros(l,1);
    for i=1:l
       x=x_in(i,:);
       Pso.result(i)=0;
       xxx = eval_movpeaks(x);
       Pso.result(i) = global_max - xxx;
       if (Pso.fe < Pso.max_fe)
           Pso.fe = Pso.fe + 1;
           Pso.FE(Pso.fe) = current_error;
           [Pso.fe Pso.itr Pso.FE(Pso.fe) mean(Pso.FE) Pso.test mean(Pso.offline_err_FE) std(Pso.offline_err_FE)] %#ok<NOPRT>
       else
           Pso.end = 1;
       end
    end