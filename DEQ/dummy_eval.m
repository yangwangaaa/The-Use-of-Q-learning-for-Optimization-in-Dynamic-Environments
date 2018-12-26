function [x] = dummy_eval(gen)
global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType;
	maximum = -100000.0;
	for i = 1:number_of_peaks
%         if ( peakType == PEAKFUNCTIONCONE )
%             dummy = Peak_Function_Cone(gen , i); 
%         elseif ( peakType == PEAKFUNCTIONSPHERE )
%             dummy = Peak_Function_Cone(gen , i); 
%         end
        dummy = Peak_Function_Cone(gen,i); 
        if (dummy > maximum)
            maximum = dummy;
        end
    end
    if (use_basis_function)
        dummy = Constant_Basis_Function(gen);
        if (maximum < dummy)
            maximum = dummy;
        end
    end
    x = maximum;
end
