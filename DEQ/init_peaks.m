function [] =init_peaks()
    global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType;
	shift = zeros(geno_size,1);
	coordinates = zeros(geno_size,1);
	covered_peaks = zeros(number_of_peaks,1);
	peak = zeros(number_of_peaks ,geno_size + 2);
	prev_movement = zeros(number_of_peaks,geno_size);

    peak = 100.0 * rand(number_of_peaks ,geno_size);
    prev_movement = rand(number_of_peaks ,geno_size) - 0.5;
    
	if (standardheight <= 0.0)
        peak(1:number_of_peaks,geno_size + 2) =(maxheight - minheight) * rand + minheight;
    else
        peak(1:number_of_peaks,geno_size + 2) = standardheight;
    end

	if (standardwidth <= 0.0)
        peak(1:number_of_peaks,geno_size+1) =(maxwidth - minwidth) * rand + minwidth;
    else
        peak(1:number_of_peaks,geno_size+1) = standardwidth;
    end

	if (calculate_average_error)
	    global_max = -100000.0;
	    for i = 1:number_of_peaks
            coordinates = peak(i,:);
        	dummy = dummy_eval(coordinates);
    		if (dummy > global_max)
        	    global_max = dummy;
            end
        end
    end
end
