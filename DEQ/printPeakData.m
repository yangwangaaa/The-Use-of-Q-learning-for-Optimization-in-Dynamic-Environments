    function [] = printPeakData()
    global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType;
	temp = getPeakHeights();
%	System.out.print( "Peak heights:\t" );
	for i = 0:size(temp,1)-1
%	    System.out.print( temp[ i ] + "\t" );
%	}
%	System.out.println();

%	//System.out.println("Current peak: " + current_peak + "\tMax peak: " + maximum_peak );
    end
    end
