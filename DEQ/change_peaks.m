function [] = change_peaks()
    global change_frequency,global movrandseed,global geno_size,global vlength,global height_severity,global width_severity,global lambda,global number_of_peaks,global use_basis_function,global calculate_average_error,global calculate_offline_performance,global calculate_right_peak,global mincoordinate,global maxcoordinate,global minheight,global maxheight,global standardheight,global minwidth,global maxwidth, global standardwidth,global recent_change,global current_peak, global maximum_peak, global current_maximum, global offline_performance,global offline_error,global avg_error,global current_error,global global_max,global evals,global peak,global shift,global coordinates,global covered_peaks,global prev_movement,global counter,global frequency,global movrand,global movnrand,global PEAKFUNCTION1,global PEAKFUNCTIONCONE,global PEAKFUNCTIONSPHERE,global peakType;
    for i = 1:number_of_peaks
	    % shift peak locations
	    sum = 0.0;
	    for j = 1:geno_size
            shift(j) = rand - 0.5;
            sum = sum + shift(j) * shift(j);
        end
	    if (sum > 0.0)
    		sum = vlength / sqrt(sum);
        else
        % only in case of rounding errors */
    		sum = 0.0;
        end
	    sum2 = 0.0;
	    for j = 1:geno_size
            shift(j) = sum * (1.0 - lambda) * shift(j) + lambda * prev_movement(i,j);
            sum2 = sum2 + shift(j) * shift(j);
        end
	    if (sum2 > 0.0)
    		sum2 = vlength / sqrt(sum2);
	    else % only in case of rounding errors */
    		sum2 = 0.0;
        end;
	    for j = 1:geno_size
    		shift(j) = shift(j) * sum2;
        	prev_movement(i,j) = shift(j);
            if ((peak(i,j) + prev_movement(i,j))< mincoordinate) 
                peak(i,j) = 2.0 * mincoordinate - peak(i,j) - prev_movement(i,j);
                prev_movement(i,j) = prev_movement(i,j) * -1.0;
            elseif ((peak(i,j) + prev_movement(i,j)) > maxcoordinate)
                peak(i,j) =	2.0 * maxcoordinate	- peak(i,j)	- prev_movement(i,j);
                prev_movement(i,j) = prev_movement(i,j) * -1.0;
            else
                peak(i,j) = peak(i,j) + prev_movement(i,j);
            end
            %/* change peak width */
        end
	    j = geno_size+1;
	    offset = randn * width_severity;
	    if ((peak(i,j) + offset) < minwidth)
    		peak(i,j) = 2.0 * minwidth - peak(i,j) - offset;
        elseif ((peak(i,j) + offset) > maxwidth)
        	peak(i,j) = 2.0 * maxwidth - peak(i,j) - offset;
        else
            peak(i,j) = peak(i,j) + offset;
        end;
	    %/* change peak height */
	    j=j+1;
	    offset = height_severity * randn;
	    if ((peak(i,j) + offset) < minheight)
    		peak(i,j) = 2.0 * minheight - peak(i,j) - offset;
        elseif ((peak(i,j) + offset) > maxheight)
        	peak(i,j) = 2.0 * maxheight - peak(i,j) - offset;
	    else
		peak(i,j) = peak(i,j) + offset;
        end
    end
	if (calculate_average_error) 
	    global_max = -100000.0;
	    for i = 1: number_of_peaks
		    coordinates = peak(i,1:geno_size);
            dummy = dummy_eval(coordinates);
            if (dummy > global_max) 
                global_max = dummy;
                maximum_peak = i;
            end
        end
    end
	recent_change = true;
end
