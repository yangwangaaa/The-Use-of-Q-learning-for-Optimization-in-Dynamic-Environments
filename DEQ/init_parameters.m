    % number of evaluations between changes. change_frequency
    % =0 means that function never changes (or only if function change_peaks is called)
    % Scenario 1: 5000
    global change_frequency;
    change_frequency = 5000;
    %
    % number of peaks in the landscape
    % Scenario 1: 5
    %
    global number_of_peaks;
    number_of_peaks = 200;

    %
    % seed for built-in random number generator
    %
    global movrandseed;
    movrandseed = 1;

    %
    % number of dimensions, or the number of  valued genes
    % Scenerio 1:  5
    %
    global geno_size;
    geno_size = 5;

    %
    % distance by which the peaks are moved, severity
    % Scenario 1: 1.0
    %
    global vlength;
    vlength = 5;

    %
    % severity of height changes, larger numbers mean larger severity
    % Scenario 1: 7.0
    %
    global height_severity;
    height_severity = 7.0;

    %
    % severity of width changes, larger numbers mean larger severity
    %
    global width_severity;
    width_severity = 1;

    %
    % lambda determines whether there is a direction of the movement, or whether
    %  they are totally random. For lambda = 1.0 each move has the same direction,
    %  while for lambda = 0.0, each move has a random direction
    %
    global lambda;
    lambda = 0;


    %
    % if set to 1, a  landscape (basis_function) is included in the fitness
    % evaluation
    %
    global use_basis_function;
    use_basis_function = false;

    % saves computation time if not needed and set to 0 %
    global calculate_average_error;
    calculate_average_error = true; % calculate_average_error = 1;

    % saves computation time if not needed and set to 0 %
    global calculate_offline_performance;
    calculate_offline_performance = true;
    % calculate_offline_performance = 1;

    % saves computation time if not needed and set to 0 %
    global calculate_right_peak;
    calculate_right_peak = true; % calculate_right_peak = 1;

    %
    % minimum and maximum coordinate in each dimension
    % Scenario 1: 0.0 and 100.0
    %
    global mincoordinate;
    global maxcoordinate;
    mincoordinate = 0.0;
    maxcoordinate = 100.0;

    %
    % minimum and maximum height of the peaks
    % Scenario 1:  30.0 and 70.0
    %
    global minheight;
    global maxheight;
    minheight = 30.0;
    maxheight = 70.0;

    %
    % height chosen randomly when standardheight = 0.0
    % Scenario 1: 50.0
    %
    global standardheight;
    standardheight = 50.0;

    %
    % Scenario 1: 0.0001
    %
    global minwidth;
    minwidth = 1;

    %
    % Scenario 1: 0.2
    %
    global maxwidth; 
    maxwidth =12;

    %
    % width chosen randomly when standardwidth = 0.0
    % Scenario 1:  0.1
    %
    global standardwidth;
    standardwidth = 0.0;

    % Peak_Function pf = new Peak_Function1();  %  Scenario 1

    % Basis_Function bf = new Constant_Basis_Function();

    %***%END OF PARAMETER SECTION ****%

    %void change_peaks();   % preliminary declaration of function change_peaks()%
    global recent_change;
    recent_change = true; % indicates that a change has just ocurred %
    global current_peak; % peak on which the current best individual is located %
    global maximum_peak; % number of highest peak %
    global current_maximum; % fitness value of currently best individual %
    global offline_performance;
    offline_performance = 0.0;
    global offline_error;
    offline_error = 0.0;
    global avg_error;
    avg_error = 0; % average error so far %
    global current_error;
    current_error = 0; % error of the currently best individual %
    global global_max; % absolute maximum in the fitness landscape %
    global evals;
    evals = 0; % number of evaluations so far %

    % data structure to store peak data %
    global peak;
    peak=zeros(0,0); % %%peak;

    global shift;
    shift = zeros(0); % %shift;

    global coordinates;
    coordinates = zeros(0); % %coordinates;

    % which peaks are covered by the population ? %
    global covered_peaks;
    covered_peaks = zeros(0); % %covered_peaks;

    % to store every peak's previous movement %
    global prev_movement;
    prev_movement=zeros(0,0); % %%prev_movement;

    %	%
    %	 %two variables needed in method movnrand().
    %	 %
    %	 boolean backup = false;
    %	  x2;

    %
    %two variables needed in method change_stepsize_linear(). Perhaps it would be
    %appropriate to put change_stepsize_linear() in its own class.
    %
    global counter;
    counter = 1;
    global frequency;
    frequency = 3.14159 / 20.0;

    global movrand;
    global movnrand;

    global PEAKFUNCTION1
    PEAKFUNCTION1 = 0;
    global PEAKFUNCTIONCONE;
    PEAKFUNCTIONCONE = 1;
    global PEAKFUNCTIONSPHERE;
    PEAKFUNCTIONSPHERE = 2;

    global peakType;
    peakType = PEAKFUNCTIONCONE;
    %
    % Constructor
    %
    