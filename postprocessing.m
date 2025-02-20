if ~exist('ephdata', 'class')
    error('ephdata class not found. Verify file location and MATLAB path.');
end

if ~exist('lambertsolver', 'class')
    error('lambertsolver class not found. Verify file location and MATLAB path.');
end

load('roundTripSolutions.mat');
% load('departureSolutions.mat');


3function porkChopForArrivalSolution(dep_epoch, arr_epoch, asteroid_name)
    path = 'C:\Users\gadma\Documents\MATLAB\ATATD-Toolbox\All_NEOS_ATA&TD_2018_2019.csv';
    neo_data = ephdata.read_data(path);

    asteroid_id = ephdata.getAsteroidId(neo_data, asteroid_name);

    % Get the Earth's position and velocity at the departure epoch
    [r_earth, v_earth] = EphSS_car(3, dep_epoch);
    % Get the asteroid's position and velocity at the arrival epoch
    [r_ast, v_ast] = EphSS_car(asteroid_id, arr_epoch);

    range = arr_epoch - dep_epoch;

    lambertsolver.plotPorkChop(dep_epoch, arr_epoch, 3, asteroid_id, range, 1, 'Sun');
end


function porkChopForDepartureSolution(dep_epoch, arr_epoch, asteroid_name)

    path = 'C:\Users\gadma\Documents\MATLAB\ATATD-Toolbox\All_NEOS_ATA&TD_2018_2019.csv';
    neo_data = ephdata.read_data(path);

    asteroid_id = ephdata.getAsteroidId(neo_data, asteroid_name);

    % Get the asteroid's position and velocity at the departure epoch
    [r_ast, v_ast] = EphSS_car(asteroid_id, dep_epoch);
    % Get the Earth's position and velocity at the arrival epoch
    [r_earth, v_earth] = EphSS_car(3, arr_epoch);

    range = arr_epoch - dep_epoch;

    lambertsolver.plotPorkChop(dep_epoch, arr_epoch, asteroid_id, 3, range, 1, 'Sun');
end