clc; clear; close all;

% Check if class exists
if ~exist('ephdata', 'class')
    error('ephdata class not found. Verify file location and MATLAB path.');
end

if ~exist('lambertsolver', 'class')
    error('lambertsolver class not found. Verify file location and MATLAB path.');
end


mainRoundTrip();


function mainRoundTrip()
    % 1) Earth->Ast solutions
    if isempty(gcp('nocreate'))
        parpool('local', 6); 
    end
    outboundSolutions = departureSolutions();  % writes to 'AsteroidSolutions.mat'

    % 2) If we want the round-trip
    if ~isempty(outboundSolutions)
        roundTrips = arrivalSolutions(outboundSolutions);

        if ~isempty(roundTrips)
            roundTripTable = struct2table(roundTrips);
            % Maybe sort by total mission time or final arrival date:
            roundTripTable = sortrows(roundTripTable, 'TotalMissionDays');
            disp(roundTripTable);

            save('RoundTripSolutions.mat','roundTrips');
        else
            disp('No round-trip solutions found given the inbound solutions!');
        end
    end
end


function roundTripSolutions = arrivalSolutions(solutions)
    % ARRIVALSOLUTIONS
    %   For each Earth->Asteroid solution in `solutions`, find a feasible
    %   Asteroid->Earth return starting from the same day we arrive
    %   (zero dwell) up to final mission day (12/31/2037).
    %
    %   Stores the short/long arc velocity vectors in the output struct.

    % We'll store round-trip solutions in a struct array, but using a cell array
    % to collect partial results in parallel.
    roundTripSolutions = struct('AstName',       {}, ...
                                'AstID',         {}, ...
                                'DepartEarth',   {}, ...
                                'ArriveAst',     {}, ...
                                'DepartAst',     {}, ...
                                'ArriveEarth',   {}, ...
                                'ArcTypeReturn', {}, ...
                                'TOF_daysReturn',{}, ...
                                'vInfReturn',    {}, ...
                                'dvAstDep',      {}, ...
                                'TotalMissionDays', {}, ...
                                ... % NEW FIELDS for velocity vectors:
                                'v1DepartVec',   {}, ...
                                'v2ArriveVec',   {});

    % If there are no inbound solutions, nothing to do
    if isempty(solutions)
        fprintf('No outbound solutions to process.\n');
        return;
    end

    finalMissionDate = date2mjd2000([2037 12 31 0 0 0]);
    step  = 5;  % smaller step for the return search?

    % Pre-allocate a cell array where each element will hold
    % the solutions found for that single inbound solution
    nSol = length(solutions);
    solutionsCell = cell(nSol,1);

    % PARFOR loop over each feasible Earth->Asteroid solution
    parfor iSol = 1:nSol
        astID      = solutions(iSol).AstID;
        arrivalAst = solutions(iSol).Arrival; % MJD2000 day we arrive at the asteroid

        % Local array to store the solutions from this single iteration
        localSolutions = [];

        % If we arrive after finalMissionDate, no time left for return
        if arrivalAst >= finalMissionDate
            solutionsCell{iSol} = localSolutions;
            continue;
        end

        % Build a +/- range so that earliest departure is arrivalAst,
        % and latest arrival is finalMissionDate
        range = finalMissionDate - arrivalAst;

        % 1) We call the updated solver for the second leg: ASTEROID -> EARTH
        [vInfShortRet, dvShortRet, vInfLongRet, dvLongRet, ...
        DepGridRet, ArrGridRet, ...
        v1ShortCell, v2ShortCell, v1LongCell, v2LongCell] = ...
            lambertsolver.findTransferSolutions( ...
                arrivalAst,        ...  % center date for departure
                finalMissionDate,  ...  % center date for arrival
                astID,             ...  % depart from the asteroid
                3,                 ...  % arrive at Earth
                range, step, 'sun');

        % 2) Constraints: dvAstDep < 0.5, vInfEarth < 1.5
        returnShortMask = (vInfShortRet < 1.5) & (dvShortRet < 0.5);
        returnLongMask  = (vInfLongRet  < 1.5) & (dvLongRet  < 0.5);

        % 3) Short-Arc Feasible Solutions
        [iD_rs, iA_rs] = find(returnShortMask);
        for k = 1:numel(iD_rs)
            idxD = iD_rs(k);
            idxA = iA_rs(k);

            tDepAst   = DepGridRet(idxD);
            tArrEarth = ArrGridRet(idxA);

            % Must not depart asteroid before arrival, and must finish by final day
            if tDepAst < arrivalAst,          continue; end
            if tArrEarth > finalMissionDate,  continue; end

            totalMissionDays = tArrEarth - solutions(iSol).Departure;

            entry = struct( ...
                'AstName',          solutions(iSol).AstName, ...
                'AstID',            astID, ...
                'DepartEarth',      solutions(iSol).Departure, ...
                'ArriveAst',        arrivalAst, ...
                'DepartAst',        tDepAst, ...
                'ArriveEarth',      tArrEarth, ...
                'ArcTypeReturn',    "Short", ...
                'TOF_daysReturn',   tArrEarth - tDepAst, ...
                'vInfReturn',       vInfShortRet(idxD, idxA), ...
                'dvAstDep',         dvShortRet(idxD, idxA), ...
                'TotalMissionDays', totalMissionDays, ...
                ... % Store the actual velocity vectors for the return
                'v1DepartVec',      v1ShortCell{idxD, idxA}, ...
                'v2ArriveVec',      v2ShortCell{idxD, idxA} ...
                );

            localSolutions = [localSolutions; entry];
        end

        % 4) Long-Arc Feasible Solutions
        [iD_rl, iA_rl] = find(returnLongMask);
        for k = 1:numel(iD_rl)
            idxD = iD_rl(k);
            idxA = iA_rl(k);

            tDepAst   = DepGridRet(idxD);
            tArrEarth = ArrGridRet(idxA);

            if tDepAst < arrivalAst,   continue; end
            if tArrEarth > finalMissionDate, continue; end

            totalMissionDays = tArrEarth - solutions(iSol).Departure;

            entry = struct( ...
                'AstName',          solutions(iSol).AstName, ...
                'AstID',            astID, ...
                'DepartEarth',      solutions(iSol).Departure, ...
                'ArriveAst',        arrivalAst, ...
                'DepartAst',        tDepAst, ...
                'ArriveEarth',      tArrEarth, ...
                'ArcTypeReturn',    "Long", ...
                'TOF_daysReturn',   tArrEarth - tDepAst, ...
                'vInfReturn',       vInfLongRet(idxD, idxA), ...
                'dvAstDep',         dvLongRet(idxD, idxA), ...
                'TotalMissionDays', totalMissionDays, ...
                ... % velocity vectors for the return
                'v1DepartVec',      v1LongCell{idxD, idxA}, ...
                'v2ArriveVec',      v2LongCell{idxD, idxA} ...
                );

            localSolutions = [localSolutions; entry];
        end

        % Store local results for this iSol
        solutionsCell{iSol} = localSolutions;
    end

    % Merge all sub-results into one struct array
    roundTripSolutions = vertcat(solutionsCell{:});

    % Show how many total round-trip solutions we found
    fprintf('\nFound %d feasible round-trip solutions (zero dwell)!\n', numel(roundTripSolutions));
end


function solutions = departureSolutions()
    % departureSolutions  Finds Earth->NEO transfer solutions in parallel,
    % storing velocity vectors from the updated findTransferSolutions.

    % 1. Paths & Data
    path = 'C:\Users\gadma\Documents\MATLAB\ATATD-Toolbox\All_NEOS_ATA&TD_2018_2019.csv';
    neo_data = ephdata.read_data(path);

    % 2. Define date window (MJD2000)
    date_dep = date2mjd2000([2033 01 01 0 0 0]);  % earliest Earth departure
    date_arr = date2mjd2000([2037 12 31 0 0 0]);  % must arrive asteroid before end-2037

    % Filter parameters
    smaRange = [0.7 1.3];  % AU
    eMax     = 0.4;
    iMax     = 10;         % deg

    % Load & filter the NEOs
    filteredNEOs = filterAsteroids(path, smaRange, eMax, iMax);

    % 3. Search Range & Step (days)
    range = date_arr - date_dep;  % e.g. 700 days if date_arr-date_dep=700
    step  = 20;   % 20-day increments

    % 4. Prepare a cell array to store each asteroid’s solutions
    total_asteroids = height(filteredNEOs);
    solutionsCell = cell(total_asteroids,1);

    % 5. Parallel loop over the filtered asteroids
    parfor iAst = 1:total_asteroids
        % Local array that will collect feasible solutions for this asteroid
        localSolutions = [];

        % Extract asteroid name
        thisName = filteredNEOs.full_name{iAst};  

        % Find asteroid ID in your data
        try
            astID = ephdata.getAsteroidId(neo_data, thisName);
        catch
            warning('Could not find ID for %s, skipping...', thisName);
            solutionsCell{iAst} = localSolutions; % Store empty result
            continue;
        end

        fprintf('\nExamining asteroid: %s\n Remaining: (%d/%d)\n', thisName, iAst, total_asteroids);

        % 6. Call the updated lambert-solver-based function that now
        %    also returns velocity vectors in cell arrays:
        [vInfShort, dvShort, vInfLong, dvLong, ...
         DepGrid, ArrGrid, ...
         v1ShortCell, v2ShortCell, v1LongCell, v2LongCell] = ...
            lambertsolver.findTransferSolutions( ...
                date_dep, date_arr, ...     % nominal center points
                3, astID, ...              % Earth->asteroid
                range, step, 'sun');

        % 7. Apply constraints: vInf < 1.5, dvRendez < 0.5
        shortMask = (vInfShort < 1.5) & (dvShort < 0.5);
        longMask  = (vInfLong  < 1.5) & (dvLong  < 0.5);

        foundAny = false;

        % --- Short Arc Feasible Solutions
        [iD_s, iA_s] = find(shortMask);
        for k = 1:numel(iD_s)
            iD = iD_s(k);
            iA = iA_s(k);

            % Build an entry struct
            entry = struct( ...
                'AstName',   thisName, ...
                'AstID',     astID, ...
                'Departure', DepGrid(iD), ...
                'Arrival',   ArrGrid(iA), ...
                'ArcType',   "Short", ...
                'vInf',      vInfShort(iD, iA), ...
                'dvRendez',  dvShort(iD, iA), ...
                'TOF_days',  (ArrGrid(iA) - DepGrid(iD)), ...
                ... % Store the actual velocity vectors in heliocentric frame
                'v1DepartVec', v1ShortCell{iD, iA}, ...
                'v2ArriveVec', v2ShortCell{iD, iA} ...
                );

            localSolutions = [localSolutions; entry];
            foundAny = true;
        end

        % --- Long Arc Feasible Solutions
        [iD_l, iA_l] = find(longMask);
        for k = 1:numel(iD_l)
            iD = iD_l(k);
            iA = iA_l(k);

            entry = struct( ...
                'AstName',   thisName, ...
                'AstID',     astID, ...
                'Departure', DepGrid(iD), ...
                'Arrival',   ArrGrid(iA), ...
                'ArcType',   "Long", ...
                'vInf',      vInfLong(iD, iA), ...
                'dvRendez',  dvLong(iD, iA), ...
                'TOF_days',  (ArrGrid(iA) - DepGrid(iD)), ...
                ... % velocity vectors
                'v1DepartVec', v1LongCell{iD, iA}, ...
                'v2ArriveVec', v2LongCell{iD, iA} ...
                );

            localSolutions = [localSolutions; entry];
            foundAny = true;
        end

        if foundAny
            fprintf('Found feasible Earth->Ast transfers for %s\n', thisName);
        else
            fprintf('No Earth->Ast transfer found for %s under constraints\n', thisName);
        end

        % Store results in the cell array
        solutionsCell{iAst} = localSolutions;
    end

    % 8. Merge all sub-results into one struct array
    solutions = vertcat(solutionsCell{:});

    % 9. Show final solutions or save them
    disp('---- SUMMARY OF ALL FEASIBLE EARTH->ASTEROID SOLUTIONS ----');
    if ~isempty(solutions)
        solutionsTable = struct2table(solutions);
        solutionsTable = sortrows(solutionsTable, 'Departure');
        disp(solutionsTable);
    else
        disp('No feasible solutions found with the given constraints!');
    end

    % 10. Save to a .mat for subsequent steps
    save('DepartureSolutions.mat','solutions');
end


function filteredTbl = filterAsteroids(csvPath, smaRange, eMax, iMax)
    % FILTERASTEROIDS filters asteroids by orbital parameters using ephdata methods.
    %
    %  INPUTS:
    %   csvPath  - Path to the CSV file with asteroid data.
    %   smaRange - 2-element vector [a_min, a_max], the allowable semi-major axis range (AU).
    %   eMax     - Maximum eccentricity limit.
    %   iMax     - Maximum inclination limit (degrees).
    %
    %  OUTPUT:
    %   filteredTbl - Table with only those asteroids that fit the specified orbital constraints.

    % 1) Read the CSV data into a table
    dataTbl = ephdata.read_data(csvPath);

    % 2) (Optional) Clean object names for clarity
    % [~, dataTbl] = ephdata.clean_object_names(dataTbl);

    % 3) Select only the columns we actually need
    %    (change column names to match exactly what's in your CSV)
    columnsNeeded = {'full_name','a','e','i'};
    subsetTbl = ephdata.select_columns(dataTbl, columnsNeeded);

    % 4) Build a logical mask to filter the table
    %    For example, a in [smaRange(1), smaRange(2)],
    %                   e in [0, eMax],
    %                   i in [0, iMax] (assuming i is in degrees)
    condition = (subsetTbl.a >= smaRange(1)) & (subsetTbl.a <= smaRange(2)) & ...
                (subsetTbl.e <= eMax) & ...
                (subsetTbl.i <= iMax);

    % 5) Use ephdata.filter_data to filter the subset
    filteredTbl = ephdata.filter_data(subsetTbl, condition);

    % Print summary
    fprintf('Filtered down to %d asteroids (from %d total).\n', ...
        height(filteredTbl), height(dataTbl));
end
