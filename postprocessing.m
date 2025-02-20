if ~exist('ephdata', 'class')
    error('ephdata class not found. Verify file location and MATLAB path.');
end

if ~exist('lambertsolver', 'class')
    error('lambertsolver class not found. Verify file location and MATLAB path.');
end

if ~exist('propagate.m', 'file')
    error('propagate.m not found. Verify file location and MATLAB path.');
end

if ~exist('roundTripSolutions.mat', 'file')
    error('roundTripSolutions.mat not found. Verify file location and MATLAB path.');
end

if ~exist('departureSolutions.mat', 'file')
    error('departureSolutions.mat not found. Verify file location and MATLAB path.');
end

load('DepartureSolutions.mat','solutions');     % Earth->Ast
load('RoundTripSolutions.mat','roundTrips');    % full round trip

mergedTable = mergeDepartureAndReturn(solutions, roundTrips);
filteredTable = filterByLayover(mergedTable, 60, 180);

entry = filteredTable(1,:);

plotRoundTrip(entry, 1000)



function filtered = filterByLayover(mergedTable, layover_min, layover_max)
    filtered = ephdata.filter_data(mergedTable, mergedTable.Layover > layover_min & mergedTable.Layover < layover_max);
end


function mergedTable = mergeDepartureAndReturn(departSol, roundTripSol)
% MERGEDEPARTUREANDRETURN
%   Merges Earth->Asteroid solutions ("departSol") with the matching 
%   Asteroid->Earth solutions ("roundTripSol") into a single table.
%
%   Two entries match if:
%       AstID is the same
%       Earth departure times match (rounded)
%       Asteroid arrival times match (rounded)
%
%   The result is a table that has both outbound and return columns.
%
%   If no matches are found, 'mergedTable' is empty.

    if isempty(departSol) || isempty(roundTripSol)
        warning('One of the solution sets is empty. Returning empty table.');
        mergedTable = table();
        return;
    end

    % 1) Convert the struct arrays to tables
    depTab = struct2table(departSol, 'AsArray', true);
    retTab = struct2table(roundTripSol, 'AsArray', true);

    % 2) Rename columns in depTab to clearly indicate "Outbound" data
    %    (so they don't clash with the return columns)
    depTab = renamevars(depTab, ...
        {'Departure','Arrival','ArcType','vInf','dvRendez','TOF_days','v1DepartVec','v2ArriveVec'}, ...
        {'EarthDepartureEpoch','AsteroidArrivalEpoch','DepartureArcType','DepartureVInf','DepartureDV','DepartureTOF','V1DepartVecEarth','V2ArriveVecAsteroid'});

    % 3) Rename columns in retTab to clearly indicate "Return" data
    %    (some fields in roundTripSolutions have similar names)
    retTab = renamevars(retTab, ...
        {'DepartEarth','ArriveAst','DepartAst','ArriveEarth','ArcTypeReturn','vInfReturn','dvAstDep','TOF_daysReturn','v1DepartVec','v2ArriveVec'}, ...
        {'EarthDepartureEpoch','AsteroidArrivalEpoch','AsteroidDepartureEpoch','EarthArrivalEpoch','ReturnArcType','ReturnVInf','ReturnDV','ReturnTOF','V1DepartVecAsteroid','V2ArriveVecEarth'});

    % 4) Create "key" columns in both tables
    %    We'll match by (AstID, OutDeparture=RetDepartEarth, OutArrival=RetArriveAst),
    %    rounding to integer days so they match exactly in join
    depTab.KeyAstID = depTab.AstID;  % same ID field, just clarifying the usage
    depTab.DepEpoch = round(depTab.EarthDepartureEpoch);
    depTab.ArrEpoch = round(depTab.AsteroidArrivalEpoch);

    retTab.KeyAstID = retTab.AstID;
    retTab.DepEpoch = round(retTab.EarthDepartureEpoch);
    retTab.ArrEpoch = round(retTab.AsteroidArrivalEpoch);

    % 5) Perform an inner join using these keys
    %    'MergeKeys' => merges them into single columns
    mergedTable = innerjoin(depTab, retTab, ...
        'LeftKeys',  {'KeyAstID','DepEpoch','ArrEpoch'}, ...
        'RightKeys', {'KeyAstID','DepEpoch','ArrEpoch'}, ...
        'LeftVariables', {
            'AstName','AstID','EarthDepartureEpoch','AsteroidArrivalEpoch','DepartureArcType','DepartureVInf',...
            'DepartureDV','DepartureTOF','V1DepartVecEarth','V2ArriveVecAsteroid'
        }, ...
        'RightVariables', {
            'AsteroidDepartureEpoch','EarthArrivalEpoch','ReturnArcType','ReturnVInf','ReturnDV','ReturnTOF','V1DepartVecAsteroid','V2ArriveVecEarth'
        });


    % The resulting table will have all columns from depTab and retTab
    % for rows that matched on the 3 key columns.

    % 6) Sort by outbound departure for convenience
    if ~isempty(mergedTable)
        mergedTable.Layover = mergedTable.AsteroidDepartureEpoch - mergedTable.AsteroidArrivalEpoch;
        mergedTable = sortrows(mergedTable, 'Layover');
    end

    fprintf('Merged table has %d matching round-trip entries.\n', height(mergedTable));
end


function plotRoundTrip(oneRow, nSteps)
    % PLOTROUNDTRIP
    %   Plots Earth->Asteroid then Asteroid->Earth in two separate figures,
    %   adding markers for departure/arrival points.

    figure; hold on;

    muSun = getAstroConstants('Sun','mu');

    %% === Outbound Leg (Earth->Asteroid) ===
    t0_out = oneRow.EarthDepartureEpoch;
    t1_out = oneRow.AsteroidArrivalEpoch;

    % Earth at departure
    [rE0, ~] = EphSS_car(3, t0_out);
    vSC0_out = oneRow.V1DepartVecEarth;  % Lambert velocity in heliocentric

    [timeOut, rSC_out] = propagate(rE0, vSC0_out, t0_out, t1_out, muSun, nSteps);

    % Earth & asteroid positions for each sample time
    rEarth_out = zeros(nSteps,3);
    rAst_out   = zeros(nSteps,3);

    for i = 1:nSteps
        [rE, ~] = EphSS_car(3, timeOut(i));
        [rA, ~] = EphSS_car(oneRow.AstID, timeOut(i));
        rEarth_out(i,:) = rE;
        rAst_out(i,:)   = rA;
    end

    % (1) Plot the main trajectories
    plot3(rSC_out(:,1), rSC_out(:,2), rSC_out(:,3), '--k','LineWidth',1.5);
    plot3(rEarth_out(:,1), rEarth_out(:,2), rEarth_out(:,3), 'b','LineWidth',1);
    plot3(rAst_out(:,1),   rAst_out(:,2),   rAst_out(:,3),   'r','LineWidth',1);

    % (2) Add markers for departure & arrival on spacecraft track
    % The spacecraft's first point is departure, last point is arrival
    % plot3(rSC_out(1,1), rSC_out(1,2), rSC_out(1,3), '^k','MarkerSize',8,'MarkerFaceColor','k',...
    %       'DisplayName','SC Depart');
    % plot3(rSC_out(end,1), rSC_out(end,2), rSC_out(end,3), 'vk','MarkerSize',8,'MarkerFaceColor','k',...
    %       'DisplayName','SC Arrive');

    % Earth departure (circle, outline + fill both blue):
    plot3(rEarth_out(1,1), rEarth_out(1,2), rEarth_out(1,3), ...
        'ob','MarkerSize',6,'MarkerFaceColor','b',...
        'DisplayName','Earth@Depart');

    % Ast departure (circle, blue outline, red fill):
    plot3(rAst_out(1,1), rAst_out(1,2), rAst_out(1,3), ...
        'sr','MarkerSize',6,'MarkerFaceColor','r',...
        'DisplayName','Ast@Depart');

    % Earth arrival (square, red outline, blue fill):
    plot3(rEarth_out(end,1), rEarth_out(end,2), rEarth_out(end,3), ...
        '^','MarkerSize',6,'MarkerFaceColor','b',...
        'DisplayName','Earth@Arrive');

    % Ast arrival (square, red outline + fill):
    plot3(rAst_out(end,1), rAst_out(end,2), rAst_out(end,3), ...
        'v','MarkerSize',6,'MarkerFaceColor','r',...
        'DisplayName','Ast@Arrive');


    % (4) Plot the Sun
    plot3(0,0,0,'yo','MarkerFaceColor','y','MarkerSize',8);

    axis equal; grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    missionTitle = ['Departure Leg: ' strtrim(oneRow.AstName{1}) ' (Earth->Ast)'];
    title(missionTitle);

    legend('Outbound SC','Earth','Asteroid','Earth@Depart','Ast@Depart','Earth@Arrive','Ast@Arrive','Sun','Location','best');

    %% === Return Leg (Asteroid->Earth) ===
    figure; hold on;

    t0_ret = oneRow.AsteroidDepartureEpoch;
    t1_ret = oneRow.EarthArrivalEpoch;

    [rAst0, ~] = EphSS_car(oneRow.AstID, t0_ret);
    vSC0_ret   = oneRow.V1DepartVecAsteroid;

    [timeRet, rSC_ret] = propagate(rAst0, vSC0_ret, t0_ret, t1_ret, muSun, nSteps);

    rEarth_ret = zeros(nSteps,3);
    rAst_ret   = zeros(nSteps,3);

    for i = 1:nSteps
        [rE, ~] = EphSS_car(3, timeRet(i));
        [rA, ~] = EphSS_car(oneRow.AstID, timeRet(i));
        rEarth_ret(i,:) = rE;
        rAst_ret(i,:)   = rA;
    end

    plot3(rSC_ret(:,1), rSC_ret(:,2), rSC_ret(:,3), '--k','LineWidth',1.5);
    plot3(rEarth_ret(:,1), rEarth_ret(:,2), rEarth_ret(:,3), 'b','LineWidth',1);
    plot3(rAst_ret(:,1),   rAst_ret(:,2),   rAst_ret(:,3),   'r','LineWidth',1);

    % Mark departure & arrival on the return leg
    % plot3(rSC_ret(1,1), rSC_ret(1,2), rSC_ret(1,3), '^k','MarkerSize',8,'MarkerFaceColor','k',...
    %       'DisplayName','SC Depart');
    % plot3(rSC_ret(end,1), rSC_ret(end,2), rSC_ret(end,3), 'vk','MarkerSize',8,'MarkerFaceColor','k',...
    %       'DisplayName','SC Arrive');

    % Mark asteroid at departure, Earth at arrival
    plot3(rAst_ret(1,1), rAst_ret(1,2), rAst_ret(1,3), 'sr','MarkerSize',6,'MarkerFaceColor','r',...
        'DisplayName','Ast@Depart');
    plot3(rEarth_ret(1,1), rEarth_ret(1,2), rEarth_ret(1,3), 'ob','MarkerSize',6,'MarkerFaceColor','b',...
        'DisplayName','Earth@Depart');

    plot3(rAst_ret(end,1), rAst_ret(end,2), rAst_ret(end,3), 'v','MarkerSize',6,'MarkerFaceColor','r',...
        'DisplayName','Ast@Arrive');

    plot3(rEarth_ret(end,1), rEarth_ret(end,2), rEarth_ret(end,3), '^','MarkerSize',6,'MarkerFaceColor','b',...
        'DisplayName','Earth@Arrive');

    % Sun
    plot3(0,0,0,'yo','MarkerFaceColor','y','MarkerSize',8);

    axis equal; grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    missionTitle = ['Return Leg:' strtrim(oneRow.AstName{1}) ' (Ast->Earth)'];
    title(missionTitle);

    legend('Return SC','Earth','Asteroid','Ast@Depart','Earth@Depart','Ast@Arrive','Earth@Arrive','Sun','Location','best');
end

