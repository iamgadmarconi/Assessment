clc; clear; close all;

if ~exist('lambertsolver', 'class')
    error('lambertsolver class not found. Verify file location and MATLAB path.');
end

load("filtered.mat");


plotDepAndRetPorkchop(filtered, 1);


function plotDepAndRetPorkchop(table, entry_idx)
    entry = table(entry_idx,:);
    range_departure = max(table.DepartureTOF);
    range_return = max(table.ReturnTOF);

    lambertsolver.plotPorkChop(entry.EarthDepartureEpoch, entry.AsteroidArrivalEpoch, 3, entry.AstID, range_departure, 1, 'Sun')
    lambertsolver.plotPorkChop(entry.AsteroidDepartureEpoch, entry.EarthArrivalEpoch, entry.AstID, 3, range_return, 1, 'Sun')

    plotRoundTrip(entry, 1000);

end


function plotRoundTrip(oneRow, nSteps)
    % Create one figure with two subplots (top: Departure Leg, bottom: Return Leg)
    figure;
    
    %% === Departure Leg (Earth->Asteroid) ===
    subplot(2,1,1); hold on;
    
    muSun = getAstroConstants('Sun','mu');
    
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
    
    % Plot the main trajectories
    plot3(rSC_out(:,1), rSC_out(:,2), rSC_out(:,3), '--k','LineWidth',1.5);
    plot3(rEarth_out(:,1), rEarth_out(:,2), rEarth_out(:,3), 'b','LineWidth',1);
    plot3(rAst_out(:,1),   rAst_out(:,2),   rAst_out(:,3),   'r','LineWidth',1);
    
    % Markers for departure & arrival
    plot3(rEarth_out(1,1), rEarth_out(1,2), rEarth_out(1,3), 'ob','MarkerSize',6,...
        'MarkerFaceColor','b','DisplayName','Earth@Depart');
    plot3(rAst_out(1,1), rAst_out(1,2), rAst_out(1,3), 'sr','MarkerSize',6,...
        'MarkerFaceColor','r','DisplayName','Ast@Depart');
    plot3(rEarth_out(end,1), rEarth_out(end,2), rEarth_out(end,3), '^','MarkerSize',6,...
        'MarkerFaceColor','b','DisplayName','Earth@Arrive');
    plot3(rAst_out(end,1), rAst_out(end,2), rAst_out(end,3), 'v','MarkerSize',6,...
        'MarkerFaceColor','r','DisplayName','Ast@Arrive');
    
    % Plot the Sun
    plot3(0,0,0,'yo','MarkerFaceColor','y','MarkerSize',8);
    
    axis equal; grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    missionTitle = ['Departure Leg: ' strtrim(oneRow.AstName{1}) ' (Earth->Ast)'];
    title(missionTitle);
    legend('Outbound SC','Earth','Asteroid','Earth@Depart',...
        'Ast@Depart','Earth@Arrive','Ast@Arrive','Sun','Location','best');
    
    
    %% === Return Leg (Asteroid->Earth) ===
    subplot(2,1,2); hold on;
    
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
    
    % Mark asteroid at departure and Earth at arrival
    plot3(rAst_ret(1,1), rAst_ret(1,2), rAst_ret(1,3), 'sr','MarkerSize',6,...
        'MarkerFaceColor','r','DisplayName','Ast@Depart');
    plot3(rEarth_ret(1,1), rEarth_ret(1,2), rEarth_ret(1,3), 'ob','MarkerSize',6,...
        'MarkerFaceColor','b','DisplayName','Earth@Depart');
    plot3(rAst_ret(end,1), rAst_ret(end,2), rAst_ret(end,3), 'v','MarkerSize',6,...
        'MarkerFaceColor','r','DisplayName','Ast@Arrive');
    plot3(rEarth_ret(end,1), rEarth_ret(end,2), rEarth_ret(end,3), '^','MarkerSize',6,...
        'MarkerFaceColor','b','DisplayName','Earth@Arrive');
    
    % Plot the Sun
    plot3(0,0,0,'yo','MarkerFaceColor','y','MarkerSize',8);
    
    axis equal; grid on;
    xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
    missionTitle = ['Return Leg: ' strtrim(oneRow.AstName{1}) ' (Ast->Earth)'];
    title(missionTitle);
    legend('Return SC','Earth','Asteroid','Ast@Depart',...
        'Earth@Depart','Ast@Arrive','Earth@Arrive','Sun','Location','best');
end
