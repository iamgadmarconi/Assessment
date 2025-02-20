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

% After you merge, filter, and pick one entry:
mergedTable = mergeDepartureAndReturn(solutions, roundTrips);
filteredTable = filterByLayover(mergedTable, 60, 180);

% Select one row (e.g. the first)
entry = filteredTable(1, :);

% Generate the PDF
pdfFile = generatePDFReport(entry, filteredTable);
open(pdfFile);  % automatically opens in system PDF viewer (on most OS)


function pdfPath = generatePDFReport(entry, filteredTable)
    % generatePDFReport Generate a PDF report, compatible with older MATLAB versions
    %
    %  pdfPath = generatePDFReport(entry, filteredTable)
    %
    %  This function:
    %   1) Creates a title page: "(AsteroidName, dd/mm/yyyy)"
    %   2) Saves a trajectory figure to PNG and embeds it
    %   3) Creates a multi-page DOM table for 'filteredTable'
    %      using basic classes (Table, TableRow, TableEntry).

    % Make sure the MATLAB Report Generator is installed
    import mlreportgen.report.*    % For Report, TitlePage, Chapter, etc.
    import mlreportgen.dom.*       % For Table, TableRow, TableEntry, Image, etc.

    % --- Convert MJD2000 -> MATLAB datenum -> display strings ---
    MJD2000_OFFSET     = 730486;  % 1-Jan-2000 = datenum(2000,1,1) = 730486
    depEpochMatlab     = entry.EarthDepartureEpoch + MJD2000_OFFSET;
    departureDateSlash = datestr(depEpochMatlab, 'dd/mm/yyyy');   % displayed in PDF
    departureDateFile  = datestr(depEpochMatlab, 'dd-mm-yyyy');   % safe for filenames

    % Asteroid name is in a cell
    astName  = strtrim(entry.AstName{1});
    pdfTitle = 'Round-Trip Mission Report';

    % --- Construct a filename without illegal chars ---
    pdfFilename = sprintf('%s_%s', astName, departureDateFile);
    % Replace any non-word or special chars with underscore
    pdfFilename = regexprep(pdfFilename, '[^\w\-\(\) ]', '_');  

    % Create the report with the given filename, in PDF format
    rpt = Report(pdfFilename, 'pdf');
    rpt.Layout.Landscape = true;

    %% Title Page
    tp = TitlePage();
    tp.Title    = Text(pdfTitle);
    tp.Title.Style = {Bold(true), FontSize('24pt')};
    add(rpt, tp);
    %% Chapter for details
    ch = Chapter('Title', 'Mission Details');


    % ---------------------------
    % 2) Build a multi-page table
    % ---------------------------
    add(ch, Text('Filtered Round-Trip Table:'));

    % Create a DOM Table object
    myTable = Table();
    myTable.Border       = 'solid';
    myTable.BorderWidth  = '1pt';
    myTable.ColSep       = 'solid';
    myTable.RowSep       = 'solid';
    myTable.Width        = '100%';
    myTable.TableEntriesStyle = {FontSize('8pt')};  % Decrease font size for many columns

    % Allow splitting across pages
    % (If your older version doesn't have this property, comment out)
    % Add a header row
    % Get variable names from the filtered table
    varNames  = filteredTable.Properties.VariableNames;
    exclusions = {'V1DepartVecEarth', 'V2ArriveVecAsteroid', 'V1DepartVecAsteroid', 'V2ArriveVecEarth'};
    
    % Determine which indices to include
    includeIdx = find(~ismember(varNames, exclusions));
    
    % Build the header row using only the included variable names
    headerRow = TableRow();
    headerRow.Style = {Bold(true)};
    for idx = includeIdx
        entryObj = TableEntry(varNames{idx});
        append(headerRow, entryObj);
    end
    append(myTable, headerRow);
    
    % Convert the table data to cell array
    rawData = table2cell(filteredTable);
    [nRows, ~] = size(rawData);
    
    % Build each data row using only the columns we want
    for r = 1:nRows
        rowObj = TableRow();
        for idx = includeIdx
            valString = formatValue(rawData{r, idx});
            entryObj  = TableEntry(valString);
            append(rowObj, entryObj);
        end
        append(myTable, rowObj);
    end


    % Add the table to the chapter
    add(ch, myTable);

    % Add the chapter to the report
    add(rpt, ch);

    % Close report
    close(rpt);
    pdfPath = rpt.OutputPath;
    fprintf('Report generated: %s\n', pdfPath);
end


function [header, body] = createTableData(T)
    % Convert your MATLAB table to the header/body for a FormalTable
    varNames = T.Properties.VariableNames;
    numCols  = length(varNames);
    numRows  = height(T);

    % Header
    header = cell(1, numCols);
    for c = 1:numCols
        header{c} = Paragraph(varNames{c});
    end

    % Body
    body = cell(numRows, numCols);
    rawData = table2cell(T);
    for r = 1:numRows
        for c = 1:numCols
            valStr = formatValue(rawData{r,c});
            body{r,c} = Paragraph(valStr);
        end
    end
end

function strOut = formatValue(val)
    % Convert numeric/logical/etc. to string for the table
    if isnumeric(val)
        strOut = num2str(val);
    elseif islogical(val)
        strOut = string(val);
    elseif ischar(val)
        strOut = val;
    elseif isstring(val)
        strOut = char(val);
    elseif iscell(val)
        strOut = strjoin(string(val), ', ');
    else
        strOut = 'N/A';
    end
end


function filtered = filterByLayover(mergedTable, layover_min, layover_max)
    filtered = ephdata.filter_data(mergedTable, mergedTable.Layover > layover_min & mergedTable.Layover < layover_max);
    filtered = sortrows(filtered, 'AstID');
    fprintf('Filtered table has %d matching round-trip entries.\n', height(filtered));
end


function mergedTable = mergeDepartureAndReturn(departSol, roundTripSol)

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
