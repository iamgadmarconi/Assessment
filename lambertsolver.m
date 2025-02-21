classdef lambertsolver

    methods (Static)
        
        
        function plotPorkChop(date_dep, date_arr, dep_body_idx, arr_body_idx, range, step, central_body)

            [DV_solutions, DepartureGrid, ArrivalGrid] = lambertsolver.findDVSolutions( ...
                date_dep, date_arr, dep_body_idx, arr_body_idx, range, step, central_body);
            
            minDV = min(min(DV_solutions)); 
            figure;
            contourf(DepartureGrid, ArrivalGrid, DV_solutions', minDV:0.1:8.5);
            colorbar;
            xlabel('Departure Date (MJD2000)');
            ylabel('Arrival Date (MJD2000)');
            title('\DeltaV Pork-Chop Plot');
        end


        function [vInfShortMat, dvRendShortMat, ...
            vInfLongMat,  dvRendLongMat, ...
            DepartureGrid, ArrivalGrid, ...
            v1ShortCell,  v2ShortCell, ...
            v1LongCell,   v2LongCell] = ...
            findTransferSolutions(date_dep, date_arr, dep_body_idx, arr_body_idx, ...
                                    range, step, central_body)

            DepartureGrid = (date_dep - range) : step : (date_dep + range);
            ArrivalGrid   = (date_arr - range) : step : (date_arr + range);
        
            nDep = numel(DepartureGrid);
            nArr = numel(ArrivalGrid);
        
            % Initialize output arrays for vInf, dv
            vInfShortMat   = NaN(nDep,nArr);
            dvRendShortMat = NaN(nDep,nArr);
            vInfLongMat    = NaN(nDep,nArr);
            dvRendLongMat  = NaN(nDep,nArr);
        
            % Initialize cell arrays for the actual velocity vectors
            v1ShortCell = cell(nDep, nArr);  % departure velocity (short arc)
            v2ShortCell = cell(nDep, nArr);  % arrival velocity   (short arc)
            v1LongCell  = cell(nDep, nArr);  % departure velocity (long arc)
            v2LongCell  = cell(nDep, nArr);  % arrival velocity   (long arc)
        
            % Precompute ephemerides for departure & arrival bodies
            rDepList = zeros(nDep,3);  vDepList = zeros(nDep,3);
            rArrList = zeros(nArr,3);  vArrList = zeros(nArr,3);
        
            for i = 1:nDep
                [rDepList(i,:), vDepList(i,:)] = EphSS_car(dep_body_idx, DepartureGrid(i));
            end
            for j = 1:nArr
                [rArrList(j,:), vArrList(j,:)] = EphSS_car(arr_body_idx, ArrivalGrid(j));
            end
        
            mu_central = getAstroConstants(central_body, 'mu');
        
            % Main parallel loop
            parfor iD = 1:nDep
                rDep = rDepList(iD,:);
                vDep = vDepList(iD,:);
                tDep = DepartureGrid(iD);
        
                for iA = 1:nArr
                    tArr = ArrivalGrid(iA);
                    TOF  = (tArr - tDep) * 86400;  % [s]
        
                    if TOF <= 0
                        % Negative or zero flight times are meaningless
                        continue;
                    end
        
                    rArr = rArrList(iA,:);
                    vArr = vArrList(iA,:);
        
                    % Lambert short arc
                    [v1_s, v2_s] = lambertsolver.lambertArcATD(rDep, rArr, TOF, +1, mu_central, 1e-3, 20);
                    if ~any(isnan(v1_s)) && ~any(isnan(v2_s))
                        % Record the heliocentric velocity vectors
                        v1ShortCell{iD, iA} = v1_s;
                        v2ShortCell{iD, iA} = v2_s;
        
                        % Also store the magnitudes of vInf, dv
                        vInfShortMat(iD,iA)   = norm(v1_s - vDep);  % departure v∞ wrt planet
                        dvRendShortMat(iD,iA) = norm(v2_s - vArr);  % arrival dv wrt target
                    end
        
                    % Lambert long arc
                    [v1_l, v2_l] = lambertsolver.lambertArcATD(rDep, rArr, TOF, -1, mu_central, 1e-3, 20);
                    if ~any(isnan(v1_l)) && ~any(isnan(v2_l))
                        v1LongCell{iD, iA} = v1_l;
                        v2LongCell{iD, iA} = v2_l;
        
                        vInfLongMat(iD,iA)   = norm(v1_l - vDep);
                        dvRendLongMat(iD,iA) = norm(v2_l - vArr);
                    end
                end
            end
        end

        function [DV_solutions, DepartureGrid, ArrivalGrid] = findDVSolutions( ...
            date_dep, date_arr, dep_body_idx, arr_body_idx, range, step, central_body)

            DepartureGrid = (date_dep - range) : step : (date_dep + range);
            ArrivalGrid   = (date_arr - range) : step : (date_arr + range);

            nDep = numel(DepartureGrid);
            nArr = numel(ArrivalGrid);

            DV_solutions = NaN(nDep, nArr);

            % Precompute ephemerides
            rDepList = zeros(nDep,3);  vDepList = zeros(nDep,3);
            rArrList = zeros(nArr,3);  vArrList = zeros(nArr,3);

            for i=1:nDep
                [rDepList(i,:), vDepList(i,:)] = EphSS_car(dep_body_idx, DepartureGrid(i));
            end
            for j=1:nArr
                [rArrList(j,:), vArrList(j,:)] = EphSS_car(arr_body_idx, ArrivalGrid(j));
            end

            mu_central = getAstroConstants(central_body, 'mu');

            % parallelize outer loop
            parfor iD = 1:nDep
                rDep = rDepList(iD,:);
                vDep = vDepList(iD,:);
                tDep = DepartureGrid(iD);

                for iA = 1:nArr
                    tArr = ArrivalGrid(iA);
                    TOF  = (tArr - tDep)*86400; % days to seconds

                    if TOF <= 0
                        continue;
                    end

                    rArr = rArrList(iA,:);
                    vArr = vArrList(iA,:);

                    [v1_short, v2_short] = lambertsolver.lambertArcATD( ...
                                        rDep, rArr, TOF, +1, mu_central, 1e-3, 20);
                    [v1_long,  v2_long ] = lambertsolver.lambertArcATD( ...
                                        rDep, rArr, TOF, -1, mu_central, 1e-3, 20);

                    if any(isnan(v1_short)) || any(isnan(v2_short))
                        DV_short = Inf;
                    else
                        DV_short = norm(v1_short - vDep) + norm(v2_short - vArr);
                    end

                    if any(isnan(v1_long)) || any(isnan(v2_long))
                        DV_long = Inf;
                    else
                        DV_long  = norm(v1_long  - vDep) + norm(v2_long  - vArr);
                    end

                    DV_solutions(iD,iA) = min(DV_short, DV_long);
                end
            end
        end


        function vrow = rowvec(v)
            % Ensures the output is a 1x3 row vector
            vrow = v(:).';  % Flatten to column, then transpose -> 1xN
        end


        function Phi = STMLambert(r, v, dt, mu)
            % Compute the 3x3 STM for small changes 
            %   in initial velocity for Lambert arcs (row-vector version).

            r = lambertsolver.rowvec(r);
            v = lambertsolver.rowvec(v);

            % call fGKepler
            [rf, vf_calc, dE] = lambertsolver.fGKepler(r, v, dt, mu, 1e-8);

            if any(isnan(rf)) || any(isnan(vf_calc))
                Phi = NaN(3,3);
                return;
            end

            a = mu / (2*mu/norm(r) - norm(v)^2);

            F = 1 - a/norm(r)*(1 - cos(dE));
            G = dt + sqrt(a^3/mu)*( sin(dE) - dE );

            dF = - sqrt(mu*a)/(norm(rf)*norm(r)) * sin(dE);
            dG = 1 - a/norm(rf)*(1 - cos(dE));

            C = a*sqrt(a^3/mu)*(3*sin(dE) - (2+cos(dE))*dE) ...
                - a*dt*(1-cos(dE));

            % compute dr, dv
            vf_calc = dF*r + dG*v;     % row vectors
            dr      = rf - r;         % 1x3
            dv      = vf_calc - v;    % 1x3

            term1 = (norm(r)/mu)*(1 - F)*((dr')*v - (dv')*r);
            term2 = (C/mu)*((vf_calc')*v);
            term3 = G * eye(3);

            Phi = term1 + term2 + term3;
        end


        function [rf, vf, dE] = fGKepler(r, v, dt, mu, tol)
            % Universal f & g function to propagate from (r,v).

            r = lambertsolver.rowvec(r);
            v = lambertsolver.rowvec(v);

            a = mu / (2*mu/norm(r) - norm(v)^2);
            n = sqrt(mu / a^3);   % mean motion
            dM = n*dt;

            sigma_0 = dot(r,v)/sqrt(mu);
            dE = dM;  % initial guess
            fE = dE - (1 - norm(r)/a)*sin(dE) - sigma_0/sqrt(a)*(cos(dE)-1) - dM;

            i_iter = 0;
            i_iter_max = 100;

            while (abs(fE) > tol) && (i_iter < i_iter_max)
                dfE = 1 - (1 - norm(r)/a)*cos(dE) + sigma_0/sqrt(a)*sin(dE);
                dE  = dE - fE/dfE;

                fE  = dE - (1 - norm(r)/a)*sin(dE) ...
                        - sigma_0/sqrt(a)*(cos(dE)-1) - dM;

                i_iter = i_iter + 1;
            end

            if abs(fE) > tol
                % warning(['fGKepler did not converge. fE=', num2str(fE)]);
                rf = [NaN, NaN, NaN];
                vf = [NaN, NaN, NaN];
                return;
            end

            F = 1 - a/norm(r)*(1 - cos(dE));
            G = dt + sqrt(a^3/mu)*(sin(dE) - dE);

            rf = F*r + G*v;  % 1x3

            if norm(rf) < 1e-12
                % warning('rf is extremely small');
                rf = [NaN, NaN, NaN];
                vf = [NaN, NaN, NaN];
                return;
            end

            dG = 1 - a/norm(rf)*(1 - cos(dE));
            vf = (1/G)*( dG*rf - r );

            rf = lambertsolver.rowvec(rf);
            vf = lambertsolver.rowvec(vf);
        end


        function [a_min, e_min, dt_min, v] = minETransfer(r1, r2, t_m, mu)
            % Quick guess for minimal-energy transfer.

            c = norm(r2 - r1);
            r1_norm = norm(r1);
            r2_norm = norm(r2);

            a_min = 1/4 * (r1_norm + r2_norm + c);

            cos_dTA = dot(r1, r2) / (r1_norm*r2_norm);
            sin_dTA = t_m * sqrt(1 - cos_dTA^2);

            p_min = (r1_norm * r2_norm) / c * (1 - cos_dTA);
            e_min = sqrt(1 - p_min/a_min);

            beta_e = 2 * asin(sqrt((2*a_min - c)/(2*a_min)));
            dt_min = sqrt(a_min^3/mu) * (pi - t_m*(beta_e - sin(beta_e)));

            F = 1 - r2_norm/p_min*(1 - cos_dTA);
            G = r2_norm*r1_norm / sqrt(mu*p_min)*sin_dTA;

            v = 1/G*(r2 - F*r1);
        end


        function [v1, v2] = lambertArcATD(r1, r2, dt_target, tm, mu, tol, max_iter)

            % Get an initial guess from minE solution:
            [~, ~, ~, v1] = lambertsolver.minETransfer(r1, r2, tm, mu);
        
            % Newton iteration at dt_target
            for iter = 1:max_iter
                [r2_s, v2_s, ~] = lambertsolver.fGKepler(r1, v1, dt_target, mu, 1e-8);
                if any(isnan(r2_s)),  v1 = NaN(1,3); v2 = NaN(1,3); return; end
        
                epsilon = r2 - r2_s;
                err_val = norm(epsilon);
        
                if err_val < tol
                    v2 = v2_s;
                    return  % converged
                end
        
                Phi = lambertsolver.STMLambert(r1, v1, dt_target, mu);
                if any(any(isnan(Phi))), v1 = NaN(1,3); v2 = NaN(1,3); return; end
        
                delta_v = Phi \ epsilon.';   % Solve for velocity correction
                v1 = v1 + delta_v.';        % Update
            end
        
            % If we exit loop, not converged
            v1 = NaN(1,3);
            v2 = NaN(1,3);
        end        
    end % methods
end
