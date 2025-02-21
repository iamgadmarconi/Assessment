function [times, rSC] = propagate(r0, v0, t0, t1, muSun, nSteps)
        % Convert days to seconds for integration
        t0_sec = t0 * 86400;
        t1_sec = t1 * 86400;
        Tspan  = t1_sec - t0_sec;
    
        % ODE initial state
        y0 = [r0(:); v0(:)];
    
        % ODE45 parameters
        odeOpts = odeset('RelTol',1e-9,'AbsTol',1e-12);
    
        % 2-body function
        twoBodySun = @(t, y) [y(4:6);
                              -muSun * y(1:3)/norm(y(1:3))^3 ];
    
        % Integrate from 0 to Tspan
        [tout,yout] = ode45(twoBodySun, [0, Tspan], y0, odeOpts);
    
        timesSec = linspace(0, Tspan, nSteps);
    
        rSC = zeros(nSteps,3);
    
        for i = 1:nSteps
            % Interpolate the ODE solution at time = timesSec(i)
            yinterp = interp1(tout, yout, timesSec(i));
            rSC(i,:) = yinterp(1:3);
        end
    
        % Convert back to MJD2000
        times = linspace(t0, t1, nSteps);
    end
    