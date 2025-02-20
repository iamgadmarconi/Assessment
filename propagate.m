function [times, rSC] = propagate(r0, v0, t0, t1, muSun, nSteps)
    % PROPAGATETRAJECTORY
    %   Integrates a 2-body trajectory about the Sun from time t0 to t1.
    % 
    % Inputs:
    %   r0, v0 : 1x3 vectors, spacecraft position & velocity in heliocentric frame at t0
    %   t0, t1 : scalar times in MJD2000
    %   muSun  : gravitational parameter of the Sun (km^3/s^2)
    %   nSteps : how many points you want in the output
    %
    % Outputs:
    %   times : 1 x nSteps, the sample times (MJD2000)
    %   rSC   : nSteps x 3, spacecraft position at each step (heliocentric)
    
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
    
        % We'll sample uniformly in "nSteps" points
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
    