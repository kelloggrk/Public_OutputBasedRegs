% The class that instantiates the indexmodel object 
classdef indexmodel
    properties
        % We specify here all of the properties we want to have
        % attached to each object
        BQQ
        BQE
        BEE
        BEF
        BQeta
        phi
        phiact
        gamma0
        Q0
        E0
        sigma
        mu00
        Nz
        Z
        SOC
        EE
        EQ
    end
    methods
        % Constructor
        function obj = indexmodel(params)
            % Define parameters that vary by application
            obj.BQQ = params.BQQ; obj.BQE = params.BQE; obj.BEE = params.BEE; 
            obj.BEF = params.BEF; obj.BQeta = params.BQeta; obj.phi = params.phi; 
            obj.phiact = params.phiact; obj.gamma0 = params.gamma0;
            obj.Q0 = params.Q0; obj.E0 = params.E0; 
            obj.sigma = params.sigma;
            
            % Impute actual regulatory intercept mu00
            obj.mu00 = obj.E0 - obj.gamma0 * obj.Q0;            
            
            % Set up normal distribution for eta: mean, std dev, and
            % nodes for integration (assume mean of eta is zero (normalization))
            obj.Nz = 1000;                  % # of nodes
            Z = linspace(0,1,obj.Nz+2);     % vector of Zs to integrate over
            obj.Z = norminv(Z',0,obj.sigma);
            obj.Z = obj.Z(2:length(Z)-1);   % drop first and last extremes (zero and infty)            
 
            % Unconstrained Q and E at Eeta
            % Use fact that at (Q0,E0,Eeta), B_E = phiact
            denom = obj.BEE + obj.BQE * (-obj.BEE*obj.gamma0-obj.BQE)/(obj.BQE*obj.gamma0+obj.BQQ);
            obj.EE = obj.E0 - obj.phiact / denom;   % unconstrained value of E at Eeta
            % Now unconrstained value of Q at Eeta
            obj.EQ = obj.Q0 + (obj.EE-obj.E0) *...
                (-obj.BEE*obj.gamma0-obj.BQE)/(obj.BQE*obj.gamma0+obj.BQQ);            
            
            % Verify that second order condition holds
            obj.SOC = obj.BQQ * obj.BEE - obj.BQE^2;
            if obj.SOC<=0
                fprintf('Second order condition violated: B_QQ*B_EE-B_QE^2 <= 0 \n')
            else
            end
        end
        
        
        
        % Function for unconstrained line of Q and E choices
        function [Q,E] = uncon(obj,eta)
            % Finds unconstrained Q and E as a function of eta
            Q = obj.EQ + -obj.BEE*obj.BQeta / obj.SOC * eta;
            E = obj.EE +  obj.BQE*obj.BQeta / obj.SOC * eta;
        end
            
        
        
        % Total welfare function - not normalized
        function W = Bnn(obj,Q,E,eta)
            % INPUTS 
            % Q: quantity or attribute value
            % E: energy use (or efficiency) value
            % eta: value of shock
            
            % Use facts that on the actual reg constraint, B_Q = -gamma0 *
            % B_E, and that at (Q0,E0,Eeta), B_E = phiact
            
            % Define welfare using Taylor expansion around (Q0, E0, Eeta)
            % First-order terms use shadow values of actual regulatory constraint
            Wp = obj.phiact*(E-obj.E0) - obj.gamma0*obj.phiact*(Q-obj.Q0)...
                + 1/2*obj.BEE*(E-obj.E0).^2 + 1/2*obj.BQQ*(Q-obj.Q0).^2....
                + obj.BQE*(Q-obj.Q0).*(E-obj.E0) + obj.BQeta*(Q-obj.Q0).*eta;
            W = Wp - obj.phi*E;       % total welfare, including externality
        end
        
        
        
        % Total welfare function, normalized to equal zero at all points on
        % unconstrained choice line
        function W = B(obj,Q,E,eta)
            % INPUTS 
            % Q: quantity or attribute value
            % E: energy use (or efficiency) value
            % eta: value of shock
            
            % Normalization: total welfare is equal to zero at
            % (Q0,E0,Eeta). Also equal to zero anywhere else
            % along the unconstrained choice line, given a value of eta
            [Qu,Eu] = uncon(obj,eta);           % unconstrained Q and E
            Wu = Bnn(obj,Qu,Eu,eta);            % welfare when unconstrained
            W = Bnn(obj,Q,E,eta) - Wu;          % normalized welfare
        end
        
        
        
        % Socially optimal Q, E, and welfare, given eta
        function [W,Q,E] = socialopt(obj,eta)
            % INPUTS 
            % Q: quantity or attribute value
            % E: energy use (or efficiency) value
            % eta: value of shock
            
            % Get unconstrained location and then shift by dQ/dphi and dE/dphi
            [Qu,Eu] = uncon(obj,eta);
            Q = Qu - obj.BQE / obj.SOC * obj.phi;
            E = Eu + obj.BQQ / obj.SOC * obj.phi;
            
            % Welfare
            W = B(obj,Q,E,eta);
        end
        
        
        
        % Determine when regulation binds
        function [flag,hat] = etahat(obj,mu0,gamma)
            % INPUTS 
            % mu0: regulation intercept
            % gamma: regulation slope
            
            % OUTPUTS:
            % flag: 0 if never binds, 1 if always binds, 2 if binds for
            % eta>=hat, 3 if binds for eta<=hat
            % hat: value of eta at which regulation just binds
            
            % Compare gamma to slope of unconstrained line (= -B_QE / B_EE)
            slope = -obj.BQE / obj.BEE;
            diff = gamma - slope;
            
            [Qu,Eu] = uncon(obj,obj.Z);           % unconstrained choices
            Ec = mu0 + gamma * Qu;                % constrained choices
            
            if abs(diff) < 1E-10        % slopes are same
                hat = inf;              % output for etahat undefined
                % Compare location to see if reg always or never binds
                d = min(Eu-Ec);                     % positive if binds
                if d>0
                    flag = 1;                       % reg always binds
                else
                    flag = 0;                       % reg never binds
                end
            else                        % reg line and unconstrained lines intersect
                if diff < 0
                    flag = 2;           % reg binds for high eta
                else
                    flag = 3;           % reg binds for low eta
                end
                % Find Qhat at which lines intersect
                Qhat = Qu(1) + (Ec(1) - Eu(1)) / (-diff);
                % Find eta for that Qhat
                hat = obj.Z(1) + (Qhat - Qu(1)) * obj.SOC / (-obj.BEE * obj.BQeta);
            end
        end
        
        
        
        % Calculate Q, E, and welfare given regulation and eta
        function [W,Q,E] = regWQE(obj,mu0,gamma,eta)
            % INPUTS 
            % mu0: regulation intercept
            % gamma: regulation slope            
            % eta: value of shock
            
            % Find where regulation binds
            [flag, hat] = etahat(obj,mu0,gamma);
            
            % Get unconstrained locations
            [Qu,Eu] = uncon(obj,eta);
            
            % Get constrained E's at same Q's as unconstrained
            Ec = mu0 + gamma * Qu;
            
            % Loop over values in eta. Whenever reg binds, shift Q and E to
            % regulation using dQ/dmu0 and dE/dmu0 derivatives
            N= length(eta);
            Q = Qu; E = Eu;             % initialize Q and E with unconstrained values
            if flag>0
                if flag==1
                    ind = 1:N;          % change all Q,E
                elseif flag==2
                    ind = find(eta>hat);   % change Q,E for high eta
                else
                    ind = find(eta<hat);   % change Q,E for low eta
                end
                % Loop over values that need to be changed
                n = length(ind);
                denom = obj.BEE*gamma^2 + 2*obj.BQE*gamma + obj.BQQ;  % denom of dQ/dmu0 and dE/dmu0
                dQdmu0 = -(obj.BEE*gamma+obj.BQE) / denom;
                dEdmu0 = (obj.BQE*gamma+obj.BQQ) / denom;
                for i = 1:n
                    delta = Ec(ind(i)) - Eu(ind(i));        % vertical change in u0
                    % Shift Q and E per a shift in mu0 equal to delta
                    Q(ind(i)) = Qu(ind(i)) + dQdmu0 * delta;
                    E(ind(i)) = Eu(ind(i)) + dEdmu0 * delta;
                end
            else        % always unconstrained; do nothing
            end
            
            % Now calculate welfare
            W = B(obj,Q,E,eta);
        end
        
        
        
        % Calculate expected welfare given reg
        function EW = Ewelfare(obj,mu0,gamma)
            % INPUTS 
            % mu0: regulation intercept
            % gamma: regulation slope  
            
            % Normalize using welfare from social optimum
            [Wo,Qo,Eo] = socialopt(obj,0);
            
            % integrand
            ifun = @(x) regWQE(obj,mu0,gamma,x) .* normpdf(x,0,obj.sigma) / Wo;
            % integral.  Bounds are 5 sd's on either side of mean eta and F
            EW = integral(ifun,-obj.sigma*5,obj.sigma*5,'AbsTol',1e-12) * Wo;
        end
        

        
        % Value of FOC for mu0, given slope gamma. Utility function.
        function Val = FOC_mu_int(obj,mu,gamma)
            % INPUTS
            % mu: regulation intercept
            % gamma: regulation slope
            % flag: 2 = reg binds for high eta; 3 = reg binds for low eta
            
            % Calculates integral that is part of FOC_mu0
            % Only applies for standard that only binds for some values of eta
            
            % Normalize by change in welfare from unconstrained to social
            % opt
            [Qu,Eu] = uncon(obj,0);
            [Wo,Qo,Eo] = socialopt(obj,0);
            Norm = Wo / (Eu-Eo);
            
            % Find hat: value of eta where reg just binds
            [flag,hat] = etahat(obj,mu,gamma);
            
            etallim = obj.Z(1); etahlim = obj.Z(obj.Nz);    % integration limits
            if flag==2
                etallim = max([hat; etallim]);   % lower limit
            else
                etahlim = min([hat; etahlim]);   % upper limit
            end
            
            % Check limits
            if etallim >= etahlim
                % mu is so high that regulation never binds within data range
                Val = obj.phi * obj.BQQ;
            else
                % Define integrand
                ifun = @(x) ((obj.BQE+obj.BEE*gamma)*obj.BQeta*(x-hat)...
                    + obj.phi*(obj.BQE*gamma+obj.BQQ)) .* normpdf(x,0,obj.sigma) / Norm;

                % Compute integral
                Val = integral(ifun,etallim,etahlim,'AbsTol',1e-9);
            end
        end
        
        
        
        % Find optimal mu0 given gamma, and calc optimal expected welfare
        function [EW, mustar] = optmu(obj,gamma)
            % INPUTS
            % gamma: regulation slope
            
            % Finds optimal regulation intercept mustar given a reg slope gamma
            
            % Obtain unconstrained choices
            [Qu,Eu] = uncon(obj,obj.Z);       
            
            % Initial guess: mustar is shift down from unconstrained line
            % by phi / B_EE
            e = Eu(1); q = Qu(1);
            es = e + obj.phi / obj.BEE;     % shift down
            guess = es - gamma * q;           
             
            % Procedure different if gamma has same slope as unconstrained choice line
            slope = -obj.BQE / obj.BEE;
            diff = gamma - slope;           
            if abs(diff) < 1E-10        % slopes are same
                mustar = guess;         % mustar is the guess
            else
                options = optimset('TolX',1e-7);
                % Define utility function that calls FOC_mu integral
                util = @(x,obj,gamma) FOC_mu_int(obj,x,gamma);
                mustar = fzero(util,guess,options,obj,gamma);
            end
            
            % optimized welfare
            EW = Ewelfare(obj,mustar,gamma);
        end
        
        
        
        % Find optimal regulation slope
        function [EW, gammastar, mustar] = optgamma(obj)
                
            % Define utility function to minimize
            util = @(x) -optmu(obj,x);
            
            % Initial guess: solve analytically assuming standard always binds
            % Roots of polynomial
            VarEta = obj.sigma^2;
            a = obj.BQE * obj.phi^2;
            b = obj.BEE * obj.BQeta^2 * VarEta + obj.BQQ * obj.phi^2;
            c = obj.BQE * obj.BQeta^2 * VarEta;
            groots = roots([a b c]);
            guess = min(groots);        % low root is the local max; high root is local min
            
            % minimize
            tol = 1E-6 * c / obj.BQQ^2;
            options = optimoptions('fminunc','OptimalityTolerance',tol);
            gammastar = fminunc(util,guess,options);
            
            % Report mustar and welfare
            [EW, mustar] = optmu(obj,gammastar);
        end
        
        
        
        % Find optimal intensity standard (mu0 = 0)
        function [EW, gammastar] = optintensity(obj)
            
            % Initial guess is the max of the unrestricted slope or the
            % line going through the unrestriced line at Eeta
            % Optimum slope is less than this guess.
            [Q,E] = uncon(obj,0);
            guess = max([-obj.BQE/obj.BEE; E/Q]);
            
            % Search for optimum manually; algorithms tend to get stuck at
            % high values of gamma where reg never binds
            % Start coarse
            grid = linspace(0,1,1000)';
            grid = grid * guess;
            Wgrid = zeros(length(grid),1);
            for i = 1:length(grid)
                Wgrid(i) = Ewelfare(obj,0,grid(i));
            end
            
            % Fine search, starting with best from coarse search
            m = max(Wgrid);
            ind = find(Wgrid==m);
            if ind==length(grid)        % initial guess was best
                grid2 = linspace(grid(ind-1),grid(ind),5000)';
            elseif ind==1               % slope of zero was best
                grid2 = linspace(grid(ind), grid(ind+1),5000)';
            else
                grid2 = linspace(grid(ind-1),grid(ind+1),10000)';
            end
            Wgrid2 = zeros(length(grid2),1);
            for i = 1:length(grid2)
                Wgrid2(i) = Ewelfare(obj,0,grid2(i));
            end   
            
            EW = max(Wgrid2);
            gammastar = grid2(Wgrid2==EW);
        end
        
        
        
        % Collect all results and plots
        function [EWo, EW0, EWi, EWs, Q, E, MU, GAMMA] = resultsall(obj)
            
            % Return expected welfare for social optimum, flat std, intenstity standard,
            % and optimally indexed standard. P
            % Return Q and E data for plots
            
            % Use 95% c.i. for plot range
            out = floor(obj.Nz * 0.025);
            range = obj.Z(out+1:obj.Nz-out);
            
            % Get unconstrained and socially optimal Q * E
            [Qu,Eu] = uncon(obj,range);
            [Wo,Qo,Eo] = socialopt(obj,range);
            EWo = mean(Wo);     % socially optimal welfare (indep of eta)
            
            % Get optimal flat standard
            [EW0, mustar0] = optmu(obj,0);
            [W0,Q0,E0] = regWQE(obj,mustar0,0,range);   
            
            % Get optimal intensity standard
            [EWi, gammastari] = optintensity(obj);
            [Wi,Qi,Ei] = regWQE(obj,0,gammastari,range);
            
            % Get optimal indexed standard
            [EWs, gammastars, mustars] = optgamma(obj);
            [Ws,Qs,Es] = regWQE(obj,mustars,gammastars,range);
            
            % Return Q's and E's
            Q = [Qu Qo Q0 Qi Qs]; E = [Eu Eo E0 Ei Es];     
            
            % Return mu's and gamms
            MU = [mustar0 0 mustars]; GAMMA = [0 gammastari gammastars];
        end
        
    end
end






























