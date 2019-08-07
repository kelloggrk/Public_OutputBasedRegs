% The class that instantiates the indexmodelF object 
% Subclass of indexmodel that also allows for fuel price uncertainty
classdef indexmodelF < indexmodel
    properties
        % We specify here all of the properties we want to have
        % attached to each object (except those already defined in the
        % superclass)
        F0
        sigmaF
        sigmaeta
        rho
        SIGMA
        NzF
        ZF
    end
    methods
        % Constructor
        function obj = indexmodelF(params)
            % Object defns from superclass
            obj = obj@indexmodel(params);
            
            % Fuel price volatility specific parameters
            obj.F0 = params.F0;             % initial (mean) fuel price
            obj.sigmaeta = params.sigmaeta; % volatility of demand shocks
            obj.sigmaF = params.sigmaF;     % volatility of fuel prices
            obj.rho = params.rho;           % covar between eta and F
            
            % Covariance matrix for eta and F
            obj.SIGMA = [obj.sigmaeta^2 obj.rho*obj.sigmaeta*obj.sigmaF;...
                obj.rho*obj.sigmaeta*obj.sigmaF obj.sigmaF^2];  
          
            % Set up normal distribution for F: mean, std dev, and
            % nodes for integration
            obj.NzF = 1000;                     % # of nodes
            ZF = linspace(0,1,obj.NzF+2);       % vector of ZFs to integrate over
            obj.ZF = norminv(ZF',obj.F0,obj.sigmaF);
            obj.ZF = obj.ZF(2:length(ZF)-1);    % drop first and last extremes (zero and infty)   
        end
        
        
        
        % eta and F ellipse that form a 95% confidence interval
        function [ploteta, plotF] = plotetaF(obj)
            % Create ellipse of eta and F that encompass a 95% confidence
            % interval in two dimensions. Ellipse follows countour of
            % iso-pdf (or iso-Mahalanobis distance)
            % See https://stackoverflow.com/questions/35222566/multivariate-normal-distribution-matlab-probability-area
            % See https://upload.wikimedia.org/wikipedia/commons/a/a2/Cumulative_function_n_dimensional_Gaussians_12.2013.pdf
            alpha = 0.05;               % confidence level
            r = sqrt(-2*log(alpha));    % Mahalanobis distance
            M = chol(obj.SIGMA);        % Choleski decomposition of SIGMA
            theta = 0:0.01:2*pi;        % radians to plot
            plotdata = r * [cos(theta)', sin(theta)'] * M + [0 obj.F0];     % [eta, F]
            ploteta = plotdata(:,1); plotF = plotdata(:,2);   
        end
        
        
        % Function for Q and E choices
        function [Q,E] = uncon(obj,eta,F)
            % Finds unconstrained Q and E as a function of eta and F
            Q = obj.EQ + (-obj.BEE*obj.BQeta*eta + obj.BQE*obj.BEF*(F-obj.F0)) / obj.SOC;
            E = obj.EE + (obj.BQE*obj.BQeta*eta - obj.BQQ*obj.BEF*(F-obj.F0)) / obj.SOC;
        end
            
        
        
        % Total welfare function - not normalized
        function W = Bnn(obj,Q,E,eta,F)
            % INPUTS 
            % Q: quantity or attribute value
            % E: energy use (or efficiency) value
            % eta: value of output shock
            % F: value of fuel price shock
            
            % Use facts that on the actual reg constraint, B_Q = -gamma0 *
            % B_E, and that at (Q0,E0,Eeta,F0), B_E = phiact
            
            % Define welfare using Taylor expansion around (Q0, E0, Eeta, F0)
            % First-order terms use shadow values of actual regulatory constraint
            Wp = obj.phiact*(E-obj.E0) - obj.gamma0*obj.phiact*(Q-obj.Q0)...
                + 1/2*obj.BEE*(E-obj.E0).^2 + 1/2*obj.BQQ*(Q-obj.Q0).^2....
                + obj.BQE*(Q-obj.Q0).*(E-obj.E0) + obj.BQeta*(Q-obj.Q0).*eta...
                + obj.BEF*(E-obj.E0).*(F-obj.F0);
            W = Wp - obj.phi*E;       % total welfare, including externality
        end
        
        
        
        % Total welfare function, normalized to equal zero at all
        % unconstrained choices
        function W = B(obj,Q,E,eta,F)
            % INPUTS 
            % Q: quantity or attribute value
            % E: energy use (or efficiency) value
            % eta: value of output shock
            % F: value of fuel price shock
            
            % Normalization: total welfare is equal to zero at
            % (Q0,E0,Eeta,F). Also equal to zero at all unconstrained 
            % choices, given values of eta and F
            [Qu,Eu] = uncon(obj,eta,F);           % unconstrained Q and E
            Wu = Bnn(obj,Qu,Eu,eta,F);            % welfare when unconstrained
            W = Bnn(obj,Q,E,eta,F) - Wu;          % normalized welfare
        end
        
        
        
        % Socially optimal Q, E, and welfare, given eta and F
        function [W,Q,E] = socialopt(obj,eta,F)
            % INPUTS 
            % Q: quantity or attribute value
            % E: energy use (or efficiency) value
            % eta: value of output shock
            % F: value of fuel price shock
            
            % Get unconstrained location and then shift by dQ/dphi and dE/dphi
            [Qu,Eu] = uncon(obj,eta,F);
            Q = Qu - obj.BQE / obj.SOC * obj.phi;
            E = Eu + obj.BQQ / obj.SOC * obj.phi;
            
            % Welfare
            W = B(obj,Q,E,eta,F);
        end

        
        
        % Calculate Q, E, and welfare given regulation, eta, and F
        function [W,Q,E] = regWQE(obj,mu0,gamma,eta,F)
            % INPUTS 
            % mu0: regulation intercept
            % gamma: regulation slope            
            % eta: value of output shock
            % F: value of fuel price shock
            
            [Qu,Eu] = uncon(obj,eta,F);     % Get unconstrained locations
            Q = Qu; E = Eu;                 % initialize Q and E with unconstrained values
            Ec = mu0 + gamma * Qu;          % Get constrained E's at same Q's as unconstrained
            
            % Identify cases where regulation binds
            Ind = find(Ec<Eu);
            
            % Change Q and E where constraint binds. Idea is to trace out
            % change in Q and E as the constraint shifts down vertically
            % from just binding at the unconstrained (Q,E) to its actual location
            if ~isempty(Ind)       % make sure not always unconstrained
                denom = obj.BEE*gamma^2 + 2*obj.BQE*gamma + obj.BQQ;    % denom of dQ/dmu0 and dE/dmu0
                dQdmu0 = -(obj.BEE*gamma+obj.BQE) / denom;             
                dEdmu0 = (obj.BQE*gamma+obj.BQQ) / denom;       
                delta = Ec - Eu;                                        % vertical change in mu0
                % Shift Q and E according to a shift in mu0 equal to delta
                Q(Ind) = Qu(Ind) + dQdmu0 * delta(Ind);
                E(Ind) = Eu(Ind) + dEdmu0 * delta(Ind);
            end
            
            % Now calculate welfare
            W = B(obj,Q,E,eta,F);
        end
        
        
        
        % Calculate expected welfare given reg
        function EW = Ewelfare(obj,mu0,gamma)
            % INPUTS 
            % mu0: regulation intercept
            % gamma: regulation slope  
            
            % Normalize using welfare from social optimum
            [Wo,~,~] = socialopt(obj,0,obj.F0);
            
            % integrand
            ifunreg = @(x,y) regWQE(obj,mu0,gamma,x,y);
            ifunpdf = @(x,y) reshape(mvnpdf([x(:) y(:)], [0 obj.F0], obj.SIGMA), size(x));
            ifun = @(x,y) ifunreg(x,y) .* ifunpdf(x,y) / Wo;

            % integral. Bounds are 5 sd's on either side of mean eta and F
            EW = integral2(ifun,-obj.sigma*5,obj.sigma*5,...
                obj.F0-obj.sigmaF*5,obj.F0+obj.sigmaF*5,'AbsTol',1e-7) * Wo;
        end
        
        
        
        % Find optimal mu0 given gamma, and calc optimal expected welfare
        function [EW, mustar] = optmu(obj,gamma)
            % INPUTS
            % gamma: regulation slope
            
            % Finds optimal regulation intercept mustar given a reg slope gamma
            
            % Upper bound intersects the unconstrained choice with highest
            % E. Use highest eta; lowest F if B_EF<0; otherwise use highest
            % eta and highest F
            if obj.BEF<0
                [Qmax,Emax] = uncon(obj,obj.sigma*5,obj.F0-obj.sigmaF*5);
            else
                [Qmax,Emax] = uncon(obj,obj.sigma*5,obj.F0+obj.sigmaF*5);
            end
            % Lower bound intersects the expected social optimum
            % Standard will pass through this point if optimal standard always binds
            [~,Qmin,Emin] = socialopt(obj,0,obj.F0);
            % Solve for min and max mu's
            mumax = Emax - gamma * Qmax; mumin = Emin - gamma * Qmin;
            
            % Do coarse grid search for optimal mu
            grid = linspace(0,1,25)';
            grid = mumin + grid * (mumax - mumin);
            Wgrid = zeros(length(grid),1);
            for i = 1:length(grid)
                Wgrid(i) = Ewelfare(obj,grid(i),gamma);     % welfare at each grid point
            end            
            
            % Precise search
            guess = grid(Wgrid==max(Wgrid));        % initial guess
            lowbnd = guess - (grid(2)-grid(1));     % lower bound
            upbnd = guess + (grid(2)-grid(1));      % upper bound
            [Wo,~,~] = socialopt(obj,0,obj.F0);     % to normalize tolerance
            tol = 1e-6 / Wo;
            options = optimset('TolFun',tol);
            util = @(x) -Ewelfare(obj,x,gamma);             % utility function to minimize
            mustar = fminbnd(util,lowbnd,upbnd,options);    % optimal mu
            
            % optimized welfare
            EW = Ewelfare(obj,mustar,gamma);
        end
        
        
        
        % Find optimal regulation slope
        function [EW, gammastar, mustar] = optgamma(obj)
                
            % Upper bound on slope is usually unconstrained dE/dQ under eta shifts
            gammamax = -obj.BQE / obj.BEE;
            
            % Do coarse grid search for optimal gamma (min guess is zero)
            ngrid = 50;
            grid = linspace(0,1,ngrid)';
            grid = grid * gammamax;
            diff = grid(2)-grid(1);                 % increment between grid points for gamma
            Wgrid = zeros(length(grid),1);
            for i = 1:length(grid)
                Wgrid(i) = optmu(obj,grid(i));      % welfare at each grid point
            end  

            % Keep searching if gammamax was best (can happen if rho<0)
            if Wgrid(ngrid)==max(Wgrid)
                diff = 2 * diff;                    % double increment
                Wold = max(Wgrid);                  % old max
                Wnew = Wold;                        % new max
                gamma = grid(ngrid);                % last gamma checked
                while Wnew>=Wold                    % loop until new max is lower
                    Wold = Wnew;                    % replace old value
                    gamma = gamma + diff;           % update gamma
                    Wnew = optmu(obj,gamma);        % welfare at new gamma
                end
                grid(ngrid+1) = gamma - diff;       % best gamma
                Wgrid(ngrid+1) = Wold;              % best welfare
            end
            
            % Precise search
            guess = grid(Wgrid==max(Wgrid));        % initial guess
            guess = guess(1);                       % make scalar
            lowbnd = guess - diff;                  % lower bound
            upbnd = guess + diff;                   % upper bound
            [Wo,~,~] = socialopt(obj,0,obj.F0);      % to normalize tolerance
            tol = 1e-5 / Wo;
            options = optimset('TolFun',tol);
            util = @(x) -optmu(obj,x);              % utility function to minimize
            gammastar = fminbnd(util,lowbnd,upbnd,options);    % optimal gamma
            
            % Report mustar and welfare
            [EW, mustar] = optmu(obj,gammastar);
        end
        
        
        
        % Find optimal intensity standard (mu0 = 0)
        function [EW, gammastar] = optintensity(obj)
            
            % upper bound is an intensity standard that barely binds
            % Get max of E/Q around the 95% c.i. of unconstrained choices
            [ploteta, plotF] = plotetaF(obj);           % eta and F for 95% c.i.
            [Qu,Eu] = uncon(obj,ploteta,plotF);         % unconstrained Q and E
            guess = max(Eu./Qu);
            
            % Do coarse grid search for optimal gamma (min guess is zero)
            grid = linspace(0,1,50)';
            grid = grid * guess;
            Wgrid = zeros(length(grid),1);
            for i = 1:length(grid)
                Wgrid(i) = Ewelfare(obj,0,grid(i));
            end
            
            % Precise search
            guess = grid(Wgrid==max(Wgrid));        % initial guess
            lowbnd = guess - (grid(2)-grid(1));     % lower bound
            upbnd = guess + (grid(2)-grid(1));      % upper bound
            [Wo,~,~] = socialopt(obj,0,obj.F0);      % to normalize tolerance
            tol = 1e-5 / Wo;
            options = optimset('TolFun',tol);
            util = @(x) -Ewelfare(obj,0,x);         % utility function to minimize
            gammastar = fminbnd(util,lowbnd,upbnd,options);    % optimal gamma
            
            % optimized welfare
            EW = Ewelfare(obj,0,gammastar);
        end
        
        
        
        % Collect all results and plots
        function [EWo, EW0, EWi, EWs, Q, E, EQ, EE, MU, GAMMA] = resultsall(obj)
            % Return expected welfare for social optimum, flat std, intenstity standard,
            % and optimally indexed standard. Return plot data.
            
            % Get 95% c.i. of eta and F for plotting
            [ploteta, plotF] = plotetaF(obj);
            
            % Get unconstrained and socially optimal Q * E
            [Qu,Eu] = uncon(obj,ploteta,plotF);
            [Wo,Qo,Eo] = socialopt(obj,ploteta,plotF);
            EWo = mean(Wo);     % socially optimal welfare (indep of eta and F)
            
            % Min and max of plotted Q's
            Qmin = min([Qu; Qo]); Qmax = max([Qu; Qo]); 
            
            % Get optimal flat standard
            [EW0, mustar0] = optmu(obj,0);
            % Plot line for standard
            Q0 = [Qmin; Qmax]; E0 = [mustar0; mustar0];
            
            % Get optimal intensity standard
            [EWi, gammastari] = optintensity(obj);
            % Plot line for standard
            Qi = [Qmin; Qmax]; Ei = [gammastari*Qmin; gammastari*Qmax];
            
            % Get optimal indexed standard
            [EWs, gammastars, mustars] = optgamma(obj);
            % Plot line for standard
            Qs = [Qmin; Qmax]; Es = mustars + [gammastars*Qmin; gammastars*Qmax];

            % Return Q's and E's
            Q.Qu = Qu; Q.Qo = Qo; Q.Q0 = Q0; Q.Qi = Qi; Q.Qs = Qs;
            E.Eu = Eu; E.Eo = Eo; E.E0 = E0; E.Ei = Ei; E.Es = Es;
            
            % Return mu's and gamms
            MU = [mustar0 0 mustars]; GAMMA = [0 gammastari gammastars];                
            
            % Return Q and E at expected eta and F
            EQ = obj.EQ; EE = obj.EE;      
        end
    end
end






























