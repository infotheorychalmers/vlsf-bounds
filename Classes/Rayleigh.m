classdef Rayleigh < SimulationObject
    properties
        L % number of diversity branches (L=1 is considered in paper but it is easily extended)
        np % number of pilot symbols
    end
    methods
        % Constructor
        function obj = Rayleigh(k, n_tot, nf, np, d_max, SNR, SNR_f, bler_target, MC_realizations)
            obj@SimulationObject(k, n_tot, nf, d_max, SNR, SNR_f, bler_target, MC_realizations);
            obj.np = np;
            obj.L=1;
        end
        % Compute information density according to Eq. 26
        % Note: s=1 in the paper
        function info_density = generateInfoDensity(obj, s)
            if nargin < 2
                s=1;
            end  
            nd = obj.n - obj.np; %number of data symbols
            info_density = zeros(obj.MC_realizations, obj.ell_m); % allocate for information density
            sigma_e = 1/(obj.np*obj.rho); %variance of fading estimation
            % Generate monte carlo estimates
            parfor k = 1:obj.MC_realizations
                second_term = nan(nd, obj.ell_m);
                for l = 1:obj.L
                    S = sqrt(obj.rho)*[1+1i, -1+1i, -1-1i, 1-1i]/sqrt(2);   % QPSK constellation
                    X = S(randi([1,4],nd, obj.ell_m));                      % X is arranged as nd x Nmax
                    H = obj.randcn(1,obj.ell_m);                                     % H is arranged as 1 x Nmax
                    H_hat = sqrt(sigma_e)*obj.randcn(1,obj.ell_m) + H;          % H_hat is arranged as 1 x Nmax
                    N = obj.randcn(nd, obj.ell_m);                              % N is arranged as nd x Nmax
                    Y = X*diag(H) + N;                                      % Y is arranged as nd x Nmax
                    temp = S.'*H_hat;                                       %used in expectation term
                    
                    first_term = -s*abs(Y - X*diag(H_hat)).^2;
                    for v = 1:obj.ell_m
                        second_term(:,v) = -log((1/4)*sum(exp( -s * abs(bsxfun(@minus, Y(:,v), temp(:,v).'))),2));
                    end
                    info_density(k,:) = info_density(k,:) + sum(first_term + second_term, 1);
                end
            end
            info_density = cumsum(info_density,2); %compute accumlated information density
        end
        
        % Generate list of valid thresholds for the feedback link
        % Smallest pcs considered is bler*1e-3
        function threshold_list = generatePossibleFeedbackThresholds(obj)
            pcs_min = obj.bler_target * 1e-3;
            thres = -log(pcs_min)*obj.nf; % thresholds required for pcs_min according to Eq. 33
            threshold_list = linspace(0, thres, 64); % Generate 64 values
        end
        
        % Generate feedback error for given nf and threshold
        function [psc, pcs] = feedbackErrorProbability(obj, thres)
            if obj.nf == 0
                pcs=0;psc=0;
            else
                pcs = exp(-thres/obj.nf); % Eq. 33
                psc = 1-exp(-thres /(obj.nf*(obj.rho_f*obj.nf + 1))); %Eq. 31
            end
        end
    end
end