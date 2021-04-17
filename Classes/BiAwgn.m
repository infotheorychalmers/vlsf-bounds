classdef BiAwgn < SimulationObject
    properties
    end
    methods
        % Constructor
        function obj = BiAwgn(k, n_tot, nf, d_max, SNR, SNR_f, bler_target, MC_realizations)
            obj@SimulationObject(k, n_tot, nf, d_max, SNR, SNR_f, bler_target, MC_realizations);
        end
        
        % Generate information density samples according to Eq. 19
        % Note: s=1 in the paper
        function info_density = generateInfoDensity(obj, s)
            if nargin < 2
                s=1;
            end
            info_dens = @(x) cumsum(log(2) - log(1+exp(-2*x)), 2);
            mu = s*obj.rho;
            Z = sqrt(s^2*obj.rho)*randn(obj.MC_realizations,obj.n * obj.ell_m) + mu;
            info_density = info_dens(Z);
            info_density = info_density(:,obj.n:obj.n:end);
        end
        
        % Generate list of valid thresholds for the feedback link
        function threshold_list = generatePossibleFeedbackThresholds(obj)
            threshold_list = linspace(-sqrt(obj.rho*obj.nf), sqrt(obj.rho*obj.nf), 50); % threshold should be between S and C messages
        end
        
        % Generate feedback error for given nf and threshold
        function [psc, pcs] = feedbackErrorProbability(obj, thres)
            if obj.nf == 0
                pcs=0;psc=0;
            else
                psc = qfunc(sqrt(obj.nf*obj.rho_f) + thres);
                pcs = qfunc(sqrt(obj.nf*obj.rho_f) - thres);
            end
        end
    end
end