% Static functions to evaluate the VLSF, FLNF, and VLNSF performance
classdef VLSFHelper
    properties
        
    end
    methods (Static)
        % Find the smallest average number of transmissios according to
        % Thm. 1.
        % param: obj - Rayleigh or BiAwgn object
        % param: info_density - matrix of information density values
        % param: pcs - probability of stop -> continue error
        % param: pcs - probability of continue -> stop error
        % param: q - time sharing value
        % output: ell - minimum average number of transmissions
        function ell = findSmallestAverageNumberOfTransmissions(obj, info_density, pcs, psc, q)
            gamma_dec_steplength = 1;
            gamma_dec = 0;
            bler = 1;
            pmf_tau = [];
            %loop to find smallest threshold such that bler target is achievable
            while bler >= obj.bler_target
                [pmf_tau, tau] = VLSFHelper.getStoppingTime(obj, info_density, gamma_dec); %compute stopping time
                if isa(obj, 'BiAwgn')
                    bler = q + (1-q)*VLSFHelper.computeBler(obj, tau, pmf_tau, info_density, pcs); % Eq. 12 using Eq. 21
                elseif isa(obj, 'Rayleigh')
                    bler = q + (1-q)*VLSFHelper.computeBlerRelaxed(obj, pmf_tau, pcs, gamma_dec); % Eq. 12 using Eq. 27
                else
                    disp('Object not recognized.')
                    ell=inf;
                    return;
                end
                gamma_dec = gamma_dec + gamma_dec_steplength;
                if ~sum(pmf_tau) % when threshold is so large that tau > ell_m w.p. 1, stop searching
                    break;
                end
            end
            if bler <= obj.bler_target
                ell = (1-q)*VLSFHelper.computeAverageNumberOfTransmissions(obj, psc, pcs, pmf_tau); % Eq. 13
            else
                ell = inf;
            end
        end
    end
    
    methods (Access = private, Static)
        % Compute stopping time for given threshold
        function [pmf, tau] = getStoppingTime(obj, info_density, threshold)
            tau = (obj.ell_m+1)*ones(1,obj.MC_realizations); %assume that all tau are more than ell_m
            [ok, idx] = max(info_density >= threshold, [], 2); %pick out the index of the first time info dens crosses gamma. ok is a mask saying that we found one index crossing gamma or not.
            tau(ok) = idx(ok); %if there was a crossing, fill in the stopping time
            k = (1:obj.ell_m)';
            pmf = sum(bsxfun(@eq, k, tau),2)/numel(ok); % create PMF
        end
        
        %Compute relaxed version of average BLER Eq. 12 using Eq. 21
        function bler = computeBler(obj, tau, pmf_tau, info_density, pcs)
            M = 2^obj.k;
            ccdf_tau = 1-cumsum(pmf_tau); %Pr[tau > k] for k = 1,..., Nmax
            pr_taubar_tau = VLSFHelper.prob_tau_taubar(obj, tau, info_density);%evaluate upper bound on Pr[tau_bar <= tau, tau_bar = j]
            nu = 1:obj.ell_m;
            xi = (1-pcs).^(nu);
            alpha = [pcs*ones(1,obj.ell_m-1), 1];
            bler = min(sum(xi.*( alpha .* ccdf_tau' + (M-1).*pr_taubar_tau)),1);
        end
        
        %Compute relaxed version of average BLER Eq. 12 using Eq. 27
        function bler = computeBlerRelaxed(obj, pmf_tau, pcs, gamma_dec)
            M = 2^obj.k;
            alpha = @(k,pcs) (1-pcs).^(k-1) .* (pcs*(k <= (obj.ell_m-1)) + (k == obj.ell_m));
            ccdf_tau = 1-cumsum(pmf_tau); %Pr[tau > k] for k = 1,..., ell_m
            alpha_sum = sum(alpha((1:obj.ell_m)', pcs).*ccdf_tau); % evaluate sum_k phi_k*Pr[tau>k]
            v=1:obj.ell_m;
            xi = sum((1-pcs).^(v-1)); % xi_v
            bler = min(xi*(M-1)*exp(-gamma_dec) + alpha_sum, 1); % Eq. 12
        end
        
        %Compute the average number of transmissions according to Eq. 11
        function ell = computeAverageNumberOfTransmissions(obj, psc, pcs, pmf_tau)
            ccdf_tau = 1-cumsum([0; pmf_tau]); %Pr[tau > k] for k = 0,..., obj.ell_m
            Gv = nan(obj.ell_m,1); % Compute Gv according to Eq. 13
            for v = 1:obj.ell_m
                k = 1:v-1;
                Gv(v) = sum(k.*(1-pcs).^(k-1) .* pcs);
                k = v:obj.ell_m-1;
                Gv(v) = Gv(v) + sum(k.* psc.^(k-v).*(1-pcs)^(v-1)*(1-psc)) + obj.ell_m * psc^(obj.ell_m-v)*(1-pcs)^(v-1);
            end
            Gv = [0; Gv; Gv(end)];
            dH = diff(Gv);
            ell = sum(dH.*ccdf_tau); % Eq. 11
        end
        
        % Compute Eq. 21
        function prob = prob_tau_taubar(obj, tau, info_density)
            %iv) upper bound Pr[tau_bar <= tau, tau_bar = j]
            K = length(tau);
            tau(tau > obj.ell_m) = nan; %Insert nans for those taus that did not cross threshold, this is beacuse there are no indices for these tau in info dens.
            indx = sub2ind(size(info_density), 1:K, tau); %we will get nan if we never cross
            indx(isnan(indx)) = [];
            tau(isnan(tau)) = [];
            S_values_after_crossing = info_density(indx)';
            prob = nan(1, obj.ell_m);
            for j = 1:obj.ell_m
                mask = bsxfun(@eq, tau', j); %find all tau equal to j
                if sum(mask)==0
                    prob(j)=0;
                else
                    prob(j) = mean(exp(-S_values_after_crossing).*mask);
                end
            end
        end
    end
end