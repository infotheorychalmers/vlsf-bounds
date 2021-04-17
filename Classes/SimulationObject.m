% This class contains all the parameters in commom for BSC, BIAWGN, and
% Rayleigh objects.
classdef SimulationObject
    properties
        k % Number of information bits
        n_tot % Total number of channel uses in a round
        nf % Number of channel uses in the feedback transmission
        n % Number of channel uses in forward transmission
        ell_m
        d_max % Maximum delay in channel uses
        SNR % SNR value
        rho % SNR in linear domain
        SNR_f % SNR value for feedback transmission
        rho_f % SNR value for feedback (linear)
        bler_target % target block-error probability
        MC_realizations % number of monte-carlo realizations
    end
    methods
        %Constructor
        function obj = SimulationObject(k, n_tot, nf, d_max, SNR, SNR_f, bler_target, MC_realizations)
            obj.k = k;
            obj.n_tot = n_tot;
            obj.nf = nf;
            obj.n = n_tot-nf;
            obj.d_max = d_max;
            obj.SNR = SNR;
            obj.SNR_f = SNR_f;
            obj.MC_realizations = MC_realizations;
            obj.bler_target=bler_target;
            obj.rho = obj.db2lin(SNR);
            obj.rho_f = obj.db2lin(SNR_f);
            obj.ell_m = d_max / n_tot;
        end
    end
    %
    methods
        % Generate N(0,1) samples
        function x = randcn(obj, rows,cols)
            x = sqrt(0.5)*(randn(rows,cols) + 1i*randn(rows,cols));
        end
        % Convert between decibels and linear
        function x_lin = db2lin(obj, x_db)
            x_lin = 10.^(x_db./10);
        end
        % Convert between linear and decibels
        function x_db = lin2db(obj, x_lin)
            x_db = 10*log10(x_lin);
        end
    end
end