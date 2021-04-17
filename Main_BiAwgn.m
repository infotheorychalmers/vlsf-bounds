% This file can be used to generate the variable-length stop-feedback 
% results on bi-awgn in the paper:
% "Short-Packet Transmission via Variable-Length Codes in the Presence of Noisy Stop Feedback"
% by J. ?stman, R. Devassy, G. Durisi, and E. G. Str?m

addpath ./Classes

% Set up parameters
k=30; % number of information bits
ntot = 40; % channel uses per transmission round (feedback and forward)
d_max=400; % maximum allowed channel usages
SNR = 0; % forward transmission SNR 
SNR_f=0; % Feedback transmission SNR
bler_target = 1e-2; % Block-error rate target
nf=4; % number of channel uses for feedback transmission
MC_realizations=1e4; % Number of monte-carlo realizations
q = 0; % time-sharing value

%%--------------------------------------------------
%           Simulation of Bi-AWGN
%---------------------------------------------------
% Assumptions:
% - Feedback decoding based on threshold detection
%%--------------------------------------------------
biawgn_obj = BiAwgn(k, ntot, nf, d_max, SNR, SNR_f, bler_target, MC_realizations); % Create bi-awgn object
info_density = biawgn_obj.generateInfoDensity(); % compute information density
fb_threshold_candidates = biawgn_obj.generatePossibleFeedbackThresholds(); % Generate valid feedback thresholds

ell_a = inf; % initial value for minimum average service time
for i = 1:length(fb_threshold_candidates)
    [psc, pcs] = biawgn_obj.feedbackErrorProbability(fb_threshold_candidates(i)); % compute ACK/NACK error probabilities
    ell = VLSFHelper.findSmallestAverageNumberOfTransmissions(biawgn_obj, info_density, pcs, psc, q); % Compute average # transmissions based on Thm. 1
    ell_a = min(ell_a, ell*ntot); % save smallest service time
end
disp(['l_a = ' num2str(ell_a) ' is achievable'])

