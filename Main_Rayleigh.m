% This file can be used to generate the variable-length stop-feedback 
% results on Rayleigh fading in the paper:
% "Short-Packet Transmission via Variable-Length Codes in the Presence of Noisy Stop Feedback"
% by J. ?stman, R. Devassy, G. Durisi, and E. G. Str?m

addpath ./Classes % include objects used for simulations

% Set up parameters
k=30; % number of information bits
ntot = 40; % channel uses per transmission round (feedback and forward)
d_max=400; % maximum allowed channel usages
SNR = 10; % forward transmission SNR 
SNR_f=10; % Feedback transmission SNR
bler_target = 1e-3; % Block-error rate target
np=6; % number of pilot symbols (if applicable)
nf=4; % number of channel uses for feedback transmission
MC_realizations=1e4; % Number of monte-carlo realizations
q = 0; % time-sharing value

%%--------------------------------------------------
%           Simulation of Rayleigh fading
%---------------------------------------------------
% Assumptions:
% - QPSK inputs
% - Pilot-assisted transmission of length np
% - Feedback decoding based on energy detection
% - One coherence block per transmission round
%%--------------------------------------------------
rayl_obj = Rayleigh(k, ntot, nf, np, d_max, SNR, SNR_f, bler_target, MC_realizations); % Create Rayleigh-fading object
info_density = rayl_obj.generateInfoDensity(); % compute information density
fb_threshold_candidates = rayl_obj.generatePossibleFeedbackThresholds(); % Generate valid feedback thresholds

ell_a = inf; % initial value for minimum average service time
for i = 1:length(fb_threshold_candidates)
    [psc, pcs] = rayl_obj.feedbackErrorProbability(fb_threshold_candidates(i)); % compute ACK/NACK error probabilities
    ell = VLSFHelper.findSmallestAverageNumberOfTransmissions(rayl_obj, info_density, pcs, psc, q); % Compute average # transmissions based on Thm. 1
    ell_a = min(ell_a, ell*ntot); % save smallest service time
end
disp(['l_a = ' num2str(ell_a) ' is achievable'])

