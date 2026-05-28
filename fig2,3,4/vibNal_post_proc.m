function [mags, nfs, phis, magsC, phisC, gains, gainsC] = vibNal_post_proc(params, vib)
% VIBNAL_POST_PROC Snippet of code from scc_tcNalAll_mjsMod.m
%
% This function performs post-processing on mscan data (from Dewey's
% analysis code). Unfortunately in scc_tcNalAll_mjsMod.m this
% post-processed data isn't saved for later analysis. So if we want to do
% analysis on tuning curves we need to do these calculations again (or
% modify scc_tcNalAll_mjsMod.m even more than I already have, which I don't
% want to do right now).
%
% [mags, nfs, phis, magsC, phisC, gains, gainsC] = vibNal_post_proc(params, vib)
%
% @param[in] params - struct loaded from the params.mat file in the
% analysis* folder for a given MScan data set. This is originally saved by
% scc_mscanNalAll_mjsMod.m.
%
% @param[in] vib - struct loaded from the vibNal.mat file in the
% analysis* folder for a given MScan data set. This is originally saved by
% scc_mscanNalAll_mjsMod.m.
% 
%
% @author mjs (original code logic by Dewey).
%
% 2023-07-02
%
% Created.

% HARDWIRED PARAMETERS FROM DEWEYS SCRIPT
nsd = 3; % std above noise floor required for clean plot
phi_ref_f = 9000;

% Parameters
f1s = params.f1.f;
L1s = params.f1.L;
f1N = length(f1s);
L1N = length(L1s);

L1Pas = 2e-5 * sqrt(2) * 10.^(L1s/20);

trial_dur = params.trial_dur;
stim_dur = params.stim_dur;
trial_n = params.trial_n;
    
% Displacment data
mags = vib.f1.mag;
nfs = vib.f1.nf + nsd .* vib.f1.nfsd;
phis = unwrap(vib.f1.phi)/(2*pi);

% Clean data
magsC = mags;
magsC(mags<nfs)=NaN;
phisC = phis;
phisC(mags<nfs)=NaN;

% Align phases to highest level at reference frequency
if nanmedian(phisC(1:8, end)) > 0
    phis = phis-1;
    phisC = phisC-1;
end

[~,phi_ref_fi] = ismember(phi_ref_f, f1s);
phi_refs = phis(phi_ref_fi, :);
phi_diffs = phi_refs - phi_refs(:,end);
phi_diffs = round(phi_diffs,0);
phis = phis - phi_diffs;
phisC = phisC - phi_diffs;

% Calculate gain
gains = 20*log10(mags./L1Pas);
gainsC = 20*log10(magsC./L1Pas);
       
end

