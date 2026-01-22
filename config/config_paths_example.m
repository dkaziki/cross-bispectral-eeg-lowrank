% config_paths_example.m
%
% Example configuration file for the repository
% Copy this file to "config_paths.m" and edit the paths below
%
% This file defines the struct "cfg" used by the simulation scripts.

cfg = struct();

% Path to the METH toolbox (contains sa_eeg.mat, mri.mat, etc.)
cfg.path_meth = '/path/to/meth';

% Path to the bsfit toolbox
cfg.path_bsfit = '/path/to/bsfit';

