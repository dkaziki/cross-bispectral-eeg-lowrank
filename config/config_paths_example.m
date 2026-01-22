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


% Path to real EEG dataset (NOT included in this repository)
cfg.path_realdata = '/path/to/data_export';


% Path to subject file list (e.g., allnames.mat; NOT included)
cfg.path_allnames = '/path/to/allnames.mat';

