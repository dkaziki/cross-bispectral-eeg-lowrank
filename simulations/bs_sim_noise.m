%% BS Main Procedure for SIM DATA with SNR Variation
% bs_sim_noise.m
% Simulation: effect of noise / SNR on low-rank cross-bispectral decomposition
%
% Manuscript: "Low-Rank Tensor Decomposition for Cross-Bispectral Analysis of EEG Data"
%
% Dependencies:
%  - METH toolbox: https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html
%  - bsfit toolbox: https://github.com/guidonolte/bsfit
%
% How to run:
%  1) Copy config/config_paths_example.m -> config/config_paths.m
%  2) Edit paths in config/config_paths.m
%  3) Run this script

% Note: projecting voxel-wise source noise to sensor space is computationally intensive.
% For quick tests, reduce n_voxels or the time length (10000).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);

thisDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(thisDir);

% Load user config (defines struct 'cfg')
cfgFile = fullfile(repoDir, 'config', 'config_paths.m');
assert(exist(cfgFile,'file')==2, ...
    'Missing config/config_paths.m. Copy config_paths_example.m to config_paths.m and edit paths.');
run(cfgFile);

% Add dependencies to MATLAB path
addpath(cfg.path_meth);
addpath(cfg.path_bsfit);

% Load forward model and MRI template (from METH toolbox)
load(fullfile(cfg.path_meth, 'examples', 'sa_eeg.mat'));
load(fullfile(cfg.path_meth, 'templates', 'mri.mat'));

% Example Parameters

sr = 256;
segleng = 256;
segshift = 128;
epleng = 512;
nchan = 61;
maxfreqbins = 51;
n_model = 2;

snr_dB_levels = -15:1:15;
num_runs = 10;

mean_err = zeros(length(snr_dB_levels),1);
std_err = zeros(length(snr_dB_levels),1);
mean_loc = zeros(length(snr_dB_levels),1);
std_loc = zeros(length(snr_dB_levels),1);

idx1 = 4624; % source 1
idx2 = 239;  % source 2

for s = 1:length(snr_dB_levels)
    snr_db = snr_dB_levels(s);
    fprintf('\n==== SNR %d dB ====\n', snr_db);
    
    run_err = zeros(num_runs,1);
    run_loc = zeros(num_runs,1);

    for r = 1:num_runs
        
        s1 = data2filtdata(randn(10000,1), sr, [10,10], sr+1);
        s2 = s1.^2; s2 = s2 - mean(s2); s2 = s2 / max(abs(s2));

        s1_amp = 300; s2_amp = 300;
        s1 = s1_amp * s1;
        s2 = s2_amp * s2;

        %v1 = randn(3,1);
        v1 = [0; 1; 0];
        v1 = v1 / norm(v1);
        
        %v2 = randn(3,1);
        v2 = [1; 0; 0];
        v2 = v2 / norm(v2);

        dipole1 = v1' .* s1;
        dipole2 = v2' .* s2;

        leadfield1 = squeeze(sa.V_medium(:, idx1, :));
        leadfield2 = squeeze(sa.V_medium(:, idx2, :));

        S1 = leadfield1 * dipole1';
        S2 = leadfield2 * dipole2';
        sources_sim = S1' + S2';

       
        
% Source Space Noise uncorrelated

% --- Source-space noise: INDEPENDENT per voxel AND per direction ---
n_voxels = size(sa.V_medium, 2);

% noise(t, voxel, direction)
source_noise = randn(10000, n_voxels, 3);

% project to sensor space
sensor_noise = zeros(10000, nchan);
for v = 1:n_voxels
    for d = 1:3
        L = squeeze(sa.V_medium(:, v, d));     
        sensor_noise = sensor_noise + source_noise(:, v, d) * L';
    end
end


alpha_noise = data2filtdata(sensor_noise, sr, [8,12], sr+1);
beta_noise  = data2filtdata(sensor_noise, sr, [20,30], sr+1);


background = 0.5 * alpha_noise + 0.5 * beta_noise;

%  ---------------------------SNR------------------------------------------

% --- scale noise to target SNR ---
        signal_power = mean(sources_sim(:).^2);
        noise_power = mean(background(:).^2);
        snr_linear = 10^(snr_db / 10);
        scale_factor = sqrt(signal_power / (snr_linear * noise_power));
        background_scaled = background * scale_factor;

        
        data = sources_sim + background_scaled;
 

        [cs, csnr, ~] = data2bs_univar(data, segleng, segshift, epleng, maxfreqbins, []);
        xx = abs(cs ./ csnr);

        max_bicoherence = 0;
        for ch = 1:nchan
            for f1 = 9:14
                for f2 = 9:14
                    b = xx(ch,f1,f2);
                    if b > max_bicoherence
                        max_bicoherence = b;
                        max_channel = ch;
                        max_f1 = f1;
                        max_f2 = f2;
                    end
                end
            end
        end

        dominant_freq = max_f1;   % or max_f2
	freqpairs = [dominant_freq dominant_freq];

        para = [];
        [bs, ~] = data2bs_event(data, segleng, segshift, epleng, freqpairs, para);

        
        [a, ~, err, ~, ~] = bsfit(bs, n_model, para);

        
        lead = sa.V_medium;
        A = mkfilt_eloreta(lead, 0.5);
        [~, nvoxel, ndum] = size(A);
        F = zeros(nvoxel, ndum, n_model);
        for i = 1:n_model
            for k = 1:ndum
                F(:, k, i) = A(:,:,k)' * a(:, i);
            end
        end
        [Fout, ~] = moca_ncomp(F);

        % --- Localization Error ---
        true_voxels = [idx1, idx2];
        loc_err = 0;
        for i = 1:n_model
            Floc = sqrt(sum(Fout(:,:,i).^2, 2));
            [~, peak_voxel] = max(Floc);
            dists = sqrt(sum((sa.grid_medium(peak_voxel,:) - sa.grid_medium(true_voxels,:)).^2, 2));
            [min_dist, ~] = min(dists);
            loc_err = loc_err + min_dist;
        end
        loc_err = loc_err / n_model;

        run_err(r) = err;
        run_loc(r) = loc_err;
    end

   
    mean_err(s) = mean(run_err);
    std_err(s) = std(run_err);
    mean_loc(s) = mean(run_loc);
    std_loc(s) = std(run_loc);
end


% For log Scale

epsilon = 1e-3;
mean_loc(mean_loc == 0) = epsilon;
std_loc(std_loc == 0) = epsilon;

%% Plot with Error Bars

figure;

% --- Model Error Plot ---

subplot(2,1,1);
errorbar(snr_dB_levels, mean_err, std_err, '-o', ...
    'Color', [0 0.45 0.74], 'MarkerFaceColor', [0 0.45 0.74], ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 8);
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Model Error', 'FontSize', 12);
set(gca, 'YScale', 'log');
title('Model Error vs. SNR', 'FontSize', 13);
set(gca, 'FontSize', 11);

% --- Localization Error Plot ---
subplot(2,1,2);
errorbar(snr_dB_levels, mean_loc, std_loc, '-o', ...
    'Color', [0.85 0.33 0.1], 'MarkerFaceColor', [0.85 0.33 0.1], ...
    'LineWidth', 2, 'MarkerSize', 6, 'CapSize', 8);
grid on;
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Localization Error (mm)', 'FontSize', 12);
set(gca, 'YScale', 'log');
title('Localization Error vs. SNR', 'FontSize', 13);
set(gca, 'FontSize', 11);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
