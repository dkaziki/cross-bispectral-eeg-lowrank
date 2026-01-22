%% Simulation example pipeline 
% run_sim_example.m
%
% Note: This script reproduces a simulation pipeline used in the manuscript.
% It requires external toolboxes (METH + bsfit).
%
% Dependencies:
%  - METH toolbox: https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html
%  - bsfit toolbox: https://github.com/guidonolte/bsfit
%
%
% Setup:
%  1) Copy config/config_paths_example.m -> config/config_paths.m
%  2) Edit paths in config/config_paths.m
%  3) Run this script


% This script provides a configurable simulation pipeline used in the manuscript.
% It contains code blocks for injecting:
%   - 2 sources (default; enabled below)
%   - 4 sources (optional; currently commented)
%   - 6 sources (optional; currently commented)

% These scenarios were used to systematically assess the robustness and
% scalability of the proposed low-rank cross-bispectral decomposition.

% How to choose the number of sources:
%   1) Set n_sources = 2, 4, or 6.
%   2) Enable the corresponding “source injection” block by uncommenting it
%      (and keep the other blocks commented).
%   3) Update the "true_voxel_indices" and "dipole_dirs" lists accordingly
%      (these are used only for visualization/evaluation).
%

% Notes:
%   - The low-rank model order n (used in bsfit) is set separately and can be
%     adjusted depending on the scenario.
%   - Noise injection is provided as an optional block (commented by default).

% The voxel indices below specify the locations of simulated dipolar sources
% in the discretized source space defined by the METH forward model.
%
% These indices correspond to rows in `sa.grid_medium` and columns in
% `sa.V_medium`. They were selected to ensure spatial separation between
% sources while remaining within anatomically plausible cortical regions.
%
% The specific voxel values are fixed to allow reproducible simulations
% across runs and across different source-complexity scenarios.
%
% Note: The exact anatomical interpretation of voxel indices depends on the
% head model and source grid provided by the METH toolbox.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);

thisDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(thisDir);

% Load user config (defines struct 'cfg')
cfgFile = fullfile(repoDir, 'config', 'config_paths.m');
assert(exist(cfgFile,'file')==2, ...
    'Missing config/config_paths.m. Copy config_paths_example.m to config_paths.m and edit paths.');
run(cfgFile);

% Add dependencies
addpath(cfg.path_meth);
addpath(cfg.path_bsfit);

% Load forward model and MRI template (from METH toolbox)
load(fullfile(cfg.path_meth, 'examples', 'sa_eeg.mat'));
load(fullfile(cfg.path_meth, 'templates', 'mri.mat')); %#ok<NASGU>  % mri not used below, but kept for completeness



%% 
    %Parameters
    sr = 256; 
    segleng = 256; 
    segshift = 128;
    epleng = 512;
    nchan = 61;
    maxfreqbins = 51;

    n_sources = 2; % number of injected sources in this example
    snr_db = 0;    % set if using the noise block (optional)
    

    
%% ___________________SIMULATED DATA____________________________________________________________________________________________________________________


%% ______________________2 sources_________________________________________

% 
% for run_idx = 1:num_runs

% source injection voxel indices
dipole_idx1 = 4624;  % alpha source
dipole_idx2 = 239;   % coupled source

leadfield1 = squeeze(sa.V_medium(:, dipole_idx1, :));  % [nchan x 3]
leadfield2 = squeeze(sa.V_medium(:, dipole_idx2, :));  % [nchan x 3]

s1=data2filtdata(randn(10000,1),sr,[10,10],sr+1);

%s2=data2filtdata(randn(10000,1),sr,[11,11],sr+1);

% coupled signal
    s2 = s1.^2;
    s2 = s2 - mean(s2);            % zero mean
    s2 = s2 / max(abs(s2));        % normalize

% % % coupled signal
%     s3 = s1.^2;
%     s3 = s3 - mean(s3);            % zero mean
%     s3 = s3 / max(abs(s3));        % normalize

% amps
 s1_amp = 300;   % microV
 s2_amp = 300;
    
 s1 = s1_amp * s1;
 s2 = s2_amp * s2;

%-------------------------------

%random directions

%  dir1 = randn(3,1); 
%  dir2 = randn(3,1); 

 
 
%fixed directions a

 dir1=[-0.5742;0.2045;0.7927];
 dir2 =[0.9029;0.1150;0.4141];

%fixed directions b
% dir1 =[-0.2451;-0.7689;0.9828];
% dir2 =[0.2451;-0.5674;0.7861];




% % Orthogonal directions
% dir1 = [0; 1; 0];
% dir2 = [1; 0; 0];

% % Parallel directions
% dir1 = [1; 1; 1];
% dir2 = [2; 2; 2];  



%Normalize
dir1 = dir1 / norm(dir1);
dir2 = dir2 / norm(dir2);


dipole1 = dir1' .* s1;  
dipole2 = dir2' .* s2;

%------------------------------

    % forward projection 
    S1 = leadfield1 * dipole1';  
    S2 = leadfield2 * dipole2';

    sources_sim=zeros(10000,nchan);
    % store
    sources_sim(:, :) = S1' + S2' ;

    %data=sources_sim;
    

%---------------------------------------------------------------------
    
voxel_power1 = zeros(size(sa.grid_medium, 1), 1);
voxel_power2 = zeros(size(sa.grid_medium, 1), 1);

% power at each voxel
power1 = sqrt(mean(dipole1.^2, 2));  
power2 = sqrt(mean(dipole2.^2, 2));

% average across 3 orientations
voxel_power1(dipole_idx1) = mean(power1);
voxel_power2(dipole_idx2) = mean(power2);

% project to cortex
out_cortex1 = spatfiltergauss(voxel_power1, sa.grid_medium, 2, sa.cortex10K.vc);
out_cortex2 = spatfiltergauss(voxel_power2, sa.grid_medium, 2, sa.cortex10K.vc);

% figure;
% 
% subplot(1, 2, 1);
% showsurface(sa.cortex10K, [], out_cortex1);
% title('Source 1');
% view([0 90]);
% 
% subplot(1, 2, 2);
% showsurface(sa.cortex10K, [], out_cortex2);
% title('Source 2');
% view([0 90]);


% % ___________________4 sources_____________________________________________
% 
% 
% % % 
% % % voxel indices for S3 (alpha) and S4 (beta)
%  dipole_idx3 = 4678; 
%  dipole_idx4 = 47;   
% %  
% % % leadfields
%  leadfield3 = squeeze(sa.V_medium(:, dipole_idx3, :));  % [nchan x 3]
%  leadfield4 = squeeze(sa.V_medium(:, dipole_idx4, :));  % [nchan x 3]
% 
% s3=data2filtdata(randn(10000,1),sr,[12,12],sr+1);
%  
% %s4=data2filtdata(randn(10000,1),sr,[12,12],sr+1);
% % 
% % coupled signal
%     s4 = s3.^2;
%     s4 = s4 - mean(s4);            % zero mean
%     s4 = s4 / max(abs(s4));        % normalize
% 
% % amps
%  s3_amp = 300;   % microV
%  s4_amp = 300;
%     
%  s3 = s3_amp * s3;
%  s4 = s4_amp * s4;
% % 
% % 
% 
% 
% %% random dir
% 
% % dir3 = randn(3,1); 
% % dir3 = dir3 / norm(dir3);  
% % 
% % dir4 = randn(3,1); 
% % dir4 = dir4 / norm(dir4); 
% 
% 
% % fixed dir b
% 
%  dir3 =[-0.2451;-0.7689;0.9828];
%  dir4 =[0.2451;-0.5674;0.7861];
% 
% % %fixed dir a
% %  dir3=[-0.5742;0.2045;0.7927];
% %  dir4 =[0.9029;0.1150;0.4141];
% 
% 
% 
% dipole3 = dir3' .* s3;  
% dipole4 = dir4' .* s4;
% % 
% % %---------------------------------------------------
% % 
% %     % forward projection 
%     S3 = leadfield3 * dipole3';  % [nchan x epleng]
%     S4 = leadfield4 * dipole4';
% 
%     % store 
%     sources_sim(:, :) = sources_sim(:,:)+S3' + S4' ;
% 
%    % data=sources_sim;
%     
% 
% 
% % 
% % % _______________________________________ VISUALIZE ALL SOURCES _______________________________________
% 
% 
% voxel_power3 = zeros(size(sa.grid_medium, 1), 1);
% voxel_power4 = zeros(size(sa.grid_medium, 1), 1);
% 
% voxel_power3(dipole_idx3) = mean(sqrt(mean(dipole3.^2, 2)));
% voxel_power4(dipole_idx4) = mean(sqrt(mean(dipole4.^2, 2)));
% 
% out_cortex3 = spatfiltergauss(voxel_power3, sa.grid_medium, 2, sa.cortex10K.vc);
% out_cortex4 = spatfiltergauss(voxel_power4, sa.grid_medium, 2, sa.cortex10K.vc);
% 
% %%all 4 sources
% 
% figure;
% 
% subplot(2, 2, 1);
% showsurface(sa.cortex10K, [], out_cortex1);
% title('Source 1');
% view([0 90]);
% 
% subplot(2, 2, 2);
% showsurface(sa.cortex10K, [], out_cortex2);
% title('Source 2');
% view([0 90]);
% 
% subplot(2, 2, 3);
% showsurface(sa.cortex10K, [], out_cortex3);
% title('Source 3');
% view([0 90]);
% 
% subplot(2, 2, 4);
% showsurface(sa.cortex10K, [], out_cortex4);
% title('Source 4');
% view([0 90]);
% % 


% %% _________________________6 Sources____________________________________________
% % % 
% 
% dipole_idx5 = 4932; 
% dipole_idx6 = 195;  
% 
%
% leadfield5 = squeeze(sa.V_medium(:, dipole_idx5, :));  % [nchan x 3]
% leadfield6 = squeeze(sa.V_medium(:, dipole_idx6, :));  % [nchan x 3]
% 
%
% s5 = data2filtdata(randn(10000,1), sr, [11,11], sr+1); 
% s6 = data2filtdata(randn(10000,1), sr, [8,8], sr+1); 
% 
% %amps
% s5_amp = 100;
% s6_amp = 100;
% s5 = s5_amp * s5;
% s6 = s6_amp * s6;
% 
% 
% % dir5 = randn(3,1); dir5 = dir5 / norm(dir5);
% % dir6 = randn(3,1); dir6 = dir6 / norm(dir6);
% 
% %fixed dir
% 
% dir6 =[-0.1520;-0.9879;-0.0299]
% dir5 =[-0.1451;0.7689; 0.9828]
% 
% dipole5 = dir5' .* s5;
% dipole6 = dir6' .* s6;
% 
% 
% S5 = leadfield5 * dipole5';
% S6 = leadfield6 * dipole6';
% 
%
% sources_sim = sources_sim + S5' + S6';
% data = sources_sim;
% 
%
% voxel_power5 = zeros(size(sa.grid_medium, 1), 1);
% voxel_power6 = zeros(size(sa.grid_medium, 1), 1);
% voxel_power5(dipole_idx5) = mean(sqrt(mean(dipole5.^2, 2)));
% voxel_power6(dipole_idx6) = mean(sqrt(mean(dipole6.^2, 2)));
% 
% out_cortex5 = spatfiltergauss(voxel_power5, sa.grid_medium, 2, sa.cortex10K.vc);
% out_cortex6 = spatfiltergauss(voxel_power6, sa.grid_medium, 2, sa.cortex10K.vc);
% 
%
% figure('Position', [100, 100, 1600, 900]);  
% 
% % 6 sources
% subplot(2, 3, 1); showsurface(sa.cortex10K, [], out_cortex1); title('Source 1'); view([0 90]);
% subplot(2, 3, 2); showsurface(sa.cortex10K, [], out_cortex2); title('Source 2'); view([0 90]);
% subplot(2, 3, 3); showsurface(sa.cortex10K, [], out_cortex3); title('Source 3'); view([0 90]);
% subplot(2, 3, 4); showsurface(sa.cortex10K, [], out_cortex4); title('Source 4'); view([0 90]);
% subplot(2, 3, 5); showsurface(sa.cortex10K, [], out_cortex5); title('Source 5'); view([0 90]);
% subplot(2, 3, 6); showsurface(sa.cortex10K, [], out_cortex6); title('Source 6'); view([0 90]);
% 
% 
% ___________________________________________________VIEW WITH 4 DIPOLES____________________________________________________________


figure;
showsurface(sa.cortex10K);
hold on;
title('True Source Dipoles', 'FontSize', 19, 'FontWeight', 'bold');


true_voxel_indices = [dipole_idx1, dipole_idx2];  
%,dipole_idx3,dipole_idx4
%,dipole_idx5,dipole_idx6

dipole_dirs = [dir1, dir2];  
%,dir3,dir4
%,dir5,dir6


arrow_length = 2;
marker_radius = 0.5;

for i = 1:n_sources
    idx = true_voxel_indices(i);
    origin = sa.grid_medium(idx, :);
    dir = dipole_dirs(:, i);
    endpoint = origin + arrow_length * dir';

   
    plot3([origin(1), endpoint(1)], ...
          [origin(2), endpoint(2)], ...
          [origin(3), endpoint(3)], ...
          'Color', 'r', 'LineWidth', 2);

    
    [sx, sy, sz] = sphere(10);
    surf(marker_radius*sx + origin(1), ...
         marker_radius*sy + origin(2), ...
         marker_radius*sz + origin(3), ...
         'FaceColor', 'r', 'EdgeColor', 'none');
end

lighting gouraud;
camlight headlight;
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');


%%  _________________________ADD NOISE_____________________________________________


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Source Space Noise
% 
% % source-level noise across all voxels
% n_voxels = size(sa.V_medium, 2);
% T        = size(sources_sim, 1);   
% source_noise = randn(T, n_voxels, 3);  % uncorrelated in source space
% 
% % project to sensor space
% sensor_noise = zeros(T, nchan);
% 
%     for d = 1:3
%         Ld = squeeze(sa.V_medium(:, :, d));  % [nchan x n_voxels] = L(c,v,d)
%         Nd = squeeze(source_noise(:, :, d));  % [T x n_voxels]    = n(t,v,d)
%         sensor_noise = sensor_noise +  Nd* Ld';  % [T x nchan], sums over v
%    
%      end
% 
% 
% alpha_noise = data2filtdata(sensor_noise, sr, [8,12], sr+1);
% beta_noise  = data2filtdata(sensor_noise, sr, [20,30], sr+1);
% 
% 
% background = 0.5 * alpha_noise + 0.5 * beta_noise;
% %  ---------------------------SNR------------------------------------------
% 
%  % --- scale noise to target SNR ---
%         signal_power = mean(sources_sim(:).^2);
%         noise_power = mean(background(:).^2);
%         snr_linear = 10^(snr_db / 10);
%         scale_factor = sqrt(signal_power / (snr_linear * noise_power));
%         background_scaled = background * scale_factor;
% 
%         
%        data = sources_sim + background_scaled;
% %-------------------------------------------------------------------------------
% 
      data= sources_sim;


%% _____________________________________START_____________________________________________________________________________________________________________________
 

[cs, csnr, nave] = data2bs_univar(data, segleng, segshift, epleng, maxfreqbins, []); 

xx = abs(cs ./ csnr); 

% frequency bins for alpha band
alpha_band = [9, 14];    

% initialize
max_bicoherence = 0;
max_channel = 0;
max_f1_idx = 0;
max_f2_idx = 0;

% find max bicoherence
for channel = 1:nchan;
  
    f1_idx_start = alpha_band(1);
    f1_idx_end = alpha_band(2);
    
    
    for f1_idx = f1_idx_start:f1_idx_end;
        for f2_idx = f1_idx_start:f1_idx_end;
           
            bicoherence = xx(channel, f1_idx, f2_idx);
            
           
            if bicoherence > max_bicoherence;
                max_bicoherence = bicoherence;
                max_channel = channel;
                max_f1_idx = f1_idx;
                max_f2_idx = f2_idx;
            end
        end
    end
end

% display
fprintf('Maximum Biocoherence: %.4f\n', max_bicoherence);
fprintf('Channel index with maximum Biocoherence: %d\n', max_channel);
fprintf('Frequency pair with maximum Biocoherence in alpha band: (f1 = %d, f2 = %d)\n', max_f1_idx, max_f2_idx);

figure;
imagesc(squeeze(xx(max_channel, :, :))); % xx for channel ...
cb = colorbar;
title('Maximum Bicoherence', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Frequency 1', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Frequency 2', 'FontSize', 16, 'FontWeight', 'bold');

set(gca, 'FontSize', 18);      % axis tick labels size
set(cb, 'FontSize', 18);       % colourbar tick labels size

hold on;
plot(max_f1_idx, max_f2_idx, 'ro', 'MarkerSize', 10, 'LineWidth', 2);


%% ____________________________________________ Calculate Cross Bispectrum _________________________________________________


    dominant_freq = max_f2_idx; 
    freqpairs(1,:) = [dominant_freq dominant_freq];
    para = []; 
    [bs, nave] = data2bs_event(data, segleng, segshift, epleng, freqpairs, para); %bs: 61x61x61 complex
    
    

%% _____________________________________Low Rank Model________________________________________

n = 2;  % example model order
para = [];
[a, d, err, err_all, bsmodel] = bsfit(bs, n, para);


%% _____________________________________________ DEMIX SOURCES (MOCA) ___________________________________________________________

lead = sa.V_medium; 
A = mkfilt_eloreta(lead, 0.5); 

[nchan, nvoxel, ndum] = size(A);
F = zeros(nvoxel, ndum, n);

for i = 1:n
    for k = 1:ndum
        F(:, k, i) = A(:,:,k)' * a(:, i);
    end
end

[Fout, wall] = moca_ncomp(F);

disp(['F size: ', mat2str(size(F))]);
disp(['Fout size: ', mat2str(size(Fout))]);

%% _____________________________________________ VISUALIZE ON CORTEX ___________________________________________________________

figure;
num_rows = ceil(sqrt(n)); 
num_cols = ceil(n / num_rows); 

titleFont = 30;      % subplot title size
cbFont    = 20;      % colourbar tick size

for i = 1:n
    Floc = sqrt(sum(Fout(:, :, i).^2, 2));  % magnitude at each voxel
    out_cortex = spatfiltergauss(Floc, sa.grid_medium, 2, sa.cortex10K.vc);

    subplot(num_cols, num_rows, i);
    showsurface(sa.cortex10K, [], out_cortex);
    
    
    % ----- Title -----
    title(sprintf('Demixed Source %d', i), ...
          'FontSize', titleFont, 'FontWeight', 'bold');

    % ----- Colourbar -----
    cb = colorbar;
    set(cb, 'FontSize', cbFont);


    
    view([0 90]);
end



% Loop through each demixed source
for i = 1:n
    Floc = sqrt(sum(Fout(:, :, i).^2, 2));  % magnitude across orientations
    [peak_val, peak_voxel] = max(Floc);
    fprintf('Source %d peak voxel index:', i, peak_voxel);
end

true_voxels = [dipole_idx1, dipole_idx2];  % S1 and S2

for i = 1:n
    Floc = sqrt(sum(Fout(:,:,i).^2, 2));
    [~, peak_voxel] = max(Floc);
    dists = sqrt(sum((sa.grid_medium(peak_voxel,:) - sa.grid_medium(true_voxels,:)).^2, 2));
    
    [min_dist, closest_idx] = min(dists);
    fprintf('Demixed Source %d is closest to true source %d (voxel %d), distance = %.2f mm\n', ...
            i, closest_idx, true_voxels(closest_idx), min_dist);
end

out_demixed{i} = spatfiltergauss(sqrt(sum(Fout(:,:,i).^2,2)), sa.grid_medium, 2, sa.cortex10K.vc);

for i = 1:n
    out_demixed{i} = spatfiltergauss(sqrt(sum(Fout(:,:,i).^2,2)), sa.grid_medium, 2, sa.cortex10K.vc);
end


% true source maps
true_maps = {out_cortex1, out_cortex2};

% demixed maps
for i = 1:n
    out_demixed{i} = spatfiltergauss(sqrt(sum(Fout(:,:,i).^2,2)), sa.grid_medium, 2, sa.cortex10K.vc);
end

% corr matrix [n_demixed x n_true]
corr_matrix = zeros(n, length(true_maps));
for i = 1:n
    for j = 1:length(true_maps)
        corr_matrix(i, j) = corr(out_demixed{i}(:), true_maps{j}(:));
    end
end


% matching one demixed source to each true map
used = []; 
for j = 1:length(true_maps)
    [sorted_vals, sorted_idx] = sort(corr_matrix(:, j), 'descend');
    for k = 1:length(sorted_idx)
        d_idx = sorted_idx(k);
        if ~ismember(d_idx, used)
            used(end+1) = d_idx;
            fprintf('true source %d best matches demixed source %d (corr = %.3f)\n', ...
                j, d_idx, corr_matrix(d_idx, j));
            break;
        end
    end 
end

%%

figure('Color', 'w', 'Position', [100, 100, 1000, 800]);

titleFont = 24;     % bigger titles
cbFont    = 18;     % slightly bigger colourbar ticks


% -------------------------------
% top row: True Source Distributions
subplot(2, 2, 1);
showsurface(sa.cortex10K, [], out_cortex1);
title('True Source 1', 'FontSize', titleFont, 'FontWeight', 'bold');
view([0 90]);

subplot(2, 2, 2);
showsurface(sa.cortex10K, [], out_cortex2);
title('True Source 2', 'FontSize', titleFont, 'FontWeight', 'bold');
view([0 90]);

% -------------------------------
% bottom row: Demixed Sources
for i = 1:2
    Floc = sqrt(sum(Fout(:, :, i).^2, 2));
    out_cortex = spatfiltergauss(Floc, sa.grid_medium, 2, sa.cortex10K.vc);

    subplot(2, 2, 2 + i);
    showsurface(sa.cortex10K, [], out_cortex);
    title(sprintf('Demixed Source %d', i), 'FontSize', titleFont, 'FontWeight', 'bold');
    view([0 90]);
end



