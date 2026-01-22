%% Real-data example pipeline (resting-state EEG)
% run_realdata_example.m
%
% Note: Raw EEG data are not included in this repository.
% This script demonstrates the processing steps used in the manuscript.
%
% Dependencies:
%  - METH toolbox: https://www.uke.de/english/departments-institutes/institutes/neurophysiology-and-pathophysiology/research/research-groups/index.html
%  - bsfit toolbox: https://github.com/guidonolte/bsfit
%
% Setup:
%  1) Copy config/config_paths_example.m -> config/config_paths.m
%  2) Edit paths in config/config_paths.m
%  3) Run this script


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
load(fullfile(cfg.path_meth, 'templates', 'mri.mat'));

 %Parameters:
 sr = 256; 
 segleng = 256; 
 segshift = 128;
 epleng = 512;
 nchan = 61;
 maxfreqbins = 51;
    
%% ______________________________________LOAD REAL DATA__________________________________________________________________________________


assert(exist(cfg.path_allnames,'file')==2, 'allnames file not found. Check cfg.path_allnames.');
load(cfg.path_allnames);  % should load allnames_K, allnames_P, etc.

nK=length(allnames_K); 
nP=length(allnames_P);
nR=length(allnames_R);


for isub = 12   %1:nK % example subject index
    igroup = 1; % 1 = K group, 0/2 = other

    if igroup == 1
      dataFile = allnames_K{isub};
    elseif igroup == 2
      dataFile = allnames_P{isub};
    else
      dataFile = allnames_R{isub};
    end
    
	
	if ~isfile(dataFile)
    	  dataFile = fullfile(cfg.path_realdata, dataFile);
	end
	assert(isfile(dataFile), 'Data file not found: %s', dataFile);

	data = readmatrix(dataFile);
	if any(isnan(data(:)))
    	  warning('readmatrix returned NaNs; falling back to dlmread.');
          data = dlmread(dataFile);
	end

%% _____________________________________START___________________________________________________________________________________________

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


% Figures:

   figure;
   imagesc(squeeze(xx(max_channel, :, :))); 
   cb = colorbar;
   title('Maximum Bicoherence', 'FontSize', 16, 'FontWeight', 'bold');
   xlabel('Frequency 1', 'FontSize', 16);
   ylabel('Frequency 2', 'FontSize', 16);

   hold on;
   plot(max_f1_idx, max_f2_idx, 'ro', 'MarkerSize', 10, 'LineWidth', 2);



    % Figure 1
    figure1 = figure;
    showtfinhead(xx, sa.locs_2D);

   
    % Figure 2
    figure2 = figure;
    imagesc(squeeze(xx(max_channel,:,:)));
    colorbar;
    
    title(sprintf('Dom Freq is %d and Dom Ch. is %d', max_f1_idx,max_channel));
  


    % Figure 3
    figure3 = figure;
    set(figure3, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 
      for i = 1:61;
        subplot(8, 8, i);
        imagesc(squeeze(xx(i,:,:)));
        colorbar;
      end

   
% Define frequency pairs

    freqpairs(1,:)=[max_f1_idx max_f2_idx]


% Calculate Cross Bispectrum 
    dominant_freq = max_f1_idx; 
    freqpairs(1,:) = [dominant_freq dominant_freq];
    para = []; 
    [bs, nave] = data2bs_event(data, segleng, segshift, epleng, freqpairs, para); 
    
   
% Low-Rank Model
    n = 2;
    para = [];
    [a, d, err, err_all, bsmodel] = bsfit(bs, n, para);
    

  
%% ________________________________________  DEMIX SOURCES ________________________________________________________________

lead=sa.V_medium;

A=mkfilt_eloreta(lead,0.5); 

% Make an inverse with MOCA
% First, find mixed sources

[nchan nvoxel ndum]=size(A);
F=zeros(nvoxel,ndum,n); % intitialze distributions for n sources
 
for i=1:n; for k=1:ndum;
F(:,k,i)=A(:,:,k)'*a(:,i);
end;end
 
figure7=figure;

num_rows = ceil(sqrt(n)); 
num_cols = ceil(n / num_rows); 

% Loop through each source

for i = 1:n
    
    Fm_un= sqrt(sum(F(:, :, i).^2, 2)); 
   
    d = 2; 
    out_cortex_un = spatfiltergauss(Fm_un, sa.grid_medium, d, sa.cortex10K.vc);
    
    
    subplot(num_cols, num_rows, i);
    para = [];  
    % para.mycolormap = 'hot'; 
    showsurface(sa.cortex10K, para, out_cortex_un);
    title(sprintf('Mixed Source %d', i));
    view([0 90]);
end


% Demix sources

[Fout,wall]=moca_ncomp(F);


disp(size(F)); % should be [nvoxel, ndum, n]
disp(size(Fout));
 
% check the strength of unmixed sources
for i = 1:n
    disp(['Source ', num2str(i), ' strength: ', num2str(norm(Fout(:, :, i)))]);
end


% View on Cortex

figure8=figure;
%set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); 

num_rows = ceil(sqrt(n)); 
num_cols = ceil(n / num_rows); 

% loop through each unmixed source
for i = 1:n
    
    Floc = sqrt(sum(Fout(:, :, i).^2, 2)); % Magnitude of source strength
   
    d = 2; 
    out_cortex = spatfiltergauss(Floc, sa.grid_medium, d, sa.cortex10K.vc);

    subplot(num_cols, num_rows, i);
    para = [];  
    % para.mycolormap = 'hot'; 
    showsurface(sa.cortex10K, para, out_cortex);
    title(sprintf('Demixed Source %d', i));
    view([0 90]);
end


Fout_sum = squeeze(sum(Fout, 2));
    
demixed_signals = Fout_sum * a';

end






