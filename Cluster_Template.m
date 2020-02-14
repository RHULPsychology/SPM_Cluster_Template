
function Cluster_Template_cars_rsa_v8_2019_rstab(subject);

%Jan 2020: Inserts rsa toolbox searchlight in place of cosmo searchlight,
%but with all same preprocessing.

%2019: I'm going to attempt to repeat rsa analysis using correct identities
%using the previous rsa first level models (which have one regressior per
%trial in fMRI trial order) and cosmo processing

%v5 attempts to fix the bugs in previous implementations. Works on new first
%level models where betas are properly sorted in trial order and trial data
%is saved in organised way. Also it explicitly and correctly corrects for
%mismatches between the expression and condition codes used in begavioural
%versus fMRI experiments.

%v4: adds parfor loops for cluster
%v3: adds expression-speciifc dm

%which_matrix == 3 means behavioral dissimilarities between all pairwise videos, 120*120
%which_matrix == 5 means 120*120 theoretical matrix
which_matrix = 3;

%add cosmo and spm
addpath(genpath('/MRIWork/MRIWork03/nf/spm12'));
addpath(genpath('/MRIWork/MRIWork03/nf/forida_begum/matlab_files'));    %will load cosmo too


%basic file management
raw_data_path = '/MRIWork/MRIWork03/nf/forida_begum/Raw_Data';
data_folders = dir([raw_data_path filesep '*2018-0*']);
num_subs = size(data_folders, 1);
first_level_folder = 'first_level04redone';
subject=str2num(subject);

SubjectDir = [raw_data_path filesep data_folders(subject).name filesep first_level_folder];
disp(SubjectDir)
hFunc(raw_data_path, subject, SubjectDir, which_matrix);

%for subject=[2:31]
%    
%    SubjectDir = [raw_data_path filesep data_folders(subject).name filesep first_level_folder];
%    Func(raw_data_path, subject, SubjectDir, which_matrix);
%
%end

disp('audi5000');


function hFunc(raw_data_path, subject, SubjectDir, which_matrix);

mask_filename = '/MRIWork/MRIWork03/nf/spm12/spm12/canonical/wc1single_subj_T1.nii';

%fMRI trial data: col1: onset, col2: cond, col3: identity, col4: expression,col5: run
load([SubjectDir filesep sprintf('firstlevel04_redone3_data_sub%02d.mat',subject)],'trial_data');
load('/MRIWork/MRIWork03/nf/forida_begum/matlab_files/rsa_matrices.mat');

% %make the important corrections
% %in behavior: 1=am 2=amac 3=ac 4=c
% %in fmri: 1=c 2=am 4=amac 3=ac
% %make like behavior:
% clear temp;
% temp = trial_data;
% temp(find(trial_data(:,2)==1),2) = 3;
% temp(find(trial_data(:,2)==2),2) = 2;
% temp(find(trial_data(:,2)==3),2) = 1;
% temp(find(trial_data(:,2)==4),2) = 4;
% trial_data = temp;

%make code table
it = 1;
for cond = [3 2 1 4];
    for exp = 1:5;
        for id=[1 3 2 5 4 6];
            code_table(it,1) = cond;
            code_table(it,2) = exp;
            code_table(it,3) = id;
            it = it+1;
        end;
    end;
end;

%recode fmri data into indices into behavioural codes to use as rsa dm 
%labels (i.e., "targets" in cosmoMVPA)
targets = [];
for trial = 1:size(trial_data,1);
    match = [];
    for code = 1:size(code_table,1);
        if sum( code_table(code,:) == trial_data(trial,[2 4 3])) == 3;
            match = code;
        end;
    end;
    targets(trial,1) = match;
end;

%I'm only doing big matrices now
if which_matrix == 1;   %expression matrix
elseif which_matrix == 2;   %simple condition average matrix;   
elseif which_matrix == 3;   %120*120 matrix of all behavioural dissimilarities
    
    %zero out diagonal so cosmo won't choke
    for i=1:size(big_behav_ave,1);
        big_behav_ave(i,i) = 0;
    end;
    dm = squareform(squareform(big_behav_ave));
    %replace NaNs with average values to avoid stats problems later
    dm(isnan(dm)) = nanmean(nanmean(dm));
    
    %     %where there's missing data in behav matrix, blank it out in control too
    %     big_behav_control(isnan(dm))=NaN;
%     regress_dm = big_behav_control;
    
    if exist('regress_dm');
        stem = 'cosmo_2019_behavDMwithControl';
    else
        stem = 'cosmo_2019_behavDMnoControl';
    end;
    
elseif which_matrix == 4;   %4*4 matrix of behavioural dimmiliaries averaged for each condition
elseif which_matrix == 5;   %120*120 big theoretical matrix
    
    dm = big_matrix;   
    stem = 'cosmo_2019_TheorDMnoControl';
 
end;    %which matrix?


%should be all files and session directory
datafiles = spm_select('FPList', SubjectDir, '^beta.*\.nii$');
% 
% %%%%%%%%%%%%%%%%%This bit requires batch cluster or will out of memory%%%%%%%%%%%%%%%%%%
% 
% disp('loading datasets');
% save('/MRIWork/MRIWork03/nf/forida_begum/loading');
% ds_it = 1;
% %only first 720 are trials, remainders are confound regressors
% for beta_file=1:size(trial_data,1);
%     
%     cosmo_warning('off');
% %     individual_datasets{ds_it} = cosmo_remove_useless_data(cosmo_fmri_dataset( spm_str_manip(datafiles(beta_file,:),'d'), 'mask',mask_filename));
%     individual_datasets{ds_it} = cosmo_fmri_dataset( spm_str_manip(datafiles(beta_file,:),'d'));
%     ds_it = ds_it + 1;
%     
% end;    %load betafiles
% 
% save('/MRIWork/MRIWork03/nf/forida_begum/about_to_stack');
% 
% ds_stacked = cosmo_stack(individual_datasets);
% clear individual_datasets;
% 
% save('/MRIWork/MRIWork03/nf/forida_begum/stacked');
% 
% ds_stacked.sa.targets = targets; %observations
% ds_stacked.sa.chunks = trial_data(:,5);
% 
% %get condition averages
% ds_aves = cosmo_fx(ds_stacked, @(x)mean(x,1), 'targets', 1)
% save('/MRIWork/MRIWork03/nf/forida_begum/averaged');
% %%%%%%%%%%%%%%%%%This bit requires batch cluster%%%%%%%%%%%%%%%%%%

trunc_dm = 120;

load('/MRIWork/MRIWork03/nf/forida_begum/averaged','ds_aves');  %use this only for testing without the cluster scheduling system - precomputed average so avoids out of memory error when loading datasets
% temp_data = ds_aves.samples(:,1:125)';
temp_data = ds_aves.samples(1:trunc_dm,1:end)';

%put together minimaloptions needed to make/visualize test dm
userOptions.analysisName = 'rsatb_searchlight01_Jan2020';
userOptions.rootPath = '/MRIWork/MRIWork03/nf/forida_begum/Raw_Data/decoding_results/rsa_toolbox/rsatb_searchlight01_Jan2020';    %I think it only saves output here
userOptions.distance = 'correlation';
userOptions.RoIColor = [0 0 0];    %let it be black by default,
userOptions.saveFigurePDF = 1;
userOptions.saveFigurePS = 0;
userOptions.saveFigurefig = 1;
userOptions.saveFiguresPDF = 1;
userOptions.saveFiguresPS = 0;
userOptions.saveFiguresfig = 1;
userOptions.displayFigures = 1;
mcolors = cool(2);
userOptions.ModelColor = mcolors(1,:);
local_options_fig = struct('filename','ROI_dms','figureNumber',1);

%make test DM struct
dms.dm = dm(1:trunc_dm,1:trunc_dm);
RDMs_model = rsa.constructModelRDMs(dms,userOptions);
rsa.figureRDMs(RDMs_model, userOptions, local_options_fig);

%load and reshape mask
mask_ds = cosmo_fmri_dataset( spm_str_manip(mask_filename,'d'));
% mask_rs = reshape(mask_ds.samples(1:125),5,5,5);
mask_rs = reshape(mask_ds.samples,79,95,79);

%more userOptions
userOptions.voxelSize = [2 2 2];
userOptions.searchlightRadius = 8;

%local options
searchlightOptions.monitor = 0;
searchlightOptions.fisher = 1;
searchlightOptions.nSessions = 1;
searchlightOptions.nConditions = trunc_dm;

save('/MRIWork/MRIWork03/nf/forida_begum/about_to_SL');
[rs, ps, ns, searchlightRDMs] = rsa.fmri.searchlightMapping_fMRI(temp_data, RDMs_model, mask_rs, userOptions, searchlightOptions);
save('/MRIWork/MRIWork03/nf/forida_begum/finished_SL');
save('/MRIWork/MRIWork03/nf/forida_begum/SL_out_data','rs','ps','ns','searchlightRDMs');




% [rs, ps, ns, searchlightRDMs.(subject)] = rsa.fmri.searchlightMapping_fMRI(singleSubjectVols, models, mask, userOptions, searchlightOptions);
% 
% %pass dataset, neighborhoods, measure and its arguments into searchlight
% neighborhood = cosmo_spherical_neighborhood( ds_aves, 'radius', 8);
% measure = @cosmo_target_dsm_corr_measure;   %pass the measure handle
% args.target_dsm = dm;  %arguments to pass into cosmo_target_dsm_corr_measure
% % args.glm_dsm = dm;
% args.center_data = true;
% % args.center_data = false;
% %     args.type = 'Kendall';
% if exist('regress_dm');
%     args.regress_dsm = regress_dm;
% end;
% %     args.metric = 'euclidean';
% results_map = cosmo_searchlight(ds_aves, neighborhood, measure, args);
% 
% results_map.samples = atanh(results_map.samples);
% 
% % save([raw_data_path filesep 'decoding_results' filesep sprintf('%s_sub%02d.mat',stem,subject)],results_map);
% 
% %write r to z transformed values
% cosmo_map2fmri(results_map, [raw_data_path filesep 'decoding_results' filesep sprintf('%s_sub%02d.nii',stem,subject)]);
% 
