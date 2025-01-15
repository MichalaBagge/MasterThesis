function [eog_l, eog_r, hypnogram_ext] = load_polysomnogram(dir_data,subject_ID)

dir_subject = fullfile(dir_data,subject_ID);
cd(dir_subject)

% Load eog (left,right) and hypnogram
load('eogl-m2.mat'); load('eogr-m2.mat'); load('vec_hypnogram.mat');
eog_l = eoglm2; eog_r = eogrm2;
clear eogrm2 eoglm2

% Create extended hypnogram
hypnogram_ext = reshape(transpose(repmat(hypnogram, 1, 256*30)), [], 1);

end