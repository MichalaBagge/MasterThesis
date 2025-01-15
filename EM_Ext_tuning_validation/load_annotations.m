function [EM_annotations, REM_annotations, SEM_annotations,start_idx,stop_idx] = load_annotations(subject_ID)

% Set directory
dir_annotations = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Data\Eye movements con-glo Julie\EM manual scorings 2019';
cd(dir_annotations)

% Load light off / on targets 
cd(dir_annotations)
% EM stop annotation
stop = sprintf('%s_Fabio_StopEMscoring_targetvector.mat',subject_ID);
load(stop);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
stop_temp = find(targetVector_R == 1);
stop_idx = stop_temp(end);
clear targetVector_R stop
% EM start annotation
start = sprintf('%s_Fabio_StartEMscoring_targetvector.mat',subject_ID);
load(start);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
start_temp = find(targetVector_R == 1);
start_idx = start_temp(1);
clear targetVector_R start


% EM REM annotation
rem = sprintf('%s_Fabio_REM_targetvector.mat',subject_ID);
load(rem);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
REM_annotations= targetVector_R(start_idx:stop_idx);
clear targetVector_R rem
% EM SEM annotation
sem = sprintf('%s_Fabio_SEM_targetvector.mat',subject_ID);
load(sem);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
SEM_annotations = targetVector_R(start_idx:stop_idx);
clear targetVector_R sem
% EM annotation vector
EM_annotations = SEM_annotations + REM_annotations;
EM_annotations(EM_annotations==2) = 1;

end