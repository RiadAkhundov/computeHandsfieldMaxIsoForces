%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  batchComputeMuscleStrengthScaleFactors                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original script by David Saxby <d.saxby@griffith.edu.au>
% Edited by Riad Akhundov <riad.akhundov@uon.edu.au>

% Batch compute fractional muscle volumes from Handsfield et al 2014 JBiomech
% This vesrion works on Rajagopal_2015 (order of muscles and muscle names model dependant)

%%%Requirements: 
%1) MATLAB 2017b or newer (2018b or newer recommended)
%2) ...

%% Works with any model if you change the fieldnames and adjust the compartmented muscles (e.g. glmax1)
% rho = 60; %Normal muscle strain value (rho) from Rajagopal = 60. Set in main script

for s=1:nSubject
    modelFile = [subject{s} modelFileDescriptor '_opt_N' num2str(N_eval_set)];
    acquisitionInfo.Subject.Height = height(s); %(m)
    acquisitionInfo.Subject.Weight = massOriginal(s); %(kg)
    acquisitionInfo.Subject.Code = subject{s}; %Subject ID (e.g., 'BA03')
    
    for i=1:nSubjectSessions{s}
        osimModel_targ_filepath = [dirScaleModels{i,s} modelFile '.osim'];
        osimModel_adjusted_filepath = [osimModel_targ_filepath(1:end-5), '_strengthAdjusted.osim'];
        if i == 1 
            
            computeHandsfieldMaxIsoForces(osimModel_targ_filepath, acquisitionInfo, rho, osimModel_adjusted_filepath);
            disp(['%% ',modelFile,' created %%']);
            
        else            
            copyfile([dirScaleModels{1,s} modelFile '_strengthAdjusted.osim'], osimModel_adjusted_filepath);
            copyfile([dirScaleModels{1,s} 'maxIsoForce_R.fig'], [dirScaleModels{i,s} 'maxIsoForce_R.fig']);
            copyfile([dirScaleModels{1,s} 'maxIsoForce_L.fig'], [dirScaleModels{i,s} 'maxIsoForce_L.fig']);
            copyfile([dirScaleModels{1,s} 'genericVsHandsfield_forces.mat'], [dirScaleModels{i,s} 'genericVsHandsfield_forces.mat']);
            
            disp(['%% ',modelFile,' copied %%']);
        end
    end
end
