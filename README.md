# computeHandsfieldMaxIsoForces

computeHandsfieldMaxIsoForces.m - 
This MATLAB function is used to adjust OpenSim model maximum isometric muscle forces in accordence with the Handsfield equation (Handsfield et al., 2014 JBiomech).

It works with any model if you change the fieldnames and adjust the compartmented muscles (e.g. glmax1), but this specific this version is set up for the 
Rajagopal_2015 model (Rajagopal et al., 2015).

------------------------------------------------------------------------------------------------------------------------------------------------------------------

Expected function inputs are:


1) modelIn - Path to input OpenSim model file (.osim) - (use '...\MyModel.osim' instead of "...\MyModel.osim" to avoid errors)

2) acquisitionInfo - A struct containing subject specific information

    acquisitionInfo.Subject.Height = height; %(m)
    acquisitionInfo.Subject.Weight = massOriginal; %(kg)
    acquisitionInfo.Subject.Code = subject; %Subject ID (e.g., 'BA03')
    
3) rho - Muscle strain value - Standard value for rho from Rajagopal et al., 2015 = 60 (just use rho = 60 unless newest literature says otherwise)

4) modelOutName - Path to output OpenSim model file (.osim) - (again, use '...\MyModel.osim' instead of "...\MyModel.osim" to avoid errors)

------------------------------------------------------------------------------------------------------------------------------------------------------------------

Call this function (after adding it and OpenSim 3.3 to the MATLAB path & PC path) with:

computeHandsfieldMaxIsoForces(osimModel_targ_filepath, acquisitionInfo, rho, osimModel_adjusted_filepath);


For more information on the CEINMS pipeline (including a CEINMS step by step guide) see bit.ly/CEINMS_Pipeline_Redone_v1_2
