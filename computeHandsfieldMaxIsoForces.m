%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        computeHandsfieldMaxIsoForces                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Riad Akhundov <Riad.Akhundov@uon.edu.au>
%
%Fractional muscle volumes (% of cumulative lower limb volumes)
%From Handsfield et al., 2014 JBiomech
%Standard value for rho from Rajagopal et al., 2015 = 60

function modelOut = computeHandsfieldMaxIsoForces(modelIn, acquisitionInfo, rho, modelOutName)
%% Assign Handsfield Muscle Volume Fractions (from Appendix A. Supplementary materials)
%Works with any model if you change the fieldnames and adjust the compartmented muscles (e.g. glmax1)
%(this version coded for the Rajagopal_2015 model)

make_plots = true; %Plot comparisson of generic vs adjusted maxIsoForces
% make_plots = false;

figure_folder = modelOutName(1:find(modelOutName == '\', 1, 'last'));

muscleVolumeFractions.addbrev = 0.0147;
muscleVolumeFractions.addlong = 0.0226;
muscleVolumeFractions.addmag = 0.0786;
muscleVolumeFractions.bflh = 0.0292;
muscleVolumeFractions.bfsh = 0.0140;
muscleVolumeFractions.edl = 0.0144; 
muscleVolumeFractions.ehl = 0.0144;  
muscleVolumeFractions.fdl = 0.0043;
muscleVolumeFractions.fhl = 0.0109;
muscleVolumeFractions.gaslat = 0.0211; 
muscleVolumeFractions.gasmed = 0.0362;
muscleVolumeFractions.glmax = 0.1193;
muscleVolumeFractions.glmed = 0.0454; 
muscleVolumeFractions.glmin = 0.0147;
muscleVolumeFractions.grac = 0.0146;
muscleVolumeFractions.iliacus = 0.0248;
muscleVolumeFractions.perbrev = 0.0183; 
muscleVolumeFractions.perlong = 0.0183;
muscleVolumeFractions.piri = 0.0061;
muscleVolumeFractions.psoas = 0.0380;
muscleVolumeFractions.recfem = 0.0379; 
muscleVolumeFractions.sart = 0.0229; 
muscleVolumeFractions.semimem = 0.0346; 
muscleVolumeFractions.semiten = 0.0260;
muscleVolumeFractions.soleus = 0.0621;
muscleVolumeFractions.tfl = 0.0089; 
muscleVolumeFractions.tibant = 0.0191; 
muscleVolumeFractions.tibpost = 0.0149;
muscleVolumeFractions.vasint = 0.0384;
muscleVolumeFractions.vaslat = 0.1166; 
muscleVolumeFractions.vasmed = 0.0606; 
 

%% Theoretical lower-limb muscle volume (unilateral) from Handsfield equations
height = acquisitionInfo.Subject.Height; %(m)
mass = acquisitionInfo.Subject.Weight; %(kg)

vTheory = (47*mass*height) + 1285; %mass=weight in kg | height in m


%% Import Model
import org.opensim.modeling.*
model = Model(modelIn);
model.initSystem;

muscles = model.getMuscles();
nMuscles = muscles.getSize();

muscleNames = cell(nMuscles, 1);
muscleForce = zeros(nMuscles,1);
muscleOptFiberLength = zeros(nMuscles,1);
maxIsoForce = zeros(nMuscles,1);


%% Import Muscle With Compartments
%Assumes right and left side maxIsoForces of model are identical (other muscle parameters can vary)
% addmag
addmagDist = muscles.get('addmagDist_r'); addmagDist = addmagDist.getMaxIsometricForce;
addmagIsch = muscles.get('addmagIsch_r'); addmagIsch = addmagIsch.getMaxIsometricForce;
addmagMid = muscles.get('addmagMid_r'); addmagMid = addmagMid.getMaxIsometricForce;
addmagProx = muscles.get('addmagProx_r'); addmagProx = addmagProx.getMaxIsometricForce;
%glmax
glmax1 = muscles.get('glmax1_r'); glmax1 = glmax1.getMaxIsometricForce;
glmax2 = muscles.get('glmax2_r'); glmax2 = glmax2.getMaxIsometricForce;
glmax3 = muscles.get('glmax3_r'); glmax3 = glmax3.getMaxIsometricForce;
%glmed
glmed1 = muscles.get('glmed1_r'); glmed1 = glmed1.getMaxIsometricForce;
glmed2 = muscles.get('glmed2_r'); glmed2 = glmed2.getMaxIsometricForce;
glmed3 = muscles.get('glmed3_r'); glmed3 = glmed3.getMaxIsometricForce;
%glmin
glmin1 = muscles.get('glmin1_r'); glmin1 = glmin1.getMaxIsometricForce;
glmin2 = muscles.get('glmin2_r'); glmin2 = glmin2.getMaxIsometricForce;
glmin3 = muscles.get('glmin3_r'); glmin3 = glmin3.getMaxIsometricForce;


%% Adjust Model
for i = 0:nMuscles-1   
    currentMuscle = muscles.get(i);
    muscleNames{i+1} = char(currentMuscle.getName());
    muscleNames_figure{i+1} = muscleNames{i+1}(1:end-2);
    muscleForce(i+1) = currentMuscle.getMaxIsometricForce();
    muscleOptFiberLength(i+1) = currentMuscle.getOptimalFiberLength()*100; % in cm
    %Muscle compartment calculations
    if contains(muscleNames{i+1},'addmagDist') %Assumes right and left side maxIsoForces of model are identical (other muscles parameters can vary)
        addmagDist_fraction = addmagDist/(addmagDist+addmagIsch+addmagMid+addmagProx);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.addmag)/muscleOptFiberLength(i+1))*addmagDist_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'addmagIsch')
        addmagIsch_fraction = addmagIsch/(addmagDist+addmagIsch+addmagMid+addmagProx);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.addmag)/muscleOptFiberLength(i+1))*addmagIsch_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'addmagMid')
        addmagMid_fraction = addmagMid/(addmagDist+addmagIsch+addmagMid+addmagProx);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.addmag)/muscleOptFiberLength(i+1))*addmagMid_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
            
    elseif contains(muscleNames{i+1},'addmagProx')      
        addmagProx_fraction = addmagProx/(addmagDist+addmagIsch+addmagMid+addmagProx);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.addmag)/muscleOptFiberLength(i+1))*addmagProx_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmax1')
        glmax1_fraction = glmax1/(glmax1+glmax2+glmax3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmax)/muscleOptFiberLength(i+1))*glmax1_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmax2')
        glmax2_fraction = glmax2/(glmax1+glmax2+glmax3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmax)/muscleOptFiberLength(i+1))*glmax2_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmax3')
        glmax3_fraction = glmax3/(glmax1+glmax2+glmax3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmax)/muscleOptFiberLength(i+1))*glmax3_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmed1')
        glmed1_fraction = glmed1/(glmed1+glmed2+glmed3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmed)/muscleOptFiberLength(i+1))*glmed1_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmed2')
        glmed2_fraction = glmed2/(glmed1+glmed2+glmed3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmed)/muscleOptFiberLength(i+1))*glmed2_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmed3')
        glmed3_fraction = glmed3/(glmed1+glmed2+glmed3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmed)/muscleOptFiberLength(i+1))*glmed3_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmin1')
        glmin1_fraction = glmin1/(glmin1+glmin2+glmin3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmin)/muscleOptFiberLength(i+1))*glmin1_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmin2')
        glmin2_fraction = glmin2/(glmin1+glmin2+glmin3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmin)/muscleOptFiberLength(i+1))*glmin2_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
    elseif contains(muscleNames{i+1},'glmin3')
        glmin3_fraction = glmin3/(glmin1+glmin2+glmin3);
        maxIsoForce(i+1) = ((rho*vTheory*muscleVolumeFractions.glmin)/muscleOptFiberLength(i+1))*glmin3_fraction;
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));
        
        
    %All other muscles    
    elseif isfield(muscleVolumeFractions, muscleNames{i+1}(1:end-2))        
        maxIsoForce(i+1) = (rho*vTheory*muscleVolumeFractions.(muscleNames{i+1}(1:end-2)))/muscleOptFiberLength(i+1);
        
        currentMuscle.setMaxIsometricForce(maxIsoForce(i+1));        
    end
end


%% Print model with new muscle strengths
model.setName([acquisitionInfo.Subject.Code, '_strengthAdjusted']);
modelOut = modelOutName;
model.print(modelOut);
disp(['The new model has been saved at ' modelOut]);


%% Comparisson Figures
if make_plots == 1
    %R 
    muscleNames_figure = muscleNames_figure';
    figure_forces_r = muscleNames_figure(1:40);
    figure_forces_r(:,2) = num2cell(muscleForce(1:40));
    figure_forces_r(:,3) = num2cell(maxIsoForce(1:40));

    figure_forces_r = figure_forces_r(all(cell2mat(figure_forces_r(:,3)) ~= 0,3),:);

    h1 = figure('WindowState', 'maximized');
    figure_forces_r_plot = [cell2mat(figure_forces_r(:,2)), cell2mat(figure_forces_r(:,3))];
    c = categorical(figure_forces_r(:,1));
    bar(c,figure_forces_r_plot)
    legend({'Original', 'Handsfield'}, 'Location', 'northeastoutside');
    ylabel('Force (N)');
    title('Maximum Isometric Muscle Forces R');

    %L
    figure_forces_l = muscleNames_figure(41:80);
    figure_forces_l(:,2) = num2cell(muscleForce(41:80));
    figure_forces_l(:,3) = num2cell(maxIsoForce(41:80));

    figure_forces_l = figure_forces_l(all(cell2mat(figure_forces_l(:,3)) ~= 0,3),:);

    h2 = figure('WindowState', 'maximized');
    figure_forces_l_plot = [cell2mat(figure_forces_l(:,2)), cell2mat(figure_forces_l(:,3))];
    c = categorical(figure_forces_l(:,1));
    bar(c,figure_forces_l_plot)
    legend({'Original', 'Handsfield'}, 'Location', 'northeastoutside');
    ylabel('Force (N)');
    title('Maximum Isometric Muscle Forces L');

    savefig(h1, [figure_folder, 'maxIsoForce_R.fig']);
    savefig(h2, [figure_folder, 'maxIsoForce_L.fig']);
    save([figure_folder, 'genericVsHandsfield_forces.mat'], 'figure_forces_r', 'figure_forces_l');
    close all;
end
