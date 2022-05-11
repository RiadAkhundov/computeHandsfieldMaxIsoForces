function modelOut = computeMuscleStrengthScaleFactors(modelIn, acquisitionInfo, rho, modelOutName)


%%
% fractional muscle volumes (% of cumulative lower limb volumes)
% from Handsfield et al 2014 JBiomech, only muscles in Rajaganopol model
% normal rho from Rajagopal = 60
%% Works only on Rajagopal_2015 (because the order of muscles is hardcoded)!

% hip muscles
glmax = 0.1193; addmag = 0.0786; glmed = 0.0454; psoas = 0.0380; iliacus = 0.0248;
sart = 0.0229; addlong = 0.0226; glmin = 0.0147; addbrev = 0.0147; grac = 0.0146;
tfl = 0.0089; piri = 0.0061;

% knee muscles
vaslat = 0.1166; vasmed = 0.0606; vasint = 0.0384; recfem = 0.0379; semimem = 0.0346; bflh = 0.0292;
semiten = 0.0260; bfsh = 0.0140;

% ankle muscles
soleus = 0.0621; gasmed = 0.0362; gaslat = 0.0211; tibant = 0.0191; peroneals = 0.0183; tibpost = 0.0149;
edl = 0.0144; ehl = 0.0144; fhl = 0.0109; fdl = 0.0043; perbrev = 0.0183; perlong = 0.0183;

%% List of Handsfield muscles that are also in Rajagonopol model
listOfHandsfieldMuscleNames = {'addbrev','addlong','addmag','bflh','bfsh','edl','ehl','fdl','fhl','gaslat',...
    'gasmed','glmax','glmed','glmin','grac','iliacus','perbrev','perlong','piri','psoas','recfem','sart','semimem',...
    'semiten','soleus','tfl','tibant','tibpost','vasint','vaslat','vasmed'};

%% List of muscles root names in the Rajaganopol model
listOfModelMusclesRootNames = {'addbrev','addlong','addmag','bflh','bfsh','edl','ehl','fdl','fhl','gaslat',...
    'gasmed','glmax','glmed','glmin','grac','iliacus','perbrev','perlong','piri','psoas','recfem','sart','semimem',...
    'semiten','soleus','tfl','tibant','tibpost','vasint','vaslat','vasmed'};

%% Import classes
import org.opensim.modeling.*
model = Model(modelIn);
model.initSystem;

%% Muscles in model
muscles = model.getMuscles();
nMuscles = muscles.getSize();

muscleNames = cell(nMuscles, 1);
muscleForce = zeros(nMuscles,1);
muscleOptFiberLength = zeros(nMuscles,1);
muscleVolume = zeros(nMuscles,1);

for i = 0:nMuscles-1
    
    currentMuscle = muscles.get(i);
    muscleNames{i+1} = char(currentMuscle.getName());
    muscleForce(i+1) = currentMuscle.getMaxIsometricForce();
    muscleOptFiberLength(i+1) = currentMuscle.getOptimalFiberLength()*100; % in cm
    muscleVolume(i+1) = (muscleForce(i+1)*muscleOptFiberLength(i+1))/rho;

end

height = acquisitionInfo.Subject.Height; %(m)
mass = acquisitionInfo.Subject.Weight; %(kg)

%% Theoretical lower-limb muscle volume (unilateral) from Handsfield equations
vTheory = (47*mass*height) + 1285; %mass=weight in kg | height in m

%% Loop through muscles and some default isometric forces per group, determine
% fractional isometric force

addmagIsoForce = zeros(4,2);
magIndex = 1;

glmaxIsoForce = zeros(3,2);
glmaxIndex = 1;

glmedIsoForce = zeros(3,2);
glmedIndex = 1;

glminIsoForce = zeros(3,2);
glminIndex = 1;

peronealsIsoForce = zeros(2,2);
peronealsIndex = 1;

extensorsIsoForce = zeros(2,2); %extensors edl/ehl were combined by Handsfield
extensorsIndex = 1;

for i = 1:length(listOfModelMusclesRootNames)
    
    lengthOfName = length(listOfModelMusclesRootNames{i});
    
    for j = 1:length(muscleNames)
        
        if strncmp(listOfModelMusclesRootNames{i}, muscleNames{j}, lengthOfName)
            
            if strncmp(muscleNames(j), 'addmag', 6)
                
                addmagIsoForce(magIndex,1) = muscleForce(j);
                magIndex = magIndex + 1;
                addmagIsoForce(magIndex,1) = muscleForce(j+1);
                magIndex = magIndex + 1;
                addmagIsoForce(magIndex,1) = muscleForce(j+2);
                magIndex = magIndex + 1;
                addmagIsoForce(magIndex,1) = muscleForce(j+3);
                
            elseif strncmp(muscleNames(j), 'glmax', 4)
                
                glmaxIsoForce(glmaxIndex,1) = muscleForce(j);
                glmaxIndex = glmaxIndex + 1;
                glmaxIsoForce(glmaxIndex,1) = muscleForce(j+1);
                glmaxIndex = glmaxIndex + 1;
                glmaxIsoForce(glmaxIndex,1) = muscleForce(j+2);
                
            elseif strncmp(muscleNames(j), 'glmed', 4)
                
                glmedIsoForce(glmedIndex,1) = muscleForce(j);
                glmedIndex = glmedIndex + 1;
                glmedIsoForce(glmedIndex,1) = muscleForce(j+1);
                glmedIndex = glmedIndex + 1;
                glmedIsoForce(glmedIndex,1) = muscleForce(j+2);
                
            elseif strncmp(muscleNames(j), 'glmin', 4)
                
                glminIsoForce(glminIndex,1) = muscleForce(j);
                glminIndex = glminIndex + 1;
                glminIsoForce(glminIndex,1) = muscleForce(j+1);
                glminIndex = glminIndex + 1;
                glminIsoForce(glminIndex,1) = muscleForce(j+2);
                
            elseif strncmp(muscleNames(j), 'perbrev', 4) || strncmp(muscleNames(j), 'perlong', 4)
                
                peronealsIsoForce(peronealsIndex,1) = muscleForce(j);
                peronealsIndex = peronealsIndex + 1;
                
            elseif strncmp(muscleNames(j), 'edl', 3) || strncmp(muscleNames(j), 'ehl', 3)
                
                extensorsIsoForce(extensorsIndex,1) = muscleForce(j);
                extensorsIndex = extensorsIndex + 1;
                           
            else % evlaute the name of the muscle and isoforce
                myMuscleIsoForce = [listOfModelMusclesRootNames{i}, 'IsoForce']; %this is for one "fibre" per msucle cases, 
                eval([myMuscleIsoForce ' = muscleForce(j);']);                   %unlike addmag with 4 [addmagDist,addmagIsch,addmagMid,addmagProx]
                
            end
            
        break;    
            
        end
        
    end
    
end

sumAddmagIsoForce = sum(addmagIsoForce(:,1));
sumGlmaxIsoForce = sum(glmaxIsoForce(:,1));
sumGlmedIsoForce = sum(glmedIsoForce(:,1));
sumGlminIsoForce = sum(glminIsoForce(:,1));
sumPeronealsIsoForce = sum(peronealsIsoForce(:,1));
sumExtensorsIsoForce = sum(extensorsIsoForce(:,1));

for i = 1:size(addmagIsoForce,1)
    
   addmagIsoForce(i,2) = addmagIsoForce(i,1)/sumAddmagIsoForce;
    
end

for i = 1:size(glmaxIsoForce,1)
    
    glmaxIsoForce(i,2) = glmaxIsoForce(i,1)/sumGlmaxIsoForce;
    
end

for i = 1:size(glmedIsoForce,1)
    
    glmedIsoForce(i,2) = glmedIsoForce(i,1)/sumGlmedIsoForce;
    
end

for i = 1:size(glminIsoForce,1)
    
    glminIsoForce(i,2) = glminIsoForce(i,1)/sumGlminIsoForce;
    
end

for i = 1:size(peronealsIsoForce,1)
    
    peronealsIsoForce(i,2) = peronealsIsoForce(i,1)/sumPeronealsIsoForce;
    
end

for i = 1:size(extensorsIsoForce,1)
    
    extensorsIsoForce(i,2) = extensorsIsoForce(i,1)/sumExtensorsIsoForce;
    
end

%% Compute new muscle forces and apply to model
for i = 1:length(listOfHandsfieldMuscleNames)
    
    for j = 1:length(listOfModelMusclesRootNames)
    
        if strcmp(listOfHandsfieldMuscleNames{i}, listOfModelMusclesRootNames{j})
                
            eval(['fracVolInstance = ' listOfModelMusclesRootNames{j};]);
            volumeFractionInstance = vTheory*fracVolInstance;
            strengthFraction = (rho*volumeFractionInstance)/muscleOptFiberLength(j);
            
            if strcmp(listOfHandsfieldMuscleNames{i}, 'addmag')
                
                for k = 1:size(addmagIsoForce,1)
                    
                    aMuscle = muscles.get(k+1); %right side
                    aMuscle.setMaxIsometricForce(strengthFraction*addmagIsoForce(k,2));
                    
                    aMuscle = muscles.get(k+41); %left side
                    aMuscle.setMaxIsometricForce(strengthFraction*addmagIsoForce(k,2));
                    
                end
                
            elseif strcmp(listOfHandsfieldMuscleNames{i}, 'glmax')
                
                for k = 1:size(glmaxIsoForce,1)
                    
                    aMuscle = muscles.get(k+13); %right side
                    aMuscle.setMaxIsometricForce(strengthFraction*glmaxIsoForce(k,2));
                    
                    aMuscle = muscles.get(k+53); %left side
                    aMuscle.setMaxIsometricForce(strengthFraction*glmaxIsoForce(k,2));
                    
                end
                
            elseif strcmp(listOfHandsfieldMuscleNames{i}, 'glmed')
                
                for k = 1:size(glmedIsoForce,1)
                    
                    aMuscle = muscles.get(k+16); %right side
                    aMuscle.setMaxIsometricForce(strengthFraction*glmedIsoForce(k,2));
                    
                    aMuscle = muscles.get(k+56); %left side
                    aMuscle.setMaxIsometricForce(strengthFraction*glmedIsoForce(k,2));
                    
                end
                
            elseif strcmp(listOfHandsfieldMuscleNames{i}, 'glmin')
                
                for k = 1:size(glminIsoForce,1)
                    
                    aMuscle = muscles.get(k+19); %right side
                    aMuscle.setMaxIsometricForce(strengthFraction*glminIsoForce(k,2));
                    
                    aMuscle = muscles.get(k+59); %left side
                    aMuscle.setMaxIsometricForce(strengthFraction*glminIsoForce(k,2));
                    
                end
                
            elseif strcmp(listOfHandsfieldMuscleNames{i}, 'perbrev') || strcmp(listOfHandsfieldMuscleNames{i}, 'perlong')
                
                newPerbrevIsoForce = strengthFraction*peronealsIsoForce(1,2);
                newPerlongIsoForce = strengthFraction*peronealsIsoForce(2,2);
                
                aMuscle = muscles.get(25); %right side
                aMuscle.setMaxIsometricForce(newPerbrevIsoForce);
                aMuscle = muscles.get(26); %right side
                aMuscle.setMaxIsometricForce(newPerlongIsoForce);
                
                aMuscle = muscles.get(65); %left side
                aMuscle.setMaxIsometricForce(newPerbrevIsoForce);
                aMuscle = muscles.get(66); %left side
                aMuscle.setMaxIsometricForce(newPerlongIsoForce);
                
            elseif strcmp(listOfHandsfieldMuscleNames{i}, 'edl') || strcmp(listOfHandsfieldMuscleNames{i}, 'ehl')
                
                newEdlIsoForce = strengthFraction*extensorsIsoForce(1,2);
                newEhlIsoForce = strengthFraction*extensorsIsoForce(2,2);
                
                aMuscle = muscles.get(8); %right side
                aMuscle.setMaxIsometricForce(newEdlIsoForce);
                aMuscle = muscles.get(9); %right side
                aMuscle.setMaxIsometricForce(newEhlIsoForce);
                
                aMuscle = muscles.get(48); %left side
                aMuscle.setMaxIsometricForce(newEdlIsoForce);
                aMuscle = muscles.get(49); %left side
                aMuscle.setMaxIsometricForce(newEhlIsoForce);
                
            else % all other muscles
                
                lengthOfMuscleInstance = length(listOfModelMusclesRootNames{j});
                
                for x = 1:length(muscleNames)
                    
                    if strncmp(muscleNames{x}, listOfModelMusclesRootNames{j}, lengthOfMuscleInstance)
                        
                        muscleName = model.getMuscles().get(muscleNames{x}); %right side
                        muscleClass = char(muscleName.getConcreteClassName);
                        eval(['myMuscle = ' muscleClass '.safeDownCast(muscleName);']);
                        myMuscle.setMaxIsometricForce(strengthFraction);
                        
                        muscleName = model.getMuscles().get(muscleNames{x+40}); %left side
                        muscleClass = char(muscleName.getConcreteClassName);
                        eval(['myMuscle = ' muscleClass '.safeDownCast(muscleName);']);
                        myMuscle.setMaxIsometricForce(strengthFraction);
                        
                        break;
                        
                    end
                end
            end
            
            break;
            
        end
    end
end

%% Print model with new muscle strengths
model.setName([acquisitionInfo.Subject.Code, '_strengthAdjusted']);
modelOut = modelOutName;
model.print(modelOut);
disp(['The new model has been saved at ' modelOut]);

% % % save(['C:\Users\admin_R\Desktop\CEINMS Pipeline Redone v1.2\Examples\ElaboratedData\HAMI01\1Calibration\scaleModels\', 'muscleForce2.mat'], 'muscleForce');
% % % 
% % %     figure_forces_r = muscleNames_figure(41:80);
% % %     figure_forces_r(:,2) = num2cell(muscleForce_original(41:80));    
% % %     figure_forces_r(:,3) = num2cell(muscleForce_error(41:80));
% % %     figure_forces_r(:,4) = num2cell(maxIsoForce(41:80));
% % % 
% % %     figure_forces_r = figure_forces_r(all(cell2mat(figure_forces_r(:,3)) ~= 0,3),:);
% % % 
% % %     h1 = figure('WindowState', 'maximized');
% % %     figure_forces_r_plot = [cell2mat(figure_forces_r(:,2)), cell2mat(figure_forces_r(:,3)), cell2mat(figure_forces_r(:,4))];
% % %     c = categorical(figure_forces_r(:,1));
% % %     h = bar(c,figure_forces_r_plot);
% % %     legend({'Original', 'Old Erroneous' ,'New Correct'}, 'Location', 'northeastoutside');
% % %     ylabel('Force (N)');
% % %     title('Maximum Isometric Muscle Forces R');
% % %     h(1).FaceColor = [211 211 211]/256;
% % %     h(2).FaceColor = [255 64 64]/256;
% % % %     h(3).FaceColor = [112 219 147]/256;
% % %     h(3).FaceColor = 'g';
    
    