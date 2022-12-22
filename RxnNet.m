%%-------------------------------------------------------------------------
% Rodrigo García-Tejera
% Contact info: rodrigo.garcia@ed.ac.uk ; rgarcia@fisica.edu.uy
% Affiliations: 
% Centre for Regenerative Medicine, University of Edinburgh, UK. 
% Created: 19-12-2022
% Updated: 20-12-2022
%%-------------------------------------------------------------------------


classdef RxnNet 
    
    properties
        stoichiometricMatrix
        rateVector
        propensityFunctions
    end

    methods

    %Constructor for BirthDeathProcess objects
    function obj= RxnNet(stoichiometricMatrix,rateVector,propensityFunctions)
            
            obj.stoichiometricMatrix=stoichiometricMatrix;
            obj.rateVector=rateVector;
            obj.propensityFunctions=propensityFunctions;
        
     end     
        



%Runs the stochastic simulation algorithm (SSA) for the chemical reaction 
%network defined in obj [1]. 
%
%Inputs:  
% obj - chemical reaction network object
% timeVector - vector of times where the molecule numbers are going to be
%              evaluated
% initialNumbers - vector of initial molecule numbers
% extinction - boolean variable to state whether there is an absorbing
% boundary at the zero molecule state. To speed up the simulation if 
% extinction=true, simulation ends whenever the sum of all molecule numbers
% is zero.  
%
% Ouput:
% trajectories - matrix with stochastic trajectories, one column per 
% species. size(trajectories) = numel(timeVector) x number of species.
%
% Note: Beware of the sampling rate introduced by timeVector. Oversampling
% can result in a slight bias.
%
% [1] Gillespie, Daniel T. "Exact stochastic simulation of coupled chemical
% reactions." The journal of physical chemistry 81.25 (1977): 2340-2361.
  

    function trajectories=SSA(obj,timeVector,initialNumbers,extinction)
            

        %number of species in the reaction network
        [numberOfSpecies,~]=size(obj.stoichiometricMatrix);

        
        %memory allocation for the stochastic trajectories
        trajectories=nan(numel(timeVector),numberOfSpecies);
        
        %initial conditions
        moleculeNumbers=initialNumbers;
        currentTime=0;
        
        %--------------ITERATION OF THE GILLESPIE ALGORITHM----------------
        
        %index for the time vector
        index=1; 
        
        while index<=numel(timeVector)

            %updates propensities  
            rawPropensities=obj.rateVector.*obj.propensityFunctions(moleculeNumbers);
                             
            %updates normalisation factor
            propensityNormalisation=sum(rawPropensities);                     

            %picks two random numbers uniformly distributed in [0,1] 
            draw=rand(1,2);    
                
            %sets time in which a reaction takes place 
            currentTime=currentTime-log(draw(1))/propensityNormalisation;
            
            %if new time exceeds evaluation time, assign current molecule 
            % numbers to the evaluation time and continue to the next step
            if currentTime > timeVector(index)
               trajectories(index,:)=moleculeNumbers;
               index=index+1;
            end
            
            %chooses which reaction is going to take place
            cumProbabilities=cumsum(rawPropensities/propensityNormalisation);
            I=find(cumProbabilities-draw(2)>0,1,'first'); 
            
            %updates molecule numbers
            %if numberOfSpecies > 1

            %keyboard

            %end
            moleculeNumbers=moleculeNumbers+obj.stoichiometricMatrix(:,I)';
            
            
            %the following is just to speed up the simulation when n=0 is an
            %absorbing boundary. 
            if extinction==true
                if sum(moleculeNumbers)==0 
                    keyboard
                    trajectories(index:end,:)=0;
                    index=numel(timeVector)+1;
                end
            end

        end
    end

              
    end
end
