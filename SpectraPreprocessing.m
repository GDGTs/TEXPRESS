% TEXPRESS v1.0
%
% James Dillon @ Brown University (2014) 
% www.github.com/GDGTs/TEXPRESS
%
% Class File: SpectraPreprocessing
%
%   I. BaselineCorrection
%       -Reference: P.H.C. Eilers, Analytical Chemistry, 75 (2003) 3631

classdef SpectraPreprocessing < handle
    
    properties
        
        BaselineData
        BaselineProperties
    end
    
    methods
        
        % Constructor method
        function obj = SpectraPreprocessing()
        
            % Initialize variable to store baseline data
            obj.BaselineData.Baseline = [];
            obj.BaselineData.BaselineCorrected = [];
            
            % Initialize variables to store input parameters
            obj.BaselineProperties.Smoothness = [];
            obj.BaselineProperties.Asymmetry = [];
        end
        
        % Asymmetric least squares baseline correction
        function obj = BaselineCorrection(obj, varargin)

            % If only input is the data vector
            if nargin >= 3

                % Assign YData
                yData = varargin{1};
                
                % Set smoothness and asymmetry to default values
                obj.BaselineProperties.Smoothness = varargin{2};
                obj.BaselineProperties.Asymmetry = varargin{3};
                
                for i = 1:length(yData(1,:))
                    
                    % Get length of y vector
                    LengthY = length(yData(:,1));
                    
                    % Initialize variables needed for calculation
                    DifferenceMatrix = diff(speye(LengthY), 2);
                    Weights = ones(LengthY, 1);
                    
                    % Pre-allocate memory for baseline
                    obj.BaselineData.Baseline(:,i) = zeros(LengthY, 1);
                    
                    % Make sure data is positive values
                    if min(yData(:,i)) >= 0 && max(yData(:,i)) > 0
                    
                        % Minimize baseline calculations with x number of iterations
                        for j = 1:10
                            
                            % Sparse diagonal matrix
                            WeightsDiagonal = spdiags(Weights, 0, LengthY, LengthY);
                            
                            % Cholesky factorization
                            CholeskyFactor = chol(...
                                WeightsDiagonal + obj.BaselineProperties.Smoothness *...
                                DifferenceMatrix' * DifferenceMatrix);
                            
                            % Left matrix divide, multiply matrices
                            obj.BaselineData.Baseline(:,i) = ...
                                CholeskyFactor \ (CholeskyFactor' \ (Weights .* yData(:,i)));
                            
                            % Reassign weights
                            Weights = ...
                                obj.BaselineProperties.Asymmetry *...
                                (yData(:,i) > obj.BaselineData.Baseline(:,i)) +...
                                (1 - obj.BaselineProperties.Asymmetry) *...
                                (yData(:,i) < obj.BaselineData.Baseline(:,i));
                        end
                        
                        % Subtract calculated baseline from data
                        obj.BaselineData.BaselineCorrected(:,i) = ...
                            minus(yData(:,i), obj.BaselineData.Baseline(:,i));
                        
                    % If data is all zeros or negative    
                    else
                        obj.BaselineData.BaselineCorrected(:,i) = zeros(length(yData),1);
                        obj.BaselineData.Baseline(:,i) = zeros(length(yData), 1);
                    end
                end
            end
        end
    end
end
            
        