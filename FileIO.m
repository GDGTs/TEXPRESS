% TEXPRESS v1.0
%
% James Dillon @ Brown University (2014) 
% www.github.com/GDGTs/TEXPRESS
%
% Class File: FileIO
%
%   I. ImportCDF

classdef FileIO 
    
    properties
        
        FileType
        UserData  
    end

    methods

        % Constructor method
        function obj = FileIO()   

            % Initialize data variables
            obj.UserData.Date = [];
            obj.UserData.TimeValues = [];
            obj.UserData.MassValues = [];
            obj.UserData.TotalIntensityValues = [];
            obj.UserData.IntensityValues = [];
            
            % Initialize file type variable
            obj.FileType = [];
                        
        end

        function obj = SelectFiles(obj, varargin)
            
            if nargin >= 1
                
                % Set file type to input argument
                obj.FileType = varargin{1};
                
                % Prompt user to select files of a valid filetype
                [FileName, FilePath] = uigetfile(obj.FileType, 'Multiselect', 'on');  

                % For single file selections
                if iscell(FileName) == 0

                    % Convert string to cell
                    FileName = {FileName};
                end

                % Set search path to the selected folder
                path(FilePath, path);

                % Determine which import function to execute
                switch obj.FileType

                    case {'.CDF'}    

                        obj = ImportCDF(obj, FileName);
                end
            end
        end

        % Import LC/MS data with '.CDF' extension
        function  obj = ImportCDF(obj, FileName)
            
            % Import selected data into a structure
            for i = 1:length(FileName(1, :))
                
                obj.UserData(i).Name = FileName{1, i}; 
                
                % File Information
                obj.UserData(i).Date = ncreadatt(FileName{1, i}, '/', 'experiment_date_time_stamp');

                % Time Values
                TimeValues = ncread(FileName{1, i}, 'scan_acquisition_time');

                TimeValuesIndex = ncread(FileName{1, i}, 'scan_index');

                obj.UserData(i).TimeValues = TimeValues / 60;
                clear TimeValues

                % Mass Values
                MassValues = ncread(FileName{1, i}, 'mass_values');

                MassValuesMax = ncread(FileName{1, i}, 'mass_range_max');

                obj.UserData(i).MassValues = 1:round(max(MassValuesMax));
                clear MassValuesMax

                % Intensity Values
                IntensityValues = ncread(FileName{1, i}, 'intensity_values');

                TotalIntensityValues = ncread(FileName{1, i}, 'total_intensity');

                obj.UserData(i).TotalIntensityValues = TotalIntensityValues;
                clear TotalIntensityValues

                % Pre-Allocate Memory
                obj.UserData(i).IntensityValues = zeros(...
                    length(TimeValuesIndex), max(round(MassValues)));

                TimeValuesIndex(length(TimeValuesIndex) + 1) = length(MassValues);

                % Transform Data
                for j = 1:length(TimeValuesIndex) - 1

                    obj.UserData(i).IntensityValues(...
                        j, round(MassValues(TimeValuesIndex(j)+1 : TimeValuesIndex(j+1)))) =...
                        (IntensityValues(TimeValuesIndex(j)+1 : TimeValuesIndex(j+1)));
                end

                clear TimeValuesIndex
                clear MassValues
                clear IntensityValues
            end
        end        
    end
end