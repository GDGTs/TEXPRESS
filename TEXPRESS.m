 % TEXPRESS Version 1.0
%
% James Dillon @ Brown University (2014)
%
% www.github.com/GDGTs/TEXPRESS

classdef TEXPRESS < handle
    
    properties
        
        UserInterface
        UserData
  
    end
    
    methods
                
        % Constructor Method
        function obj = TEXPRESS()
            
            % Initialize class for user interface
            obj.UserInterface = ConstructUI();
            
            % Initialize class for data import/export
            obj.UserData = FileIO(); 
            
            % Initialize callbacks for UIMenu
            obj = InitializeUIMenu(obj);
        end
        
        function obj = InitializeUIMenu(obj, varargin)
            
            set(obj.UserInterface.UIMenu.LoadCDF,...
                'Callback', @obj.LoadCDF);
        end
        
        function IntializeUITable(obj, varargin)
            
            % Set table data and list boxes
            obj.UserInterface.InitializeTableData();
            obj.UserInterface.InitializeSampleListBox();
            obj.UserInterface.InitializeIonListBox();

        end
        
        function InitializeUIAxes(obj, varargin)
            
            % Initialize axes with TIC data
            obj.UserInterface.InitializeAxes() 
        end
    
        function LoadCDF(obj, varargin)
            
            % Import LC/MS data with the '.CDF' extension
            obj.UserData = obj.UserData.SelectFiles('.CDF');
            
            % Transfer imported data to user interface
            obj.UserInterface.PlotData = obj.UserData.UserData;
            
            IntializeUITable(obj);
            InitializeUIAxes(obj);
        end
    end
end

