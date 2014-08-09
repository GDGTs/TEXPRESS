% TEXPRESS v1.0
%
% James Dillon @ Brown University (2014) 
% www.github.com/GDGTs/TEXPRESS
%
% Class File: PeakProxy
%

classdef PeakProxy
    
    properties
        
        % Class data
        id
        options
        
        % Available data
        ions
        area
        
        % Available proxy status
        status
        
        % Available branched GDGT proxies
        MBT
        CBT        
        MAAT
        pH
        
        % Available isoprenoid GDGT proxies
        TEX86
        SST
        
        % Available mixed GDGT proxies
        BIT
        
        % Calibration models
        proxy
        models
        index
        
    end

    methods
   
        % Constructor method
        function obj = PeakProxy()
           
            % Store data
            obj.id = [];
            obj.options = [];
            
            % Available ions to perform functions on 
            obj.ions= [];
            obj.area = [];
            
            % Availability object
            obj.status = [];
            
            % Branched GDGT proxies
            obj.MBT = [];
            obj.CBT = [];
            obj.MAAT = [];
            obj.pH = [];
            
            % Isoprenoid GDGT proxies
            obj.TEX86 = [];
            obj.SST = [];
            
            % Mixed GDGT proxies
            obj.BIT = [];
            
            % Store available calibration models
            obj.proxy = [];
            obj.models = [];
            obj.index = [];
            
            % Initialize status of each proxy and calibration models
            obj = InitializeStatus(obj);
            obj = InitializeModels(obj);
            obj = InitializeIndex(obj);
        end
        
        % Input function to determine which values to calculate
        function obj = CalculateProxy(obj, varargin)
            
            % Input requires UITable
            if nargin == 4
                
                % Input order is ions, areas, options
                obj.ions = varargin{1};
                obj.area = varargin{2};
                obj.options = varargin{3};
                
                % Check formatting of input values
                obj = CheckInput(obj);
                        
                % Update availability status of proxies
                obj = CheckStatus(obj);
                
                % Determine which proxies to calculate
                switch obj.options
                
                    % If option input was branched GDGTs
                    case 'branched'
                        obj = BranchedGDGT(obj);
                        
                    case 'isoprenoid'
                        obj = IsoprenoidGDGT(obj);
                        
                    case 'all'
                        obj = BranchedGDGT(obj);
                        obj = IsoprenoidGDGT(obj);
                end
            end
        end
        
        % Branched GDGT proxies
        function obj = BranchedGDGT(obj, varargin)
            
            % Calculate CBT
            if strcmp(obj.status.CBT, 'on') == 1
                
                % Calculate numerator and denominator for CBT
                totalB = obj.area(obj.id.Ib) + obj.area(obj.id.IIb);
                totalA = obj.area(obj.id.Ia) + obj.area(obj.id.IIa);
                
                % Calculate CBT index
                obj.CBT = -log10(totalB/totalA);
            end
            
            % Calculate MBT
            if strcmp(obj.status.MBT, 'on') == 1
                
                switch obj.proxy.MAAT
                    
                    case 'weijers2007'
                
                        % Calculate numerator and denominator for MBT
                        totalI = obj.area(obj.id.Ia) + obj.area(obj.id.Ib) + obj.area(obj.id.Ic); 
                        totalII = obj.area(obj.id.IIa) + obj.area(obj.id.IIb) + obj.area(obj.id.IIc); 
                        totalIII = obj.area(obj.id.IIIa) + obj.area(obj.id.IIIb) + obj.area(obj.id.IIIc); 
               
                        % Calculate MBT
                        obj.MBT = totalI / (totalI + totalII + totalIII);
                    
                    case 'peterse2012'
                        % Calculate numerator and denominator for MBT (exclude IIc and IIIc)
                        totalI = obj.area(obj.id.Ia) + obj.area(obj.id.Ib) + obj.area(obj.id.Ic); 
                        totalII = obj.area(obj.id.IIa) + obj.area(obj.id.IIb) + obj.area(obj.id.IIc);
                        totalIII = obj.area(obj.id.IIIa);
                        
                        % Calculate MBT
                        obj.MBT = totalI / (totalI + totalII + totalIII);
                end
            end
            
            % Calculate BIT
            if strcmp(obj.status.BIT, 'on') == 1
               
                % Calculate numerator and denominator for BIT
                totalA = obj.area(obj.id.Ia) + obj.area(obj.id.IIa) + obj.area(obj.id.IIIa);
                
                % Calculate BIT
                obj.BIT = totalA / (obj.area(obj.id.iCren) + totalA);
            end
        
            % Calculate pH
            if strcmp(obj.status.pH, 'on') == 1
                
                switch obj.proxy.pH
                
                    case 'weijers2007'
                        % Weijers et. al. GCA 71 (2007) 703.
                        obj.pH = (3.22 - obj.CBT) / 0.38;
                        
                    case 'tierney2010'
                        % Tierney et. al. GCA 74 (2010) 4902.
                        obj.pH = 10.32 - 3.30 * obj.CBT;
                        
                    case 'peterse2012'
                        % Peterse et.al. GCA 96 (2012) 215.
                        obj.pH = 7.90 - 1.97 * obj.CBT;
                end
            end
            
            % Calculate MAAT
            if strcmp(obj.status.MAAT, 'on') == 1
                
                switch obj.proxy.MAAT
                    
                    case 'weijers2007'
                        % Weijers et. al. GCA 71 (2007) 703.
                        obj.MAAT = (obj.MBT - 0.122 - 0.187 * obj.CBT) / 0.02;
                        
                    case 'tierney2010'
                        % Tierney et. al. GCA 74 (2010) 4902.
                        obj.MAAT = 11.84 - 9.32 * obj.CBT + 32.54 * obj.MBT;
                        
                    case 'peterse2012'
                        % Peterse et.al. GCA 96 (2012) 215.
                        obj.MAAT = 0.81 - 5.67 * obj.CBT + 31.0 + obj.MBT;
                end
            end
        end

        % Isoprenoid GDGT proxies
        function obj = IsoprenoidGDGT(obj, varargin)

            % Calculate TEX86
            if strcmp(obj.status.TEX86, 'on') == 1

                switch obj.proxy.SST
                    
                    case 'schouten2002'
                        % Calculate numerator and denominator for TEX86
                        numerator = obj.area(obj.id.i2) + obj.area(obj.id.i3) + obj.area(obj.id.iCren);
                        denominator = obj.area(obj.id.i1) + numerator;

                        % Calculate TEX86
                        obj.TEX86 = numerator / denominator; 
                        
                    case 'kim2008'
                        % Calculate numerator and denominator for TEX86
                        numerator = obj.area(obj.id.i2) + obj.area(obj.id.i3) + obj.area(obj.id.iCren);
                        denominator = obj.area(obj.id.i1) + numerator;

                        % Calculate TEX86
                        obj.TEX86 = numerator / denominator; 
                        
                    case 'liu2009'
                        % Calculate numerator and denominator for TEX86
                        numerator = obj.area(obj.id.i2) + obj.area(obj.id.i3) + obj.area(obj.id.iCren);
                        denominator = obj.area(obj.id.i1) + numerator;

                        % Calculate TEX86
                        obj.TEX86 = numerator / denominator; 
                        
                    case 'kim2010a'
                        % Calculate numerator and denominator for TEX86
                        numerator = obj.area(obj.id.i2);
                        denominator = obj.area(obj.id.i1) + obj.area(obj.id.i2) + obj.area(obj.id.i3);

                        % Calculate TEX86
                        obj.TEX86 = log10(numerator / denominator); 
                        
                    case 'kim2010b'
                        % Calculate numerator and denominator for TEX86
                        numerator = obj.area(obj.id.i2) + obj.area(obj.id.i3) + obj.area(obj.id.iCren);
                        denominator = obj.area(obj.id.i1) + numerator;

                        % Calculate TEX86
                        obj.TEX86 = log10(numerator / denominator); 
                        
                    case 'powers2010'
                        % Calculate numerator and denominator for TEX86
                        numerator = obj.area(obj.id.i2) + obj.area(obj.id.i3) + obj.area(obj.id.iCren);
                        denominator = obj.area(obj.id.i1) + numerator;

                        % Calculate TEX86
                        obj.TEX86 = numerator / denominator;
                end
            end

            % Calculate SST
            if strcmp(obj.status.SST, 'on') == 1

                switch obj.proxy.SST

                    case 'schouten2002'
                        % Schouten et. al. EPSL 204 (2002) 265. (r2 = 0.92, n = 43) 
                        obj.SST = (obj.TEX86 - 0.28) / 0.015;

                    case 'kim2008'
                        % Kim et. al. GCA 72 (2008) 1154. (r2 = 0.94, n = 223)                        
                        obj.SST = obj.TEX86 * 56.2 - 10.78;
                        
                    case 'liu2009'
                        % Liu et. al. Science 323 (2009) 1187. (r2 = 0.82, n = 287)
                        obj.SST = 50.475 - 16.332 * (1 / obj.TEX86);
                        
                    case 'powers2010'
                        % Powers et. al. Org. Geochem. 41 (2010) 404 (r2 = 0.86, n = 12)
                        obj.SST = obj.TEX86 * 55.2 - 14.0;
                        
                    case 'kim2010a'
                        % Kim et. al. GCA 74 (2010) 4639. (r2 = 0.86, n = 396)                        
                        obj.SST = obj.TEX86 * 67.5 + 49.9;
                        
                    case 'kim2010b'
                        % Kim et. al. GCA 74 (2010) 4639. (r2 = 0.86, n = 396)                        
                        obj.SST = obj.TEX86 * 68.4 - 38.6;
                end
            end
        end

        % Set initial status of proxy availability
        function obj = InitializeStatus(obj, varargin)
           
            % Set initial status to 'off' for all branched GDGT proxies
            obj.status.CBT = 'off';
            obj.status.MBT = 'off';
            obj.status.MAAT = 'off';
            obj.status.pH = 'off';
            
            % Set initial status to 'off' for all isoprenoid GDGT proxies
            obj.status.TEX86 = 'off';
            obj.status.SST = 'off';
            
            % Set initial status to 'off' for all mixed GDGT proxies
            obj.status.BIT = 'off';
        end
        
        % Check if ions are available to calculate proxies
        function obj = CheckStatus(obj, varargin)
            
            % Find common branched GDGT ions (if available)
            obj.id.Ia = find(obj.ions == 1022,1);
            obj.id.Ib = find(obj.ions == 1020,1);
            obj.id.Ic = find(obj.ions == 1018,1);
            
            obj.id.IIa = find(obj.ions == 1036,1);
            obj.id.IIb = find(obj.ions == 1034,1);
            obj.id.IIc = find(obj.ions == 1032,1);
            
            obj.id.IIIa = find(obj.ions == 1050,1);
            obj.id.IIIb = find(obj.ions == 1048,1);
            obj.id.IIIc = find(obj.ions == 1046,1);
            
            % Find common isoprenoid GDGT ions (if available)
            obj.id.i1 = find(obj.ions == 1300,1);
            obj.id.i2 = find(obj.ions == 1298,1);
            obj.id.i3 = find(obj.ions == 1296,1);
            obj.id.iCren = find(obj.ions == 1292,1);
            
            % Determine if area values for CBT are available
            if isempty(obj.id.Ia) == 0 && isempty(obj.id.Ib) == 0 && ...
                    isempty(obj.id.IIa) == 0 && isempty(obj.id.IIb) == 0
                
                % Set status to 'on'
                obj.status.CBT = 'on';
                obj.status.pH = 'on';
            end
            
            % Determine if area values for MBT are available
            if strcmp(obj.status.CBT, 'on') == 1 && isempty(obj.id.Ic) == 0 && ...
                    isempty(obj.id.IIc) == 0 && isempty(obj.id.IIIa) == 0 && ...
                    isempty(obj.id.IIIb) == 0 && isempty(obj.id.IIIc) == 0
                
                % Set status to 'on'
                obj.status.MBT = 'on';
                obj.status.MAAT = 'on';
            end
            
            % Determine if area values for MBT are available
            if isempty(obj.id.Ia) == 0 && isempty(obj.id.IIa) == 0 && ...
                    isempty(obj.id.IIIa) == 0 && isempty(obj.id.iCren) == 0
                
                % Set status to 'on'
                obj.status.BIT = 'on';
            end
            
            % Determine if area values for TEX86 are available
            if isempty(obj.id.i1) == 0 && isempty(obj.id.i2) && ...
                    isempty(obj.id.i3) && isempty(obj.id.iCren) == 0
                
                % Set status to 'on'
                obj.status.TEX86 = 'on';
                obj.status.SST = 'on';
            end
        end
        
        % Check format of input values
        function obj = CheckInput(obj, varargin)
           
            % Check format of obj.area
            if iscell(obj.area) == 1

                % Check for empty cells
                for i = 1:length(obj.area)
                    if isempty(obj.area{i}) == 1
                        if iscell(obj.ions) == 1
                            obj.ions{i} = [];
                        end
                    end
                end

                % Check if end is missing values
                if length(obj.ions) > length(obj.area)
                    
                    % Number of values to remove
                    missing = length(obj.ions) - length(obj.area) - 1;
                    
                    % Remove ion values
                    obj.ions(end-missing:end) = [];
                end
                    
                % Convert cell to matrix
                obj.area = cell2mat(obj.area);
            end
            
            % Check format of obj.ions
            if iscell(obj.ions) == 1

                % Convert cell to matrix
                obj.ions = cell2mat(obj.ions);
            end
            
            % Check format of obj.options
            if strcmp(obj.options, 'branched') == 0 && ...
               strcmp(obj.options, 'isoprenoid') == 0 && ...
               strcmp(obj.options, 'all') == 0
           
                obj.options = 'all';
            end
        end
        
        % Initialize calibration models to apply
        function obj = InitializeModels(obj, varargin)
        
            % Branched GDGT calibration defaults
            obj.proxy.pH = 'weijers2007';
            obj.proxy.MAAT = 'weijers2007';
            
            % Isoprenoid GDGT calibration defauls
            obj.proxy.SST = 'kim2008';
            
            % Branched GDGT models
            obj.models(1).name = 'weijers2007';
                obj.models(1).type  = 'pH';
                obj.models(1).parent = 'branched';
                obj.models(1).children = 'CBT';
                obj.models(1).environment = 'soil';
                obj.models(1).equation = 'pH = (3.22 - CBT) / 0.38';
                obj.models(1).reference = 'Weijers et. al. GCA 71 (2007) 703.';
                
            obj.models(2).name = 'weijers2007';
                obj.models(2).type  = 'MAAT';
                obj.models(2).parent = 'branched';
                obj.models(2).children = 'MBT, CBT';
                obj.models(2).environment = 'soil';
                obj.models(2).equation = 'MAAT = 11.84 - 9.32 * CBT + 32.54 * MBT';
                obj.models(2).reference = 'Weijers et. al. GCA 71 (2007) 703.';
                
            obj.models(3).name = 'tierney2010';
                obj.models(3).type  = 'pH';
                obj.models(3).parent = 'branched';
                obj.models(3).children = 'CBT';
                obj.models(3).environment = 'lake';
                obj.models(3).equation = 'pH = 10.32 - 3.30 * CBT';
                obj.models(3).reference = 'Tierney et. al. GCA 74 (2010) 4902.';
                
            obj.models(4).name = 'tierney2010';
                obj.models(4).type = 'MAAT';
                obj.models(4).parent = 'branched';
                obj.models(4).children = 'MBT, CBT';
                obj.models(4).environment = 'lake';
                obj.models(4).equation = 'MAAT = 11.84 - 9.32 * CBT + 32.54 * MBT';
                obj.models(4).reference = 'Tierney et. al. GCA 74 (2010) 4902.';
                
            obj.models(5).name = 'peterse2012';
                obj.models(5).type  = 'pH';
                obj.models(5).parent = 'branched';
                obj.models(5).children = 'CBT';
                obj.models(5).environment = 'soil';
                obj.models(5).equation = 'pH = 7.90 - 1.97 * CBT';
                obj.models(5).reference = 'Peterse et.al. GCA 96 (2012) 215.';
                
            obj.models(6).name = 'peterse2012';
                obj.models(6).type  = 'MAAT';
                obj.models(6).parent = 'branched';
                obj.models(6).children = 'MBT, CBT';
                obj.models(6).environment = 'soil';
                obj.models(6).equation = 'MAAT = 0.81 - 5.67 * CBT + 31.0 + MBT''';
                obj.models(6).reference = 'Peterse et.al. GCA 96 (2012) 215.';
             
            % Isoprenoid GDGT models
            obj.models(7).name = 'schouten2002';
                obj.models(7).type  = 'SST';
                obj.models(7).parent = 'isoprenoid';
                obj.models(7).children = 'TEX86';
                obj.models(7).environment = 'ocean';
                obj.models(7).equation = 'SST = (TEX86 - 0.28) / 0.015';
                obj.models(7).range_low = 0;
                obj.models(7).range_high = 30;
                obj.models(7).reference = 'Schouten et. al. EPSL 204 (2002) 265.';    
                
            obj.models(8).name = 'kim2008';
                obj.models(8).type  = 'SST';
                obj.models(8).parent = 'isoprenoid';
                obj.models(8).children = 'TEX86';
                obj.models(8).environment = 'ocean';
                obj.models(8).equation = 'SST = TEX86 * 56.2 - 10.78';
                obj.models(8).range_low = 5;
                obj.models(8).range_high = 30;
                obj.models(8).reference = 'Kim et. al. GCA 72 (2008) 1154.';
                
            obj.models(9).name = 'liu2009';
                obj.models(9).type  = 'SST';
                obj.models(9).parent = 'isoprenoid';
                obj.models(9).children = 'TEX86';
                obj.models(9).environment = 'ocean';
                obj.models(9).equation = 'SST = 50.475 - 16.332 * (1 / TEX86)';
                obj.models(9).range_low = -3;
                obj.models(9).range_high = 30;
                obj.models(9).reference = 'Liu et. al. Science 323 (2009) 1187.';
                
            obj.models(10).name = 'kim2010a';
                obj.models(10).type  = 'SST';
                obj.models(10).parent = 'isoprenoid';
                obj.models(10).children = 'TEX86_L';
                obj.models(10).environment = 'ocean';
                obj.models(10).equation = 'SST = TEX86_L * 67.5 + 49.9';
                obj.models(10).range_low = -3;
                obj.models(10).range_high = 30;
                obj.models(10).reference = 'Kim et. al. GCA 74 (2010) 4639.';   
                
            obj.models(11).name = 'kim2010b';
                obj.models(11).type  = 'SST';
                obj.models(11).parent = 'isoprenoid';
                obj.models(11).children = 'TEX86_H';
                obj.models(11).environment = 'ocean';
                obj.models(11).equation = 'SST = TEX86_H * 68.4 - 38.6';
                obj.models(11).range_low = 5;
                obj.models(11).range_high = 30;
                obj.models(11).reference = 'Kim et. al. GCA 74 (2010) 4639.';   
                
            obj.models(12).name = 'powers2010';
                obj.models(12).type  = 'SST';
                obj.models(12).parent = 'isoprenoid';
                obj.models(12).children = 'TEX86';
                obj.models(12).environment = 'lake';
                obj.models(12).equation = 'SST = TEX86 * 55.2 - 14.0';
                obj.models(12).range_low = 5;
                obj.models(12).range_high = 30;
                obj.models(12).reference = 'Powers et. al. Org. Geochem. 41 (2010) 404.';
    end
        
        % Initialize index models to apply
        function obj = InitializeIndex(obj, varargin)
        
            % Branched GDGT index defaults
            obj.index(1).name = 'MBT';
                obj.index(1).type  = 'MAAT';
                obj.index(1).parent = 'branched';
                obj.index(1).equation = '(Ia+Ib+Ic) / (Ia+Ib+Ic + IIa+IIb+IIc + IIIa+IIIb+IIIc)';
                
            obj.index(2).name = 'MBT''';
                obj.index(2).type  = 'MAAT';
                obj.index(2).parent = 'branched';
                obj.index(2).equation = '(Ia+Ib+Ic) / (Ia+Ib+Ic + IIa+IIb+IIc + IIIa)';
                
            obj.index(3).name = 'CBT';
                obj.index(3).type  = 'pH';
                obj.index(3).parent = 'branched';
                obj.index(3).equation = '-log((Ib+IIb) / (Ia+IIa))';

            obj.index(4).name = 'BIT';
                obj.index(4).type  = 'Soil OM';
                obj.index(4).parent = 'branched';
                obj.index(4).equation = '(Ia+IIa+IIIa) / (Ia+IIa+IIIa + Cren.)';

            % Isoprenoid GDGT index defaults
            obj.index(5).name = 'TEX86';
                obj.index(5).type  = 'SST';
                obj.index(5).parent = 'isoprenoid';
                obj.index(5).equation = '(2+3 + Cren.'') / (1+2+3 + Cren.'')';

            obj.index(6).name = 'TEX86_L';
                obj.index(6).type = 'SST';
                obj.index(6).parent = 'isoprenoid';
                obj.index(6).equation = 'log((2) / (1+2+3))';

            obj.index(7).name = 'TEX86_H';
                obj.index(7).type  = 'SST';
                obj.index(7).parent = 'isoprenoid';
                obj.index(7).equation = 'log((2+3 + Cren.'') / (1+2+3 + Cren.''))';
        end
    end
end