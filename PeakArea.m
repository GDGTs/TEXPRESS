% TEXPRESS v1.0
%
% James Dillon @ Brown University (2014) 
% www.github.com/GDGTs/TEXPRESS
%
% Class File: PeakArea
%
%   I. Exponentially Modified Gaussian

classdef PeakArea
    
    properties
        
        % Input data
        xData
        yData
        
        % Peak properties
        peakStart
        peakEnd
        peakTime
        peakWidth
        peakHeight
        peakArea
        peakExtra

        % Fit data
        fitStart
        fitEnd
        fitX
        fitY
        fitError
        fitErrorValue
    end
    
    methods
   
        % Constructor method
        function obj = PeakArea()
           
            % Initialize input data
            obj.xData = [];
            obj.yData = [];
            
            % Initialize peak properties
            obj.peakStart = [];
            obj.peakEnd = [];
            obj.peakTime = [];
            obj.peakWidth = [];
            obj.peakHeight = [];
            obj.peakArea = [];
            obj.peakExtra = [];
            
            % Initialize fit data
            obj.fitStart = [];
            obj.fitEnd = [];
            obj.fitX = [];
            obj.fitY = [];
            obj.fitError = [];
            obj.fitErrorValue = [];
            
        end
        
        % Prepare input values
        function obj = CalculateArea(obj, varargin)
           
            % Function requires an input value
            if nargin == 5
                
                % Assign inputs to variables
                obj.xData = varargin{1};
                obj.yData = varargin{2};
                obj.peakTime = varargin{3};
                PeakModel = varargin{4};
                
                % Calculate peak width
                obj = PeakWidth(obj);

                % Find x index start/end
                obj.fitStart = find(...
                    obj.xData >= obj.peakTime - 3 * obj.peakWidth, 1);                
                obj.fitEnd = find(...
                    obj.xData >= obj.peakTime + 3 * obj.peakWidth, 1);

                % Pre-allocate space for peak fit
                obj.fitX = linspace(...
                    obj.xData(obj.fitStart),...
                    obj.xData(obj.fitEnd),...
                    obj.fitEnd - obj.fitStart);
                
                obj.fitX = transpose(obj.fitX);
                
                % Determine which peak model to send data to
                switch PeakModel
                    
                    % Exponentially Modified Gaussian (EMG)
                    case 'EMG'
                        
                        % Step 1: Find peak limits
                        obj = PeakLimits(obj);
    
                        % Step 2: Find optimal peak width and exponential factor
                        obj = ExponentialGaussianParameters(obj);
                        
                        % Step 3: Calculate y values of fitted peak
                        obj = ExponentialGaussian(obj);
                        
                        % Step 4: Calculate error of fit
                        obj = PeakFitError(obj);
                end
            end
        end
        
        % Exponentially Modified Gaussian (EMG)
        function obj = ExponentialGaussian(obj, varargin)
  
            % Time values
            x = transpose(obj.xData(obj.fitStart:obj.fitEnd));
            % Peak retention time
            u = obj.peakTime;
            % Peak height
            h = obj.peakHeight;
            % Peak width
            w = obj.peakWidth;
            % Peak exponential factor
            e = obj.peakExtra;
            
            % Predicted y values using the exponential modified gaussian equation
            obj.fitY = ...
                h * exp(-((u-x).^2) ./ (2*(w^2))) .* ...
                (w/e) .* ((pi/2)^0.5) .* ...
                erfcx((1/(2^0.5)) .* (((u-x)./w) + (w/e)));
            
            % Transpose vector
            obj.fitY = transpose(obj.fitY);

            % Calculate area
            EMG = @(x,h,u,w,e) ...
                h * exp(-((u-x).^2) ./ (2*(w^2))) .* ...
                (w/e) .* ((pi/2)^0.5) .* ...
                erfcx((1/(2^0.5)) .* (((u-x)./w) + (w/e)));

            clear x
            
            % Integral of EMG function
            obj.peakArea = integral(@(x) EMG(x,h,u,w,e),...
                obj.xData(obj.fitStart),obj.xData(obj.fitEnd));

            clear x y u h w
        end
        
        % Exponentially Modified Gaussian (EMG) simplex optimization
        function obj = ExponentialGaussianParameters(obj, varargin)
            
            % Time values
            x = transpose(obj.xData(obj.peakStart:obj.peakEnd));
            % Experimental y values
            y = transpose(obj.yData(obj.peakStart:obj.peakEnd));
            % Peak retention time
            u = obj.peakTime;
            % Peak height
            h = obj.peakHeight;
            % Peak width
            w = obj.peakWidth;
                        
            % Sum of squares using the exponential modified gaussian equation
            SumSquares = @(p) sum((y - ...
                (h * exp(-((u-x).^2) ./ (2*(p(2)^2))) .* ...
                (p(2)/p(1)) .* ((pi/2)^0.5) .* ...
                erfcx((1/(2^0.5)) .* (((u-x)./p(2)) + (p(2)/p(1)))))).^2);
            
            % Optimize the following parameters: exponential factor, peak width
            p = [0.5, w];
            
            % Simplex based optimization of parameters by minimizing sum of squares
            p = fminsearch(SumSquares, p);
                        
            % Reassign optimized parameters
            obj.peakExtra = p(1);
            obj.peakWidth = p(2);
                                    
            clear x y u h w
            clear p
            clear SumSquares
        end
        
        % Determine peak width at half height
        function obj = PeakWidth(obj, varargin)
        
            % Determine row index of peak height
            peakIndex = find(obj.xData >= obj.peakTime, 1);
        
            % Set peak height variable
            obj.peakHeight = obj.yData(peakIndex);
        
            % Find peak width at half height
            HalfHeight = obj.peakHeight / 2;
        
            % Get index of right and left half height values
            rightIndex = peakIndex + find(obj.yData(peakIndex:end) <= HalfHeight, 1);
            leftIndex = peakIndex - find(flipud(obj.yData(1:peakIndex)) <= HalfHeight, 1);
        
            % Get distance from center to left/right side of peak
            rightDistance = obj.xData(rightIndex) - obj.xData(peakIndex);
            leftDistance = obj.xData(peakIndex) - obj.xData(leftIndex);
        
            % Calculate peak asymmetry at half height
            peakAsymmetry = rightDistance / leftDistance;
        
            % If asymmetry is within a valid range
            if peakAsymmetry < 2 && peakAsymmetry > 0.5
            
                % Calculate peak width
                obj.peakWidth = rightDistance + leftDistance;
            
            % If asymmetry is invalid due to left side
            elseif leftDistance * 2 > rightDistance
            
                % Assume left distance is invalid
                obj.peakWidth = rightDistance * 2;
            
                % If asymmetry is invalid due to right side
            elseif rightDistance * 2 > leftDistance
            
                % Assume right distance is invalid
                obj.peakWidth = leftDistance * 2;
            end
        
            clear peakIndex 
            clear peakAsymmetry 
            clear halfHeight
            clear rightIndex leftIndex
            clear rightDistance leftDistance
        end
        
        % Determine peak limits for integration
        function obj = PeakLimits(obj, varargin)
        
            % Peak values
            LeftIndex = obj.fitStart;
            RightIndex = obj.fitEnd;
            
            % Derivative of signal vector
            dy = diff(obj.yData(LeftIndex:RightIndex));

            % Find index values of derivative max/min
            [~, dyMaxIndex] = max(dy);
            [~, dyMinIndex] = min(dy);
            
            % Look for peak edge in the following vector
            LeftSearch = flipud(dy(1:dyMaxIndex));
            RightSearch = dy(dyMinIndex:end);
            
            % Find where peak starts by finding dy @ 0
            dyLeftIndex = dyMaxIndex - find(LeftSearch <= 1, 1);
            dyRightIndex = dyMinIndex + find(RightSearch >= -1, 1);
            
            % Find peak heights at peak limits
            LeftHeight = obj.yData(dyLeftIndex + LeftIndex);
            RightHeight = obj.yData(dyRightIndex + LeftIndex);
            
            % Reassign new peak limits
            obj.peakStart = LeftIndex + dyLeftIndex;
            obj.peakEnd = LeftIndex + dyRightIndex;
            
            % Account for any coelution
            if RightHeight > LeftHeight
                
                % Set left index where peak height is equivalent to right
                obj.peakStart = LeftIndex + find(obj.yData(LeftIndex:RightIndex) >= RightHeight, 1);
            
            elseif LeftHeight > RightHeight
                
                % Set left index where peak height is equivalent to right
                obj.peakEnd = RightIndex - find(flipud(obj.yData(LeftIndex:RightIndex)) >= LeftHeight, 1); 
            end
            
            clear dy
            clear LeftSearch RightSearch
            clear dyLeftIndex dyRightIndex 
            clear dyMaxIndex dyMinIndex
            clear LeftIndex RightIndex 
            clear LeftHeight RightHeight
        end
        
        % Find fit error of calculated peak
        function obj = PeakFitError(obj, varargin)
           
            % Subtract calculated from experimental
            obj.fitError = obj.yData(obj.fitStart:obj.fitEnd) - obj.fitY;
            
            % Reshape vectors
            obj.fitY = obj.fitY(1:end-1);
            obj.fitError = obj.fitError(1:end-1);
             
            % Calculate root mean square deviation
            ResidualsSquared = obj.fitError.^2;
            ResidualsMean = mean(ResidualsSquared);
            RMSD = sqrt(ResidualsMean / length(ResidualsSquared));
            
            % Normalized RMSD
            obj.fitErrorValue = (RMSD / (obj.peakHeight - min(obj.fitY))) * 100;
        end
    end
end