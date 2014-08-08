# TEXPRESS v1.0

TEXPRESS is an open-source MATLAB toolbox for processing GDGT LC/MS data.

## Download

[MATLAB File Exchange](http://www.mathworks.com/matlabcentral/fileexchange/47135-texpress-toolbox-v1-0)

## Features

Check out the [TEXPRESS v1.0 wiki](https://github.com/GDGTs/TEXPRESS/wiki) for screenshots of the graphical user interface.

#### Supported File Types
* '.CDF' [[Agilent]](https://github.com/GDGTs/TEXPRESS#agilent-chemstation-macro)

#### Visualization
* TEXPRESS automatically finds common branched and isoprenoid GDGT ion chromatograms during data import

* Quickly plot any combination of samples/total ion chromatograms (TICs) or samples/extracted ion chromatograms (EICs)

#### Analysis
* Global baseline correction using an asymmetric least squares approach

* Peak detection algorithm uses the relative retention order of GDGTs to know where to look first

* Peak integration based on exponentially modified gaussian curve fitting

#### Export
* '.MAT' [MATLAB]

* '.CSV'

#### System Requirements
* MATLAB [>2013a]

* OSX, Windows

## Installation

Unzip TEXPRESS into any directory on the MATLAB search path. To add a new directory to the search path, enter the following in the MATLAB command window:

````matlab
>> savepath folderName/TEXPRESS.m
````

To run TEXPRESS:
````matlab
>> TEXPRESS
````

## Known Issues

To report a bug in the software or to request a new feature, go to: https://github.com/GDGTs/TEXPRESS/issues.

#### Performance Tips

* The 'separate axes' option creates a new axes for each sample selected in the 'Select View' listbox. Keep the number of onscreen samples to <50 for best performance.

* Some users may have difficulty rendering graphics with OpenGL. The default graphics renderer can be changed by using the MATLAB command line:

````matlab
% Option 1: Render w/ 'ZBuffer'
>> set(ans.UserInterface.UIFigure, 'Renderer', 'ZBuffer')

% Option 2: Render w/ 'Painters'
>> set(ans.UserInterface.UIFigure, 'Renderer', 'Painters')
````

## Command Line Options

Run TEXPRESS functions on the MATLAB command line:

#### Import

* Import LC/MS data ('.CDF') into the MATLAB workspace.

````matlab
% Initialize the FileIO class
>> obj = FileIO();

% Import selected data to the workspace
>> LCMSData = obj.SelectFiles('CDF');
````

#### Baseline Correction

* Perform baseline correction of a chromatogram using asymmetric least squares smoothing.

````matlab
% Initialize the SpectraPreprocessing class
>> obj = SpectraPreprocessing();

% Y = chromatographic signal, Smoothness = 10^3-10^8, Asymmetry = 0.1-0.0001
>> LCMSData = obj.BaselineCorrection(y, smoothness, asymmetry);
````

#### Peak Area

* Determine the area of a peak using the exponential modified gaussian model (EMG).

````matlab
% Initialize the PeakArea class
>> obj = PeakArea();

% Find peak area of a known peak
>> PeakData = obj.CalculateArea(xData, yData, retentionTime, 'EMG');
````

## Notes

Converting Agilent '.D' LC/MS data to the '.CDF' file extension:

#### Agilent ChemStation Macro

* Save the following as *ExportCDF.mac* to your *C\...\Chem32\Core directory*:

````cpp
Name ExportCDF

 MSSaveAIA "C:\"
 
endMacro
````

* Call the macro on the Agilent ChemStation command line to convert the currently selected file:

````cpp
>> Macro ExportCDF, go
````

* Or, call the macro automatically after each run using the 'Run Time Checklist':

````cpp
Method --> Run Time Checklist... --> Post-Run Command / Macro --> ExportCDF
````
