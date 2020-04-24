
Analysis scripts of calcium imaging in locomoting mice 
==============================

The repository contains custom written Python scripts used to extract and analyze calcium fluorescent traces during forced locomotion.
Calcium imaging was performed before and after application of a mGluR anatagonist.  The here published Python scripts implement 
the steps and calculation to obtain **locomotion fluorescence** as depicted in Fig. 3  of the publication below. The auxiliary graphics
functions assembled in ```generatePublicationFigure``` are used by multiple scripts below. 

For more details, please refer to :

Jinn Bao, Michael Graupner, Guadalupe Astorga, Thibault Collin, Abdelali Jalil, Dwi Wahyu Indriati, Jonathan Bradley, 
Ryuichi Shigemoto and Isabel Llano **(2020)**.
*Synergistic action of metabotropic and ionotropic glutamate receptors in cerebellar molecular layer interneurons in vivo.* 

Raw data can be obtained upon request. 

Individual analsyis steps and scripts 
-----------

Experiments were done on **four different animals** in **five recording sessions**. All animal-specfic settings and analysis 
parameters are stored in the ```animalSettings.py``` file. The recording to be analyzed can be specified in the header of the scripts 
below or given as an input arguments. Possible recordings are : ```animal#1, animal#3, animal#2, animal#4, animal#1_2``` , e.g. 
```python alignImages.py animal#1```. 

### Alignment of images : ```alignImages.py```

This scripts uses the average images from before and after drug application, as well as from before drug application and the recording 
at 820 nm laser stimulation wavelength. Image pairs are aligned using the openCV 
[```warpAffine```](https://docs.opencv.org/2.4/modules/imgproc/doc/geometric_transformations.html?highlight=warpaffine) function. 
Alignment output is presented in a figure, e.g. [ImageAlignment_animal#2](figureOutput/ImageAlignment_animal%232.pdf). 

```python
python alignImages.py
```
### ROI matching : ```matchROIs.py```

Using the aligned fluorescent images, individual ROIs are identified to represent the same cell body before and 
 after drug application based on their fractional overlap.
ROI matching output is presented in a figure, e.g.  [ROImatching_animal#2](figureOutput/ROImatching_animal%232.pdf). 

```python
python matchROIs.py
```

### Obtain timing information of the 2-photon imaging stacks : ```readExtractTimeStamps.py```

The timing information of individual fluorescent frames recorded with ScanImage are read and 
saved in pickle files for further analysis. 

```python
python readExtractTimeStamps.py
```

### Extracting locomotion fluorescence : ```analyzeFluoTraces.py```

This script contains the core fluorescent trace analysis. Baseline fluorescence is determined based on the pre-motorization
activity. Delta F over F baseline is calculated per ROI. 
The fluorescence during the forced locomotion period is extracted  before and after drug application.
The decay in fluorescence across locomotion sessions is removed based on the recordings before drug application. See 
[FluorescenceTraces_animal#2](figureOutput/FluorescenceTraces_animal%232.pdf) for an example output figure. 

Data across recordings is collected and summarized using ```analyzeAlexaFluoTraces.py```. See 
[SummaryFluorescenceChanges](figureOutput/SummaryFluorescenceChanges.pdf) figure. 

```python
python analyzeFluoTraces.py
python analyzeAlexaFluoTraces.py 
```

### Comparing fluorescence dynamics during 910 nm and 820 nm runs : ```analyzeFluoTraces820nm.py```

In this script the locomotion fluorescence dynamics is compared between imaging runs at 910 nm and 820 nm excitation laser 
wavelength. GCaMP fluorescence is not sensitive to calcium at 820 nm which allows to assess movement artifacts at this 
excitation wavelength. 
See [FluorescenceTraces820_animal#2](figureOutput/FluorescenceTraces820_animal%232.pdf) for an example output figure. 

Data across recordings is collected and summarized using ```summarize820Fluctuations.py```. See 
[Summary820Fluctuations](figureOutput/Summary820Fluctuations.pdf) figure. 

```python
python analyzeFluoTraces820nm.py
python summarize820Fluctuations.py
```

### Analyzing the fluorescence increase from Alexa 594 delivered with drug  : ```analyzeAlexaFluoTraces.py```

In this script the raw red fluorescence from Alexa 594 is extracted from the calcium imaging data 
 as function of drug solution application time. Alexa 594 was added
  to the drug containing ACSF solution in order to assess drug access to the imaged area.  
See [Alexa594FluorescenceTraces_animal#2](figureOutput/Alexa594FluorescenceTraces_animal%232.pdf) for an example output figure. 

Data across recordings is collected and summarized using ```summarizeAlexaFluoChanges.py```. See 
[SummaryAlexa594FluorescenceChanges](figureOutput/SummaryAlexa594FluorescenceChanges_5Experiments.pdf) figure. 

```python
python analyzeAlexaFluoTraces.py
```

### Run scripts in a batch : ```runAllScripts.py```

All scripts can be run in a batch using the ```runAllScripts.py``` file. The recording(s) to be analyzed is (are) specified
 through a list  in that file. 

```python
python analyzeFluoTraces.py
```

### Generation of the publication figure : ```generatePublicationFigure.py```

The figure illustrating the forced locomotion experiment, depicting example data and summarizing the results obtained from all 
recordings is generated with that script. 

```python
python generatePublicationFigure.py
```

Requires
-----------
Standard python packages such as **numpy**, **scipy**, **pylab**, **time**, **os**,  **sys** and **matplotlib** are required.

License
-----------
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License, version 3
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://opensource.org/licenses/GPL-3.0>.

