
Analysis scripts of calcium imaging in locomoting mice 
==============================

The repository contains custom written Python scripts used to extract calcium transients during forced locomotion. 

The here published python scripts implement the steps and calculation to obtain **locomotion fluorescence** as depicted in Fig. 3 
 of the publication below. 


For more details, please refer to :

Jinn Bao, Michael Graupner, Guadalupe Astorga, Thibault Collin, Abdelali Jalil, Dwi Wahyu Indriati, Jonathan Bradley, 
Ryuichi Shigemoto and Isabel Llano **(2020)**.
*Synergistic action of metabotropic and ionotropic glutamate receptors in cerebellar molecular layer interneurons in vivo.* 


Individual analsyis steps and scripts 
-----------

Experiments were done on **four different animals** in **five different recording sessions**. All animal specfic settings and analysis 
parameters are stored in the ```animalSettings.py``` file. The recording to be analyzed can be specified in the header of the scripts
or given as an input arguments. Possible recordings are : ```animal#1, animal#3, animal#2, animal#4, animal#1_2``` . 

### Alignment of images : ```alignImages.py```

This scripts uses the average images from before and after drug application, as well as from before drug application and the recording 
at 820 nm laser stimulation wavelength. Image pairs are aligned using the openCV ```warpAffine``` function. Alignment ouptut is 
presented in a figure (e.g. figureOutput/ImageAlignment_animal#2.pdf). 

```python
python alignImages.py
```
### ROI matching : ```matchROIs.py```

Using the algined fluorescent images, individual ROIs are indentified to represent the same cell body based on their overlap fraction. 

```python
python matchROIs.py
```

### Obtain timing information of the 2-photon imaging stacks : ```readExtractTimeStamps.py```

The timing information of individal fluorescent frames recorded with ScanImage are read and saved in pickle files for further analysis. 

```python
python readExtractTimeStamps.py
```

### Extracting locomotion fluorescence : ```analyzeFluoTraces.py```

This script contains the core fluorescent trace analysis. Baseline fluorescence is determined based on the pre-motorization
baseline. Delta F over F baseline is calculated per ROI. The fluorescence during the forced locomotion period is extracted. 
The decay in fluorescence across locomotion sessions is removed based on the recordings before drug application. See 
figureOutput/FluorescenceTraces_animal#2.pdf for an example output figure. 

```python
python analyzeFluoTraces.py
```

### Comparing fluorescence dynamics during 910 nm and 820 nm runs : ```analyzeFluoTraces820nm.py```

In this script the locomotion fluorescence dynamics is compared between imaging runs at 910 nm and 820 nm laser wavelength 
stimulation. GCaMP fluorescence is not sensitive to calcium at 820 nm which allows to assess movement artifacts at this 
stimulation wavelength. 
See [FluorescenceTraces820_animal#2](figureOutput/FluorescenceTraces820_animal%232.pdf) for an example output figure. 

```python
python analyzeFluoTraces820nm.py
```

### Analyzing the fluorescence increase from Alexa 594 delivered with drug  : ```analyzeAlexaFluoTraces.py```

In this script the raw red fluorescence from Alexa 594 is extracted as function of drug solution application time. Alexa 594 was added
  to the drug containing ACSF solution in order to assess drug access to the imaged area.  
See [Alexa594FluorescenceTraces_animal#2](figureOutput/Alexa594FluorescenceTraces_animal%232.pdf) for an example output figure. 

```python
python analyzeAlexaFluoTraces.py
```

### Run scripts in a batch : ```runAllScripts.py```

All scripts can be run in a batch using the ```runAllScripts.py``` file. The recording(s) to be analyzed is specified through a list 
in that file. 

```python
python analyzeFluoTraces.py
```

Requires
-----------
Standard python packages such as **numpy**, **scipy**, **pylab**, **time**, **os**,  **sys** and **matplotlib** are required.

License
-----------
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

