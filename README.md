# dendra2_Analysis
Code to analyse 3D multichannel (XYZC) dendra2 image series in MATLAB.  


## Data Format
- All image files for an experiment (recordings across multiple days, one XYZC format file per recording) should be saved in a single folder.
    Files for different experiments should be saved in separate folders.

- Folder names should be specified in a batch file, provided as a columnwise list in a spreadsheet (see "batch_example.xlsx").

- Images files should be saved in OME-TIFF format (1,2). If the raw images are saved in another format, try opening in Fiji software
    (https://imagej.net/software/fiji/) (3) and saving as an OME-TIFF using the Bio-Formats Exporter (https://imagej.net/formats/bio-formats) (4).

- Image files should follow a consistent naming convention, with the experiment day specified after a space at the end of the filename 
    immediately before the ".ome.tif" suffix (e.g., "Experiment 1\_day 1.ome.tif").

- On the first day of the experiment ('day 1') two recordings are performed - one before dendra2 photoconversion, and one after photoconversion. 
    The pre-photoconversion file should be named "day 1" (e.g., "Experiment 1\_day 1.ome.tif"),
    and the post photoconversion should be named "day 1.1" (e.g., "Experiment 1\_day 1.1.ome.tif").


## Requirements
- "Bioformats Image Toolbox" for MATLAB
   (see https://www.mathworks.com/matlabcentral/fileexchange/129249-bioformats-image-toolbox) (5)

- "Efficient subpixel image registration by cross-correlation" toolbox for MATLAB 
   (see https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation) (6)


## Installation
Download this repository and add the relevant directories (including Requirements) to the MATLAB path.


## Usage
All functions can be run from "dendra2_Master.m", with an input of the full batch file name (as a string). E.g.:

```matlab
dendra2_Master('C:\dendra2_Analysis\batch_example.xlsx')
```

## Licence and Acknowledgement
The code is provided "As Is" under the BSD-3-Clause license.

If you use the code in academic research, please cite this GitHub repository: https://github.com/WittonLab/dendra2_Analysis/ 


## References
(1) Goldberg I, et al. (2005). Genome Biol. 6:R47. DOI: 10.1186/gb-2005-6-5-r47 <br>
(2) Besson S, et al. (2019). Lecture Notes in Computer Science, vol 11435. DOI: 10.1007/978-3-030-23937-4_1 <br>
(3) Schindelin et al. (2012). doi:10.1038/nmeth.2019 <br>
(4) Linkert M, Et al. (2010). J. Cell Biol. 189(5), 777-782. DOI: 10.1083/jcb.201004104 <br>
(5) Tay JW. (2026). Bioformats Image Toolbox (https://github.com/Biofrontiers-ALMC/bioformats-matlab/releases/tag/v1.2.3), GitHub. Retrieved August 8, 2026. <br>
(6) Guizar-Sicairos M, et al. (2008). Opt. Lett. 33, 156-158. DOI: 10.1364/ol.33.000156
