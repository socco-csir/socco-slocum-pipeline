# GEOMAR FB1-PO Matlab Slocum glider processing toolbox

- **Authors:** <br>
Gerd Krahmann

- **Software Title:**  <br>
GEOMAR FB1-PO Matlab Slocum glider processing toolbox

- **Short Description:** <br>
This Matlab toolbox contains code to process data from Teledyne Webb Research Slocum underwater gliders.

- **License:** <br>
https://git.geomar.de/open-source/geomar_glider_toolbox/-/blob/main/LICENSE

- **Requirements:** <br>
This is a toolbox for [Mathwork's Matlab](https://www.mathworks.com). It requires Mathwork's Signal Processing and Statistics and Machine Learning Toolboxes. It has been tested with Matlab version R2018a but should work with any
reasonably new Matlab version. To process raw data you also need executable files by the glider manufacturer Teledyne Webb Research.

- **Installation Instructions:** <br>
To install the toolbox, copy all files in the folder `glider` and in its subfolders to your computer. In Matlab add the path to the folder `glider/glider_toolbox` to Matlab's search path. You also need the [Gibbs-SeaWater Toolbox](https://teos-10.org/software.html). A [zip archive](https://git.geomar.de/open-source/geomar_glider_toolbox/-/blob/main/glider/ifm13_depl07_example.zip) contains an example dataset.
The processing requires the reformatting of the raw glider data into something readable by Matlab. A number of proprietary programs are available from Teledyne Webb Research. These programs need to be located in the folder `glider/glider_toolbox/proprietary`.

- **Usage Instructions:** <br>
Are described in detail in https://git.geomar.de/open-source/geomar_glider_toolbox/-/blob/main/glider/README.md
