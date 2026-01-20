# SMELDit
A semi-automated pipeline for liposome recognition, segregation, analysis and display

## License
This project is licensed under the MIT License.

## Liposome Recognition
- To generate liposome recognitions, use the `Save_Recognized_liposomes.m` and follow the instructions for setting the input directory, number of experiments and number of images per experiment to process
- The script will take the wide field input images and automatically detect the liposomes in them and create a subfolder to store each one as a cropped image with a unique ID and save its quantities of interest
- After the automated recognition process, the recognized liposomes are displayed as a 2d histogram. The user can select different metrics for each axis through two dropdown menus
- Due to their size, image datasets and liposome recognitions are not provided

## Working with an Existing Set of Liposome Recognitions
- To load a dataset that has already been preprocessed, use the `Load_Recognized_liposomes.m` and set the directory to the folder that contains the cropped liposome images. It will load in the data set and open up the standard user interface. 

## User defined ROIs and Data Set Selections

- The user is prompted by the GUI to define their own Regions of Interest (ROIs), these can be saved using the `Save ROI` button
- New ROIs can be created by clicking the `New ROI` button
- If any ROIs have been created previously, they will be available under the `Load ROI` button
- The subset of data present inside a ROI can be exported as its own separate data set with the `Export Data Selection` button
- To inspect these selection, the user can import them using the `Import Data Selection` button. This will also allow the user to draw new ROIs specific to the subselection should they want to.


## Image display by IMDISP
Finally the liposomes are also displayed as a montage of the cropped images. For this we used IMDISP, which was originally developed by Oliver Woodford on 10 Dec 2008 (Updated 09 Jun 2010) as a standalone package for MATLAB. We have included the code here for future compatibility reasons as it is no longer maintained by the original author. Oliver Woodford (2024). imdisp (https://www.mathworks.com/matlabcentral/fileexchange/22387-imdisp), MATLAB Central File Exchange. Retrieved November 6, 2024.

## Developer
Mats van Tongeren, with modifications by Federico Ramirez Gomez.

