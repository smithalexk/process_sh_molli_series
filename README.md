# Process (sh)MOLLI data

Ben George

University of Oxford

MATLAB and python combination to process (sh)MOLLI MRI data. python fitting by default, will fallback to (slower) MATLAB curve fitting if python not available.

## Instructions:

1. Place both files in the same directory

2. In order to use the python fitting (which is faster), when you start Matlab, run the command

`> pyversion`

This should find your python installation and let Matlab use it. If this command returns nothing, you might need to setup python.

If you’re using Wolfgang and pyversion does’t do anything, try running:

`> pyversion 'C:\Python35\python.exe'`

3. Call the script with the command:

`> T1MAP = process_sh_molli_series;`

A directory popup window will be displayed, allowing you the select the path to the files.

The files will be read in and you will be prompted to select from the TiggerTime or InversionTime tags (for backwards compatibility). One option will either be blank or have all the same number. You want to select the option with varying values. Normally TriggerTime is the correct option.

4. A window will pop up allowing you to select the ROI to process. Click multiple points to trace out the region. Double click on the start point to finish your selection. Pressing ESC will process the whole image. This can take a while.

5. The magic happens using parallel processing and python fitting to speed things up. There is a fancy process indicator which will try and predict when it will finish.

6. Once the processing is complete, a 2x2 window will pop up so you can view the data. Top left is the input image. Top right is the fitting error. Bottom left is the processed data. Bottom right is the data and curve fit. Select different points to see the fitted data. Click outside any subplot to finish.

7. The returned value ‘T1MAP’ (or whatever) is a Matlab struct which has loads of information in, including: directory of the files used, input data, output data, curve fits, ROI mask and others.

8. You can call up the interactive ‘review’ window from any processed data using

`> proces_sh_molli_series(T1MAP);`

9. You can pass the directory path straight to script by running

`> T1MAP = process_sh_molli_series(path_to_dcm_files);`

10. To make things easier, you can forego the final interactive review window by running

`> T1MAP = process_sh_molli_series(path_to_dcm_files, false);`

11. You can even skip the pre-processing stages by running

`> T1MAP = process_sh_molli_series(path_to_dcm_files, false, false);`

This will assume that TriggerTime is the correct tag and process the whole image.

12. It wouldn’t be too hard to write a small wrapper to chug through loads of data in the background and save all the results. Assuming a series of subdirectories each containing DICOM files, that could would look something like:
```
D = dir;
for i = 1:length(D)
	if D(i).isdir
		T1MAP = process_sh_molli_series([D(i).path filesep D(i).name, false, false);
		save([pwd filesep D(i).name '_T1MAP.mat'], T1MAP);
	end
end
```

13. At the moment, the data is saved to a Matlab struct. You can access the output image from `T1MAP.output_LL_corrected` to display, save etc. It doesn’t write the saved data to DICOM or anything. This is the next thing to do (Dan McGowan was asking for this functionality). To save the final image as a PNG, run:

```
> imwrite(uint16(T1MAP.output_LL_corrected), 'T1MAP.png');
```

14. The script assumes that the directory you give it contains only the image files in DICOM format. It will likely break if this is not the case.

15. In theory, the script should work with both MOLLI or shMOLLI data and be capable of processing 3D shMOLLI volumes (if that is a thing).