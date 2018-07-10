# Process (sh)MOLLI data

Ben George

University of Oxford

MATLAB and python combination to process (sh)MOLLI MRI data. python fitting by default, will fallback to (slower) MATLAB curve fitting if python not available.

## Install
Copy files to MATLAB path. Both files must be in the same directory.

*Usage 1:*

Create a MATLAB struct with the processed data.

Flags to optionally plot the data at the end and to allow user interaction to specify ROI

`T1MAP = process_sh_molli_series(directory,(plot),(user selection));`

*Usage 2:*

Lauch the viewer with some already processed data

`process_sh_molli_series(T1MAP);`

