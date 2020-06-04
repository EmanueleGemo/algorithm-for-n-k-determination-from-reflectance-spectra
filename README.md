# nk-determination-algorithm
algorithm that takes a FP cavity reflection spectrum and extract the n &amp; k cavity properties by comparison to literature data
MAIN file calls the functions as intended
materials parameters must be stored in separated .mat files, and must contain a nx3 array variable of the same name, with the wavelength listed on the first column, the refractive index on the 2nd column, and the extinction coefficient on the third column.
reflectance spectra must be similarly stored in .mat files, under a single nx2 variable containing the wavelength on the 1st column, and the reflectance (abs) on the 2nd column.
each .material and .file field must contain a cell array, with each cell containing a char string of the desired material or reflectance .mat file.
