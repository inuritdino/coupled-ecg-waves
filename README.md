A Matlab/Octave function `get_coupled_QRTwaves` for finding the
coupled ECG waves annotated using WFDB toolbox of PhysioNet. Besides
the coupled wave outcome, the function accounts for the contiguous
regions of the coupled annotated waves, giving the interruption
locations and sizes of the interruptions.

The example code is demonstrated in `example_QT_DB.m` where the
coupled waves are found using the annotations of the PhysioNet QT
Database.

Currently the main application of the function is the mutual QT-RR
interval dynamics of the heart.

