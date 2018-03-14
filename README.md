A Matlab/Octave function `get_coupled_QRTwaves` for finding the
coupled ECG waves annotated using WFDB toolbox of PhysioNet. The main
application of the function is the mutual QT-RR interval dynamics in
the heart.

The examplar use is demonstrated in `example_LTST_DB.m` where the
PhysioNet Long-Term ST Database is annotated using `ecgpuwave` and the
coupled waves are found using the annotations. Do not run this example
as a whole, the DB contains 86 ~24h long ECG recordings! This would
take days to annotate and read the annotations. Use the file as a
template for other signals. The function `get_coupled_QRTwaves`
however is quite fast (~10 sec per ~24h long recording), although the
implementation is pretty straightforward.
