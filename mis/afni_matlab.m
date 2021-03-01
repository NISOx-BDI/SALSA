
SubID='A00028150'
SesID='DS2'

runafni=['sh /users/nichols/scf915/bin/FILM2/mis/run_Null3dREMLfit.sh ROCKLAND 1.4 404 ' SubID ' ' SesID ' 5 0 0 /well/nichols/users/scf915/ROCKLAND/R.PW/ROCKLAND_1400_404_3dREMLfit_AR-1_MA-1_FWHM5_poly_gsr0_aroma0'];
disp(runafni)
flag = system(runafni)

flag
