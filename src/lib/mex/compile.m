mex -c ga.c
mex -c utils.c
mex -c fitness.c
mex gaoptim.c ga.obj fitness.obj utils.obj
mex getcrossdist.c
mex getnetsinr.c