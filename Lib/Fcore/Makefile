F2PY = f2py
FVENDOR = gnu95


all: f_rdf f_rdf2 mat_ops


%: %.f90
	${F2PY} -c $< --fcompiler=${FVENDOR} -m $@


clean:
	rm -rf *.mod
	rm -rf *.so
	rm -rf *.o
