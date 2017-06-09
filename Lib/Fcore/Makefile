F2PY = f2py


all: f_rdf mat_ops


%: %.f90
	${F2PY} -c $< -m $@


clean:
	rm -rf *.mod
	rm -rf *.so
	rm -rf *.o
