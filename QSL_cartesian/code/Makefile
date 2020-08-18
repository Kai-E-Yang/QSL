FC   = gfortran
Flag = -fopenmp -ffree-form -ffree-line-length-none -O3 -fcheck=bounds

files = mod_param.f90 mod_operator.f90 mod_io.f90 mod_solver.f90 main.f90

QSL: 	$(files) mod_io.mod mod_param.mod mod_operator.mod mod_solver.mod
	$(FC) ${Flag} -o $@ $(files)

clean:
	rm -f QSL *.mod *.o

mod_io.mod: mod_param.f90 mod_operator.f90 mod_io.f90 
	$(FC) -c $^
mod_operator.mod: mod_param.f90 mod_operator.f90
	$(FC) -c $^
mod_solver.mod: mod_param.f90 mod_operator.f90 mod_solver.f90 
	$(FC) -c $^
mod_param.mod: mod_param.f90
	$(FC) -c $^

%.o: %.f90
	$(FC) -c $<