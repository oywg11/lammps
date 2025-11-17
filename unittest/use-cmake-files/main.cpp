
#ifdef USE_MPI
#include <mpi.h>
#else
#define MPI_COMM_WORLD 0
#endif
#include <lammps/lammps.h>
#include <lammps/input.h>
#include <lammps/atom.h>
#include <exception>
#include <iostream>

int main(int argc, char **argv)
{
        MPI_Init(&argc, &argv);
        try {
                auto *lmp = new LAMMPS_NS::LAMMPS(argc, argv, MPI_COMM_WORLD);
                lmp->input->file();
                delete lmp;
        } catch (std::exception &) {
        }
        return 0;
}

