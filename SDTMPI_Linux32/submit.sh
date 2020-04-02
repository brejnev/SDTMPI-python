#PBS -q UCTlong
#PBS -l nodes=1:series400M:ppn=4
#PBS -N SDTMPI_Linux32
#PBS -e /home/brejnev/SDT_project/SDT_project_uplaod/SDTMPI_Linux32/SDTMPI_Linux32.err
#PBS -o /home/brejnev/SDT_project/SDT_project_uplaod/SDTMPI_Linux32/SDTMPI_Linux32.out

cd /home/brejnev/SDT_project/SDT_project_uplaod/SDTMPI_Linux32

mpiexec -machinefile ${PBS_NODEFILE} /opt/exp_soft/python-2.7.2/lib/python2.7/site-packages/mpi4py/bin/python-mpi SDTMPI_Linux32.py test.fas muscle


