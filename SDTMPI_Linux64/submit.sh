#PBS -q UCTlong
#PBS -l nodes=1:series400M:ppn=4
#PBS -N SDTMPI_Linux64
#PBS -e /home/brejnev/SDT_project/SDTMPI_Linux64/SDTMPI_Linux64.err
#PBS -o /home/brejnev/SDT_project/SDTMPI_Linux64/SDTMPI_Linux64.out

cd /home/brejnev/SDT_project/SDTMPI_Linux64

mpiexec -machinefile ${PBS_NODEFILE} /opt/exp_soft/python-2.7.2/lib/python2.7/site-packages/mpi4py/bin/python-mpi SDTMPI_Linux64.py test.fas muscle


