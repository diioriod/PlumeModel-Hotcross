#!/bin/bash
#SBATCH --job-name=Info-Hotcross_s2ALT
#SBATCH --partition=aquari_p
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 
#SBATCH --mem=10gb
#SBATCH --time=15-00:00:00
#SBATCH --output=%x.%j.out

module load PGI/17.9
module load netCDF/4.1.3-PGI-17.9

cd /home/daniela/ForGitHub

pgf90 -mp -r8 -Mprof -Mextend -fast -Kieee -Minline -Mdclchk -O3 -I/apps/eb/netCDF/4.1.3-PGI-17.9/include put_block.f put_block1.f put_block2.f modules80-revised2014.f hotcross80-revised2014.f source80.f smagorinsky_split_mixing-revised2014.f potential_density80.f initialize_profiles80.f hs3crt_composite80.f mixing80-revised2014.f body_force80.f sponge80.f mpdata-hotcross80.f write_snapshota.f netcdf_augmentation.f view_3d.f view_2d.f stat_2d.f stat_3d.f -L/apps/eb/netCDF/4.1.3-PGI-17.9/lib -lnetcdf -lnetcdff -L/usr/bin -lcurl -o Hotcross-SAPELO2 -mcmodel=medium

echo
echo "Job ID: $SLURM_JOB_ID"
echo "QUEUE: $SLURM_JOB_PARTITION"
echo "Cores: $SLURM_CPUS_PER_TASK"
echo


export LD_LIBRARY_PATH=/apps/eb/netCDF/4.1.3-PGI-17.9/lib:${LD_LIBRARY_PATH}
export OMP_NUM_THREADS=8
cp site_stratification_fromobs.inc /scratch/daniela/DanteExp1/
cp cfl.out /scratch/daniela/DanteExp1/
cp temp.diagnostics /scratch/daniela/DanteExp1/
cp salt.diagnostics /scratch/daniela/DanteExp1/
cp tracer.diagnostics /scratch/daniela/DanteExp1/
cp Hotcross-SAPELO2 /scratch/daniela/DanteExp1/
cp DanteExp1.case /scratch/daniela/DanteExp1/
cd /scratch/daniela/DanteExp1


./Hotcross-SAPELO2 > DanteExp1.out
