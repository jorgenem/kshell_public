from typing import List

def betzy(
    shell_filename: str,
    sigma2_project_name: str,
    sigma2_n_days: int,
    sigma2_n_hours: int,
    sigma2_n_minutes: int,
    n_nodes: int,
    n_tasks_per_node: int,
    n_cpus_per_task,
    sigma2_user_email: str,
    type_of_betzy_job: str,
    omp_num_threads: int
    ) -> str:
    """
    Generate SLURM commands for Betzy.
    """
    job_commands = '#!/bin/bash \n'
    job_commands += f'#SBATCH --job-name={shell_filename[:-3]} \n'
    job_commands += f'#SBATCH --account={sigma2_project_name} \n'
    job_commands += '## Syntax is d-hh:mm:ss \n'
    job_commands += f'#SBATCH --time={sigma2_n_days}-{sigma2_n_hours:02d}:{sigma2_n_minutes:02d}:00 \n'
    job_commands += f'#SBATCH --nodes={n_nodes}\n'
    job_commands += f'#SBATCH --ntasks-per-node={n_tasks_per_node} \n'
    job_commands += f'#SBATCH --cpus-per-task={n_cpus_per_task} \n'
    job_commands += '#SBATCH --mail-type=ALL \n'
    job_commands += f'#SBATCH --mail-user={sigma2_user_email} \n'
    
    if type_of_betzy_job != "normal":
        job_commands += f'#SBATCH --qos={type_of_betzy_job} \n'
    
    job_commands += 'module --quiet purge  \n'
    job_commands += 'module load intel/2020b \n'
    job_commands += 'module load Python/3.8.6-GCCcore-10.2.0 \n'
    job_commands += 'set -o errexit  \n'
    job_commands += 'set -o nounset \n'
    
    if omp_num_threads is not None:
        job_commands += f'export OMP_NUM_THREADS={omp_num_threads} \n'
    
    return job_commands

def fram(
    shell_filename: str,
    sigma2_project_name: str,
    sigma2_n_days: int,
    sigma2_n_hours: int,
    sigma2_n_minutes: int,
    n_nodes: int,
    n_tasks_per_node: int,
    n_cpus_per_task,
    sigma2_user_email: str,
    type_of_fram_job: str,
    ) -> str:
    """
    Generate SLURM commands for Fram.
    """
    job_commands = '#!/bin/bash \n'
    job_commands += f'#SBATCH --job-name={shell_filename[:-3]} \n'
    job_commands += f'#SBATCH --account={sigma2_project_name} \n'
    job_commands += '## Syntax is d-hh:mm:ss \n'
    job_commands += f'#SBATCH --time={sigma2_n_days}-{sigma2_n_hours:02d}:{sigma2_n_minutes:02d}:00 \n'
    job_commands += f'#SBATCH --nodes={n_nodes}\n'
    job_commands += f'#SBATCH --ntasks-per-node={n_tasks_per_node} \n'
    job_commands += f'#SBATCH --cpus-per-task={n_cpus_per_task} \n'
    job_commands += '#SBATCH --mail-type=ALL \n'
    job_commands += f'#SBATCH --mail-user={sigma2_user_email} \n'
    
    if type_of_fram_job != "normal":
        job_commands += f'#SBATCH --qos={type_of_fram_job} \n'
    
    job_commands += 'module --quiet purge  \n'
    # job_commands += 'module load foss/2017a \n'
    job_commands += 'module load intel/2020b \n'
    job_commands += 'module load Python/3.8.6-GCCcore-10.2.0 \n'
    job_commands += 'set -o errexit  \n'
    job_commands += 'set -o nounset \n'

    return job_commands

def no_scheduler() -> str:
    job_commands = '#!/bin/sh \n'
    job_commands += '# export OMP_STACKSIZE=1g\n'
    job_commands += 'export GFORTRAN_UNBUFFERED_PRECONNECTED=y\n'
    job_commands += '# ulimit -s unlimited\n\n'
    
    return job_commands

def pjm_default(n_nodes: int) -> str:
    job_commands = '#!/bin/sh \n'
    job_commands += '#PJM -L "rscgrp=debug"\n'
    job_commands += '#PJM -L "node=' + str(n_nodes) + '"\n'
    job_commands += '# #PJM -L "elapse=24:00:00"\n\n'
    # job_commands += 'cd ' + os.getcwd() +'\n\n'

    return job_commands

def cx400(n_nodes: int) -> str:
    job_commands = "#!/bin/sh\n"
    job_commands += '#PJM -L "rscgrp=XXXXXXXX"\n'
    job_commands += f'#PJM -L "vnode={n_nodes}"\n'
    job_commands += '#PJM -L "vnode-core=28"\n'
    job_commands += '# #PJM --mpi "rank-map-bynode"\n'
    job_commands += '#PJM -P "vn-policy=abs-unpack"\n'
    job_commands += '#PJM -L "elapse=01:00:00"\n'
    job_commands += '#\n\n'
    job_commands += 'source /center/local/apl/cx/intel/composerxe/bin/compilervars.sh intel64\n'
    job_commands += 'source /center/local/apl/cx/intel/impi/4.1.1.036/bin64/mpivars.sh\n'
    job_commands += 'source /center/local/apl/cx/intel/mkl/bin/mklvars.sh intel64\n\n'
    job_commands += 'export I_MPI_PIN_DOMAIN=omp\n'
    job_commands += '# export OMP_NUM_THREADS=28\n'
    job_commands += 'export I_MPI_HYDRA_BOOTSTRAP=rsh\n'
    job_commands += 'export I_MPI_HYDRA_BOOTSTRAP_EXEC=/bin/pjrsh\n'
    job_commands += r'export I_MPI_HYDRA_HOST_FILE=${PJM_O_NODEINF}' + "\n"
    job_commands += 'export FORT90L=-Wl,-Lu\n'
    # + 'cd ' + os.getcwd() +'\n\n' \

    return job_commands

def ofp_flat(n_nodes: int) -> str:
    job_commands = '#!/bin/sh \n'
    job_commands += '#PJM -L "rscgrp=debug-flat"\n'
    job_commands += '#PJM -L "node=' + str(n_nodes) + '"\n'
    job_commands += '# #PJM --mpi "proc=' + str(n_nodes) + '"\n'
    job_commands += '#PJM --omp thread=272\n'
    job_commands += '#PJM -L "elapse=00:30:00"\n'
    job_commands += '#PJM -g XXXXXX\n'

    return job_commands

def ofp(n_nodes: int) -> str:
    job_commands = '#!/bin/sh \n'
    job_commands += '#PJM -L "rscgrp=debug-cache"\n'
    job_commands += '#PJM -L "node=' + str(n_nodes) + '"\n'
    job_commands += '#PJM --omp thread=272\n'
    job_commands += '#PJM -L "elapse=00:30:00"\n'
    job_commands += '#PJM -g XXXXXXXX\n'

    return job_commands

def k_large(
    stgin_filenames: List,
    stgout_filenames: List,
    n_nodes: int
    ) -> str:
    outstg = '#PJM --stgin "'
    for fn in stgin_filenames: outstg += './' + fn + ' '
    outstg += './"\n'
    outstg += '#PJM --stgout "'
    for fn in stgout_filenames: outstg += './' + fn + ' '
    outstg += './"\n\n'
    job_commands = '#!/bin/sh \n'
    job_commands += '#PJM -L "rscgrp=large"\n'
    job_commands += '#PJM -L "node=' + str(n_nodes) + '"\n'
    job_commands += '#PJM -L "elapse=06:00:00"\n'
    job_commands += outstg
    job_commands += '. /work/system/Env_base\n\n'
    job_commands += 'lfs setstripe -s 100m -c 12 .\n\n'

    return job_commands

def k_small(
    stgin_filenames: List,
    stgout_filenames: List,
    n_nodes: int
    ) -> str:
    outstg = '#PJM --stgin "'
    for fn in stgin_filenames: outstg += './'+fn+' '
    outstg += './"\n'
    outstg += '#PJM --stgout "'
    for fn in stgout_filenames: outstg += './'+fn+' '
    outstg += './"\n\n'
    job_commands = '#!/bin/sh \n'
    job_commands += '#PJM -L "rscgrp=small"\n'
    job_commands += '#PJM -L "node=' + str(n_nodes) + '"\n'
    job_commands += '#PJM -L "elapse=00:30:00"\n'
    job_commands += outstg
    job_commands += '. /work/system/Env_base\n\n'
    job_commands += 'lfs setstripe -s 100m -c 12 .\n\n'

    return job_commands

def k_micro(n_nodes: int) -> str:
    job_commands = '#!/bin/sh \n'
    job_commands += '#PJM -L "rscgrp=micro"\n'
    job_commands += f'#PJM -L "node={str(n_nodes)}\n'
    job_commands += '#PJM -L "elapse=00:30:00"\n'
    job_commands += '#PJM -g "XXXXXXXX"\n\n'
    job_commands += '. /work/system/Env_base\n\n'
    # + 'cd ' + os.getcwd() +'\n\n' \

    return job_commands

def coma(shell_filename: str, n_nodes: int) -> str:
    job_commands = '#!/bin/sh \n'
    job_commands += '#SBATCH -J ' + shell_filename[:-3] + '\n'
    job_commands += '#SBATCH -p mixed\n'
    job_commands += '#SBATCH -N ' + str(n_nodes) + '\n'
    job_commands += '#SBATCH -n ' + str(n_nodes) + '\n'
    job_commands += '# #SBATCH -t 01:00:00\n'
    job_commands += '#SBATCH --cpus-per-task=20\n'
    job_commands += '#SBATCH -o stdout\n'
    job_commands += '#SBATCH -e stderr\n\n'
    job_commands += 'export OMP_NUM_THREADS=20\n\n'
    job_commands += 'module load mkl intel intelmpi \n'
    # + 'cd ' + os.getcwd() +'\n\n' \
    # cd $SLURM_SUBMIT_DIR
    # export OMP_NUM_THREADS=16

    return job_commands