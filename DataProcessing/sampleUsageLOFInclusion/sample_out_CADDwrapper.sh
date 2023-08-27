#!/bin/bash -e 
jid0=($(sbatch chromosomeRename.sh))
echo "jid0 ${jid0[-1]}" >> slurm_ids
jid1=($(sbatch --dependency=afterok:${jid0[-1]} sample_out_cadAnno.sh))
echo "jid1 ${jid1[-1]}" >> slurm_ids
jid2=($(sbatch  --dependency=afterok:${jid0[-1]} sample_out_annovarAnno.sh))
echo "jid2 ${jid2[-1]}" >> slurm_ids
jid3=($(sbatch --dependency=afterok:${jid1[-1]} sample_out_LOFs.sh))
echo "jid3 ${jid3[-1]}" >> slurm_ids
jid4=($(sbatch --dependency=afterok:${jid3[-1]} 12.37_filtering.sh))
jid5=($(sbatch --dependency=afterok:${jid3[-1]} 20_filtering.sh))
