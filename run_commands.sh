export TMPDIR=/local/mrahman/magician/tmp && export TMP=/local/mrahman/magician/tmp && export TEMP=/local/mrahman/magician/tmp

python run_magician.py /home/mrahman/magician/data/from_real_v3/sample_distributions.tsv \
  --cores 110 \
  --snake_flags " \
  --directory /local/mrahman/magician/results/from_real_v3 \
  --rerun-incomplete"