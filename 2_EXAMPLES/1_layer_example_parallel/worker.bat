unzip data.zip

cd data
dir
python condor_single_run.py %1 %2

python zip_results.py %1
