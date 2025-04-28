

bash mcr.input.bash
bash mcr.remove.bash


Rscript abc.a.01_make_list.R abc.01_config.R
bash abc.b_run_parallel.sh
Rscript abc.c_make_data.R abc.01_config.R
Rscript abc.d_exec_calc.R abc.01_config.R


Rscript --vanilla --slave ./test2.3_eval.R & 


Rscript abc.a.02_make_list.R abc.02_config.R
bash abc.b_run_parallel.sh
Rscript abc.c_make_data.R abc.02_config.R
Rscript abc.d_exec_calc.R abc.02_config.R


Rscript --vanilla --slave ./test2.4_eval.R 


