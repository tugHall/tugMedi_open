python ./src/prepare_observed_data.py -i ../../inst/extdata/tugMedi.1.0/Test1/Input/Samples/Processes/samples.02.txt  -c ./samples/sample_Test1/config/config.conf -m get-tumor-specific -o ./samples/sample_Test1/output_tumor_specific

python ./src/prepare_observed_data.py -i ./samples/sample_Test1/output_tumor_specific/samples.02_tumor_specific.txt  -c ./samples/sample_Test1/config/config.conf -m get-vaf -o ./samples/sample_Test1/output_vaf

python ./src/prepare_observed_data.py -i ./samples/sample_Test1/output_vaf/samples.02_tumor_specific_vaf.txt  -c ./samples/sample_Test1/config/config2.conf -m get-driver-genes -o ./samples/sample_Test1/output_driver

python ./src/prepare_observed_data.py -i ./samples/sample_Test1/output_vaf/samples.02_tumor_specific_vaf.txt  -c ./samples/sample_Test1/config/config2.conf -m get-passenger-genes -o ./samples/sample_Test1/output_passenger
# samples.Rother.04_chr1.2_PASS.txt = sample_Test1/output_passenger/samples.02_tumor_specific_vaf.Rother.txt

python ./src/prepare_observed_data.py -i ../../inst/extdata/tugMedi.1.0/Test1/Input/Samples/Processes/samples.02.txt  -c ./samples/sample_Test1/config/config.conf  -o ./samples/sample_Test1/output_all
# samples.Rint.txt = samples.Rint.05_APC1.txt = sample_Test1/output_all/samples.02_tumor_specific_vaf.Rint.txt
# samples.Rother.txt = samples.Rother.05_gname_vaf0.5.txt = sample_Test1/output_all/samples.02_tumor_specific_vaf.Rother.txt
