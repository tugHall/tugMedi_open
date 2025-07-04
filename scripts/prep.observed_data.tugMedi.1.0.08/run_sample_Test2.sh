python ./src/prepare_observed_data.py -i ../../inst/extdata/tugMedi.1.0/Test2/Input.TCGA_7903/Samples/Processes/TCGA.LUAD.55-7903.01.txt  -c ./samples/sample_Test2/config/config_vaf.conf -m get-vaf -o ./samples/sample_Test2/output_vaf

python ./src/prepare_observed_data.py -i ./samples/sample_Test2/output_vaf/TCGA.LUAD.55-7903.01_vaf.txt  -c ./samples/sample_Test2/config/config_driver.conf -m get-driver-genes -o ./samples/sample_Test2/output_driver

python ./src/prepare_observed_data.py -i ./samples/sample_Test2/output_vaf/TCGA.LUAD.55-7903.01_vaf.txt  -c ./samples/sample_Test2/config/config_passenger.conf -m get-passenger-genes -o ./samples/sample_Test2/output_passenger
# TCGA.LUAD.55-7903.03_Rother.chr1-3_PASS.txt: sample_Test2/output_passenger/TCGA.LUAD.55-7903.01_vaf.Rother.txt

python ./src/prepare_observed_data.py -i ./samples/sample_Test2/output_vaf/TCGA.LUAD.55-7903.01_vaf.txt  -c ./samples/sample_Test2/config/config_passenger2.conf -m get-passenger-genes -o ./samples/sample_Test2/output_passenger2
# TCGA.LUAD.55-7903.04_Rother.gname_vaf0.43.txt: sample_Test2/output_passenger2/TCGA.LUAD.55-7903.01_vaf.Rother.txt

python ./src/prepare_observed_data.py -i ../../inst/extdata/tugMedi.1.0/Test2/Input.TCGA_7903/Samples/Processes/TCGA.LUAD.55-7903.01.txt  -c ./samples/sample_Test2/config/config.conf -o ./samples/sample_Test2/output_all
# samples.Rint.txt = sample_Test2/output_all/TCGA.LUAD.55-7903.01_tumor_specific_vaf.Rint.txt
# sample.Rother.txt = sample_Test2/output_all/TCGA.LUAD.55-7903.01_tumor_specific_vaf.Rint.txt
