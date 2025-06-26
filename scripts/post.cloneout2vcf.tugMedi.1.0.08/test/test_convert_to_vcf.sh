
# No time specified, all time points will be used
echo "[INFO] start test both mode"

mkdir -p TEST_both/testdata/

./cloneout2vcf.py \
    --time max \
    --cloneout ./testdata/0003/Output/cloneout.txt \
    --pointMutations ./testdata/0003/Output/Mutations/pointMutations.txt \
    --VAF ./testdata/0003/Output/VAF/VAF.txt \
    --fasta ./testdata/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --maf ./testdata/TCGA-AZ-6608-01A-11D-1835-10_may.maf \
    --seed \
    --output TEST_both

# Copy generated VCF files as original files
for vcf_file in TEST_both/output_*.vcf; do
    if [ -f "$vcf_file" ]; then
        cp "$vcf_file" "${vcf_file}.ori"
        echo "[INFO] Backup created: ${vcf_file}.ori"
    fi
done
