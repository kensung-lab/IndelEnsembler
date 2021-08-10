manta_sv=$1
smoove_sv=$2
surv_sv=$3
transurv_sv=$4
out_dir=$5
max_is=$6

python=python2
ref_genome=~/references/tair10/TAIR10_chr_all.fa
repeats=~/references/tair10/TAIR10_chr_all.trf.txt


# DELETIONS

cat $surv_sv | $python surv2sv.py | grep DEL | cut -d" " -f-8 > $out_dir/DEL.sv

bcftools view $manta_sv | $python manta2sv.py | grep DEL > $out_dir/manta-DEL.sv
./compare-del $out_dir/manta-DEL.sv $out_dir/DEL.sv $repeats $ref_genome $max_is | grep NONE | cut -d" " -f1 > $out_dir/manta-only.ids
grep -w -f $out_dir/manta-only.ids $out_dir/manta-DEL.sv >> $out_dir/DEL.sv
rm $out_dir/manta-DEL.sv
rm $out_dir/manta-only.ids

bcftools view $smoove_sv | $python lumpy2sv.py | grep DEL > $out_dir/smoove-DEL.sv
./compare-del $out_dir/smoove-DEL.sv $out_dir/DEL.sv $repeats $ref_genome $max_is | grep NONE | cut -d" " -f1 > $out_dir/smoove-only.ids
grep -w -f $out_dir/smoove-only.ids $out_dir/smoove-DEL.sv >> $out_dir/DEL.sv
rm $out_dir/smoove-DEL.sv
rm $out_dir/smoove-only.ids


# INSERTIONS

bcftools view $manta_sv | $python manta2sv.py | awk '$8=="INS" && NF==9 || $8=="DUP"' > $out_dir/INS.sv

cat $transurv_sv | $python trans2sv.py | $python trans2ins.py $ref_genome | sed 's/SURVEYOR/TRANSURVEYOR/g' > $out_dir/surv-INS.sv
cat $surv_sv | $python surv2sv.py | grep DUP | awk '$3<$6' >> $out_dir/surv-INS.sv

./compare-ins $out_dir/surv-INS.sv $out_dir/INS.sv $repeats $ref_genome $max_is | grep NONE | cut -d" " -f1 > $out_dir/surv-only.ids
grep -w -f $out_dir/surv-only.ids $out_dir/surv-INS.sv >> $out_dir/INS.sv
rm $out_dir/surv-only.ids
rm $out_dir/surv-INS.sv
