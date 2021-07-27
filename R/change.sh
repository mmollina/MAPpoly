find ./ -type f -exec sed -i 's/$ploidy /$ploidy /g' {} \;
find ./ -type f -exec sed -i 's/dosage.p1/dosage.p11/g' {} \;
find ./ -type f -exec sed -i 's/dosage.p2/dosage.p12/g' {} \;
find ./ -type f -exec sed -i 's/$genome.pos/$genome.pos/g' {} \;
find ./ -type f -exec sed -i 's/$chrom/$chrom/g' {} \;
find ./ -type f -exec sed -i 's/$genome.pos/$genome.pos/g' {} \;
