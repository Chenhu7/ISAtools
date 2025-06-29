#!/bin/bash

get_seeded_random() {
    seed="$1"
    openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt </dev/zero 2>/dev/null
}

SEED=123

for file in *.gtf; do
    input_gtf=$file
    bn=$(basename "$file" .gtf)

    echo "Processing $input_gtf"

    grep -w transcript "$input_gtf" | \
    awk -F'\t' '{
        split($9,a,";");
        for(i in a) {
            if(a[i] ~ /transcript_id/) {tid=a[i]; gsub(/.*transcript_id "/,"",tid); gsub(/"$/,"",tid)}
            if(a[i] ~ /gene_id/) {gid=a[i]; gsub(/.*gene_id "/,"",gid); gsub(/"$/,"",gid)}
        }
        print tid, gid
    }' OFS='\t' > transcript_gene_pairs.txt

    awk '{count[$2]++; lines[$2] = lines[$2] $0 "\n"} END {
        for (g in count) {
            if (count[g] > 1) {
                split(lines[g], l, "\n");
                for (i in l) if (l[i] != "") print l[i];
            }
        }
    }' transcript_gene_pairs.txt > multi_transcript_pairs.txt

    total=$(wc -l < transcript_gene_pairs.txt)
    total_multi=$(wc -l < multi_transcript_pairs.txt)
    delete_count20=$((total * 20 / 100))
    delete_count40=$((total * 40 / 100))

    echo "Total transcripts: $total, from multi-transcript genes: $total_multi"
    echo "Will delete: 20% = $delete_count20; 40% = $delete_count40"

    for percent in 20 40; do
        output="remove_list${percent}.txt"
        delete_count=$((total * percent / 100))

        awk -v seed="$SEED" -v target="$delete_count" '
        BEGIN {
            srand(seed)
        }
        {
            tids[$2][++count[$2]] = $1
        }
        END {
            for (g in count) {
                genes[++gn] = g
            }

            if (target >= gn) {
                for (i = 1; i <= gn; i++) {
                    g = genes[i]
                    n = int(rand() * count[g]) + 1
                    removed[tids[g][n]] = 1
                    print tids[g][n]
                    deleted++
                }
                remain = target - deleted
                if (remain > 0) {
                    for (i = 1; i <= gn && remain > 0; i++) {
                        g = genes[i]
                        for (j = 1; j <= count[g] && remain > 0; j++) {
                            tid = tids[g][j]
                            if (!(tid in removed)) {
                                print tid
                                removed[tid] = 1
                                remain--
                            }
                        }
                    }
                }
            } else {
                for (i = gn; i >= 1; i--) {
                    j = int(rand() * i) + 1
                    tmp = genes[i]; genes[i] = genes[j]; genes[j] = tmp
                }

                for (i = 1; i <= target; i++) {
                    g = genes[i]
                    n = int(rand() * count[g]) + 1
                    print tids[g][n]
                }
            }
        }' multi_transcript_pairs.txt > "$output"
    done

    grep -F -v -f remove_list20.txt "$input_gtf" > "$bn.kept20.gtf"
    grep -F -v -f remove_list40.txt "$input_gtf" > "$bn.kept40.gtf"

    rm multi_transcript_pairs.txt remove_list*.txt transcript_gene_pairs.txt
 
    echo "Finished processing $input_gtf"
done
