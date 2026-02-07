g_l=($(ls ..//simulated_treeseq/ | grep 500Kb | cut -f1 -d'.'| grep admixture | grep rate))

mkdir -p ../inferred_genetrees
mkdir -p /dev/shm/admixture/inferred_genetrees

for g in "${g_l[@]}"
do
    echo $g
    if [[ -f ../inferred_genetrees/$g.gtrees ]]; then
        continue
    else
      for i in {1..125}
      do  
          echo "python ../extract_genes.py \
              -i /dev/shm/admixture/seq/$g.fasta \
              --group-size 8 --group-index $i --gene-index 2 --gene-length 500 \
              > /dev/shm/admixture/inferred_genetrees/$g-$i.fasta && iqtree2 -s /dev/shm/admixture/inferred_genetrees/$g-$i.fasta -m GTR+G4 -abayes -nt 1 >/dev/null 2>&1" >> commands.txt
          echo "Inferring gene trees for ${i}..."
      done
      wait

      cat /dev/shm/admixture/inferred_genetrees/${g}-*.treefile > /dev/shm/inferred_genetrees/$g.gtrees
      rm /dev/shm/admixture/inferred_genetrees/$g-*.fasta
      rm /dev/shm/admixture/inferred_genetrees/$g-*.bionj
      rm /dev/shm/admixture/inferred_genetrees/$g-*.mldist
      rm /dev/shm/admixture/inferred_genetrees/$g-*.treefile
      rm /dev/shm/admixture/inferred_genetrees/$g-*.fasta.ckp.gz 
      rm /dev/shm/admixture/inferred_genetrees/$g-*.uniqueseq.phy

      mkdir -p ../logs/logs_iqtree-${g}/

      mv /dev/shm/admixture/inferred_genetrees/$g-*.iqtree ../logs/logs_iqtree-${g}/
      mv /dev/shm/admixture/inferred_genetrees/$g-*.log ../logs/logs_iqtree-${g}/
      tar -czvf ../logs/logs_iqtree-$g.tar.gz ../logs/logs_iqtree-${g}/ >/dev/null 2>&1
      rm -r ../logs/logs_iqtree-${g}/
    fi
done
