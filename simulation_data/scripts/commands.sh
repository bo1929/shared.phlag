# Simulated emission distribution 10X population size
seq 1 50 | xargs -t -P 12 -I{} bash -c "bash compute_nullDist.sh \
  ../species_trees/estimated/SR201_10X_population/{}_1X_S201_0_1e6_haploid.caster-pair \
  ../gene_trees/SR201_10X_population/{}_1X_S201_0_1e6_haploid.gtrees \
  ../qqs-SR201/10X_population-{}"

# Emissions 10X population size
seq 1 50 | xargs -t -P 12 -I{} bash -c "bash compute_freqQuad.sh \
  ../qqs-SR201/10X_population-{}/labelled_cu_tree.tree \
  ../gene_trees/SR201_10X_population/{}_1X_S201_0_1e6_haploid.gtrees \
  ../qqs-SR201/10X_population-{}"

# Emissions default population size
seq 1 50 | xargs -t -P 12 -I{} bash -c "bash compute_nullDist.sh \
  ../species_trees/estimated/SR201_default_condition/{}_1X_S201_0_haploid.caster-pair \
  ../gene_trees/SR201_default_condition/{}_1X_S201_0_haploid.gtrees \
  ../qqs-SR201/default_condition-{}"

# Simulated emission distribution default population size
seq 1 50 | xargs -t -P 12 -I{} bash -c "bash compute_freqQuad.sh \
  ../qqs-SR201/default_condition-{}/labelled_cu_tree.tree \
  ../gene_trees/SR201_default_condition/{}_1X_S201_0_haploid.gtrees \
  ../qqs-SR201/default_condition-{}"
