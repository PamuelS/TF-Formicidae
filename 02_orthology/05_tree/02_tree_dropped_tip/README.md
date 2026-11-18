# Pruned tree
Sono state droppate 5 specie dall'albero precedentemente originato. Questo per evitare bias legati alla ripetizione di specie.
Il seguente codice è stato utilizzato per rimuovere le specie in eccesso.

```bash
library(ape)


input_tree_path <- "tree_rooted.nwk"
output_tree_path <- "tree_pruned.nwk"

tree <- read.tree(input_tree_path)

cat("numero iniziale di specie:", Ntip(tree), "\n")
head(tree$tip.label)

species_to_drop <- c("Lasfla", "Myrcro", "Caresp", "Lasful", "Psearg")

not_found <- setdiff(species_to_drop, tree$tip.label)
if (length(not_found) > 0) {
  "attenzione queste etichette non sono corrette"
}

tree_dropped <- drop.tip(tree, species_to_drop)

cat("numero finale di specie:", Ntip(tree_dropped), "\n")

write.tree(tree_dropped, file = output_tree_path)
```
