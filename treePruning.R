#!/usr/local/pacerepov1/R/3.3.2/bin/Rscript

### Usage
## This script remove missing species from a tree and unroot this tree for PAML
# treePruning tree.newick missingSpps.txt
# The "tree.newick" file needs to be newick format
# The "missingSpps.txt" file needs to be in a format like spp1\tspp2\tspp3\n

zz <- file("all.Rout", open="wt")
sink(file = zz, type="message")

### Load 'ape' package
require(ape)
# Read input file with missing species & tree
args = commandArgs(trailingOnly=TRUE)
tree.before <- read.tree(text=args[1])
missing <- scan(args[2], what=list(""))[[1]]
# Prune the tree by trimming missing branches
tree.after = unroot(drop.tip(tree.before, tip=missing))
# Output the new tree
sink()
unlink("all.Rout")
write.table(write.tree(tree.after), col.names = F, row.names = F, quote = F)