# Data provenance notes

`lr_pairs_human` is a compact LogicComm-format human ligand-receptor seed set.
It is stored in `data/lr_pairs_human.rda` and documented in `R/lr_database.R` and
`man/lr_pairs_human.Rd`.

The table is intended as a small example/default database. For production use,
users can convert a full CellChatDB-like object with
`as_logiccomm_lr_db_from_cellchat()` and pass the resulting table as `lr_db`.

Required LogicComm schema:

- `lr_pair`: unique interaction identifier
- `ligand_genes`: list column of ligand gene/subunit vectors
- `receptor_genes`: list column of receptor gene/subunit vectors

Recommended metadata columns are `ligand`, `receptor`, `pathway`, and
`annotation`. Optional CellChatDB-derived modulator list columns are
`agonist_genes`, `antagonist_genes`, `co_A_receptor_genes`, and
`co_I_receptor_genes`.
