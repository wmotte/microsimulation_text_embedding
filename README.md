# Word Embeddings Simulation Study

**How much of a known generating structure does a word embedding retain?**

[![R](https://img.shields.io/badge/R-4.6-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE.md)

This repository contains the code, the generated corpus and all results for:

> Otte WM, van Wieringen ALHM, Koet BJ. *Word embeddings from text corpora: a simulation study on the representation of underlying structures.* Submitted to PLOS ONE (PONE-D-26-11851), revised version 2026.

## Overview

For natural text, the structure that produced a corpus is unknown, so we cannot check how much of it a word embedding retains. Here the structure is known exactly. Simulated shoppers walk through a virtual supermarket of 40 shelves (120 products in four departments, connected in a fixed network), and the products they pick up form the corpus. We train GloVe embeddings on that corpus and compare them with the network.

The corpus is not natural language. Its "sentences" are lists of purchased products, without syntax or function words. The study is a positive control: it measures how much of a simple, known structure survives in the co-occurrence statistics and in the embedding, and what limits that. It does not show that embeddings understand or are grounded in what words refer to.

## Main results

All numbers are for the main corpus (80,000 sentences, 804,422 tokens) unless stated otherwise.

| Measure | Result |
|---|---|
| Spearman ρ, cosine distance vs shortest-path distance (7,140 pairs) | 0.79 (Mantel p = 0.001) |
| Neighbour retrieval (chance 0.067) | 0.61 |
| Ten independent simulations (corpus + model) | ρ = 0.80 ± 0.01, retrieval 0.62 ± 0.02 |
| Baselines: raw counts / PPMI / SVD(PPMI) / word2vec | ρ 0.90 / 0.81 / 0.66 / 0.80; retrieval 0.40 / 0.50 / 0.29 / 0.41 |
| Controls: shuffled across corpus / within sentences | ρ −0.04 / 0.14; retrieval 0.09 / 0.40 |
| Context window 1 / 5 / 10 | ρ 0.93 / 0.80 / 0.71 |
| UMAP map of type I products | ρ = 0.92, trustworthiness (k = 5) 0.96 |
| PAM clustering of all 120 products vs departments | ARI 0.82 (fruit and extras pure; vegetables and beverages mixed) |
| Product frequency | R² of cosine on distance 0.60, with log frequencies 0.80 |

## Running the analysis

Requires R ≥ 4.6. Run the scripts from the repository root, in this order:

```bash
Rscript 00__sim.R          # simulation and corpus
Rscript 01__emb.R          # GloVe embedding
Rscript 03__superheat.R    # embedding and similarity heatmaps
Rscript 05__tcm.R          # co-occurrence matrix and Tables 3-4
Rscript 06__cluster.R      # PAM clustering
Rscript 02__process.R      # nearest neighbours (Table 6), UMAP maps
Rscript 07__recovery.R     # quantitative recovery, frequency analysis, Fig 8
Rscript 08__robustness.R   # 10 seeds, baselines, controls, sweeps (slow: several hours)
Rscript 09__fig_robustness.R   # summaries and Fig 10
```

`04__igraph.R` draws the network (Fig 1) and does not depend on the other scripts. `functions.R` holds the shared functions: a vectorised simulator for the repetitions, co-occurrence counting, GloVe, PPMI, SVD, word2vec and the recovery measures.

### Packages

| Package | Version | Used for |
|---|---|---|
| text2vec | 0.6.6 | GloVe, co-occurrence matrix |
| rsparse | 0.5.3 | GloVe backend |
| word2vec | 0.4.1 | skip-gram baseline |
| umap | 0.2.10.0 | two-dimensional maps |
| cluster | 2.1.8.2 | PAM, silhouette |
| igraph | 2.3.3 | shortest-path distances |
| mclust | | adjusted Rand index |
| vegan | 2.7-6 | Mantel test, Procrustes |
| irlba | | SVD baseline |
| superheat | GitHub `rlbarter/superheat` | heatmaps |
| ggplot2, ggrepel, readxl, wordspace | | figures, input, similarity |

The exact versions of every run are in `out.01.emb/session_info.txt`, `out.07.recovery/session_info.txt` and `out.08.robustness/session_info.txt`.

### Random seeds and reproducibility

| Step | Seed |
|---|---|
| Walk of shopper *i* (`00__sim.R`) | `set.seed(1 + i)` |
| Product choices (`00__sim.R`) | `set.seed(4321)` |
| GloVe (`01__emb.R`) | 444, on one thread (`n_threads = 1`) |
| UMAP, main maps (`02__process.R`) | `random_state = 123` |
| Repetitions *s* = 1–10 (`08__robustness.R`) | corpus 1000 + *s*, GloVe *s*, shuffles 2000 + *s*, UMAP *s* |

With these seeds, the scripts reproduce the results exactly. GloVe must run on one thread, because multi-threaded updates are not reproducible from a seed. The one exception is word2vec: the R package has no seed argument, so repeated fits differ slightly (by less than 0.01 in ρ and retrieval on the main corpus).

The repetitions use `simulate_corpus_fast()` in `functions.R`, a vectorised version of the same model. It reproduces the main simulation: 35.2 ± 0.1% of shoppers still inside at step 75 (main run 35.3%), and 805,445 ± 549 tokens (main run 804,422).

## Where each table and figure comes from

| Paper | Script | Output |
|---|---|---|
| Fig 1 | `04__igraph.R` | `out.04.igraph/Figuur_1.png` |
| Fig 2 | `00__sim.R` | `out.00.sim/network_flow.png` |
| Fig 3 | `00__sim.R` | `out.00.sim/states_before_exit.png` |
| Fig 4 | `00__sim.R` | `out.00.sim/product_list_length.png` |
| Fig 5 | `05__tcm.R` | `out.05.tcm/tcm_I.png` |
| Fig 6 | `03__superheat.R` | `out.03.superheat/wv_embedding_type_I__items.png`, `..._sections.png` |
| Fig 7 | `03__superheat.R` | `out.03.superheat/wv_embedding_type_I__similarity_matrix.png` |
| Fig 8 | `07__recovery.R` | `out.07.recovery/fig8_truth_vs_umap.png` |
| Fig 9 | `02__process.R` | `out.02.process/map_embedding_type_II.png`, `..._III.png` |
| Fig 10 | `09__fig_robustness.R` | `out.08.robustness/fig10_robustness.png` |
| Table 1, S1 Table | `00__sim.R` | `out.00.sim/transition_matrix.tsv` (input: `doc/transition_matrix.xlsx`) |
| Table 2 | `01__emb.R` | `out.01.emb/vocab_summary.tsv` |
| Tables 3–4 | `05__tcm.R` | `out.05.tcm/table_3_and_4_cooccurrence.tsv` |
| Table 5, S3 Table | `06__cluster.R` | `out.06.cluster/clustering_*.tsv`, `silhouette_by_k_*.tsv`, `graph_partition_vs_departments.tsv` |
| Table 6 | `02__process.R` | `doc/Tabel_6.tsv` |
| Table 7, S5 Table | `08__robustness.R`, `09__fig_robustness.R` | `out.08.robustness/robustness_summary.tsv`, `umap_summary.tsv` |
| S2 Table | `00__sim.R`, `01__emb.R`, `08__robustness.R` | `out.00.sim/filter_counts.tsv`, `out.01.emb/glove_settings.tsv`, `out.08.robustness/corpus_info_summary.tsv` |
| S4 Table | `07__recovery.R`, `06__cluster.R` | `out.07.recovery/neighbour_precision_per_item.tsv`, `out.06.cluster/clustering_details_*.tsv` |
| S1 Fig | `07__recovery.R` | `out.07.recovery/cosine_by_graph_distance.png` |
| S2 Fig | `07__recovery.R` | `out.07.recovery/product_frequency_by_shelf.png` |
| Recovery numbers in the text | `07__recovery.R` | `out.07.recovery/recovery_metrics.tsv`, `umap_metrics_type_*.tsv`, `frequency_effects.tsv`, `department_structure.txt` |
| Corpus counts in the text | `00__sim.R` | `out.00.sim/filter_counts.tsv`, `example_sentences.txt` |

The corpus itself is `out.00.sim/plain_text.txt.gz` (one sentence per line), and the fitted embedding is `out.01.emb/saved_glove.RData`.

## Citation

```bibtex
@article{otte2026embeddings,
  title   = {Word embeddings from text corpora: a simulation study on the representation of underlying structures},
  author  = {Otte, Willem M. and van Wieringen, Archibald L. H. M. and Koet, Bart J.},
  journal = {PLOS ONE},
  year    = {2026},
  note    = {Under revision (PONE-D-26-11851)}
}
```

## License

This project is licensed under the MIT License. See [LICENSE.md](LICENSE.md).

## Authors

- **Willem M. Otte** - *Principal Investigator* - Utrecht University & UMC Utrecht
- **Archibald L.H.M. van Wieringen** - Tilburg University
- **Bart J. Koet** - Tilburg University

---

*For questions or collaborations, please contact: w.m.otte@umcutrecht.nl*
