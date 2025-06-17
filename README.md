# Word Embeddings Simulation Study

**Investigating whether word embeddings can capture genuine structural relationships**

[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Overview

This repository contains the complete implementation of a microsimulation study that investigates whether word embeddings derived from text corpora can capture genuine structural relationships or merely represent sophisticated statistical artifacts. 

The study addresses the fundamental **grounding problem** in computational linguistics: *Can distributional patterns in text reliably represent real-world semantic structures?*

## Research Question

**Can word embedding algorithms reconstruct underlying spatial and conceptual structures solely from statistical co-occurrence patterns in text?**

To answer this, we created a controlled experiment using a virtual "supermarket" with a known, deterministic layout and generated 80,000 unique shopping sequences to train word embeddings.

## Key Findings

✅ **Successfully reconstructed spatial layout**: GloVe embeddings perfectly recovered the original 2D supermarket structure  
✅ **Identified semantic clusters**: Products correctly grouped by department (fruits, vegetables, beverages, extras)  
✅ **Preserved proximity relationships**: Products from adjacent shelves showed high similarity scores  
✅ **Robust across variants**: Results replicated across different product types (I, II, III)  

These results provide evidence that word embeddings can capture genuine structural relationships, not just superficial statistical patterns.

## Methodology

### 1. Microsimulation (Virtual Supermarket)
- **40 shelves** organized in 4 departments
- **120 unique products** (3 variants per shelf)
- **80,000 simulated shoppers** with probabilistic movement
- **Stochastic product selection** (15% chance per shelf visit)

### 2. Text Corpus Generation
- Generated **80,000 unique sentences** from shopping sequences
- Total of **804,422 product selections**
- Average **10.1 products per shopper**

### 3. Word Embedding Training
- **GloVe algorithm** with 50-dimensional vectors
- **Co-occurrence matrix** with context window of 5 words
- **Cosine similarity** for measuring relatedness

### 4. Validation & Analysis
- **k-medoid clustering** for department identification
- **UMAP dimensionality reduction** for spatial reconstruction
- **Similarity analysis** for proximity relationships

## Repository Structure

```
├── 00__sim.R           # Microsimulation of supermarket shoppers
├── 01__emb.R           # GloVe word embedding generation
├── 02__process.R       # Embedding analysis and UMAP visualization
├── 03__superheat.R     # Similarity matrix heatmaps and clustering
├── 04__igraph.R        # Network visualization of supermarket layout
├── 05__tcm.R           # Term co-occurrence matrix visualization
├── 06__cluster.R       # k-medoid clustering analysis
├── doc/
│   ├── transition_matrix.xlsx    # Supermarket transition probabilities
│   └── ...
└── README.md
```

## Requirements

### R Packages
```r
# Core analysis
install.packages(c("text2vec", "umap", "cluster"))

# Visualization
install.packages(c("ggplot2", "ggrepel", "superheat", "igraph"))

# Data processing
install.packages(c("readxl", "dplyr", "wordspace", "readr"))

# Additional utilities
install.packages(c("reshape2", "stringr", "scales", "matrixStats"))
```

## Usage

### Quick Start
Run the scripts in numerical order:

```bash
# 1. Generate simulation data
Rscript 00__sim.R

# 2. Train word embeddings
Rscript 01__emb.R

# 3. Analyze embeddings
Rscript 02__process.R

# 4. Create visualizations
Rscript 03__superheat.R
Rscript 04__igraph.R
Rscript 05__tcm.R
Rscript 06__cluster.R
```

### Output Directories
Each script creates its own output directory:
- `out.00.sim/` - Simulation results and corpus
- `out.01.emb/` - Word embeddings and vocabulary
- `out.02.process/` - UMAP projections and similarity analysis
- `out.03.superheat/` - Heatmaps and clustering results
- `out.04.igraph/` - Network visualizations
- `out.05.tcm/` - Co-occurrence matrices
- `out.06.cluster/` - Clustering analysis

## Key Results

### Spatial Reconstruction
The 2D UMAP projection of the 50-dimensional embeddings perfectly reconstructed the original supermarket layout, demonstrating that statistical co-occurrence patterns preserved genuine spatial relationships.

### Semantic Clustering
k-medoid clustering identified 4 clusters corresponding to the supermarket departments:
- **Fruits** (11 products) → Lime cluster
- **Vegetables** (8 products) → Onion cluster  
- **Beverages** (7 products) → Lemonade cluster
- **Extras** (14 products) → GummyBears cluster

### Similarity Patterns
Products from the same shelf showed similarity scores > 0.9, while distant products scored < 0.3, reflecting the spatial structure of the supermarket.

## Theoretical Implications

This study provides empirical evidence for the **distributional hypothesis** and suggests that:

1. **Information compression** via embeddings can preserve structural relationships
2. **Optimal data description** aligns with genuine semantic representation
3. **Distributional patterns** can ground computational semantics in real-world structures

## Citation

If you use this code or methodology in your research, please cite:

```bibtex
@article{otte2025embeddings,
  title={Word embeddings from text corpora: a simulation study on the representation of underlying structures},
  author={Otte, Willem M. and van Wieringen, Archibald L.H.M. and Koet, Bart J.},
  journal={[Journal Name]},
  year={2025},
  note={Manuscript in preparation}
}
```

## Contributing

We welcome contributions! Please feel free to:
- Report bugs or issues
- Suggest improvements to the methodology
- Extend the simulation to other domains
- Add alternative embedding algorithms

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details.

## Authors

- **Willem M. Otte** - *Principal Investigator* - Utrecht University & UMC Utrecht
- **Archibald L.H.M. van Wieringen** - Tilburg University
- **Bart J. Koet** - Tilburg University

## Acknowledgments

This research was inspired by Shannon's information theory and the distributional hypothesis of Harris and Firth. The microsimulation methodology draws from decision models in health economics.

---

*For questions or collaborations, please contact: w.m.otte@umcutrecht.nl*
