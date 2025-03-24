# GatingTree NEWS

## Changes in Version 0.2.0 (2025-03-23)
- **New Features**:
  - **GatingTreeRandomForest**: New function to apply Random Forest analysis to GatingTree node data.
  - **predictGatingTreeRandomForest**: New function to perform predictions using models built with GatingTreeRandomForest.
  - **AnalyzeNodeOverlaps**: New function to analyze the percentage of cells shared between nodes, aiding in the unsupervised analysis of GatingTree node data.
  - **ExtractTopNodes**: New function designed to extract top nodes from node clusters, enhancing the unsupervised approach in node analysis.

## Changes in Version 0.1.1 (2024-11-17)
- **Enhancements**:
  - Removed the `significant` option from **PlotDeltaEnrichment**.
  - Introduced **PlotDeltaEnrichmentPrunedTree**, a new function designed to provide a more powerful assessment of the impacts of important markers in pruned gating trees.

## Changes in Version 0.1.0 (2024-11-12)
- **Bug Fixes**:
  - **PruneGatingTree Function**: Fixed an issue where the pruning logic did not properly adhere to the initial pruning criteria. The function now correctly retains parent nodes when any of their child nodes meet retention conditions, respecting the hierarchical structure of the gating tree.

## Changes in Version 0.1.0 (2024-11-08)
- **Fixed**:
  - Resolved a NAMESPACE issue that was affecting documentation visibility for certain functions.
