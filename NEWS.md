## Changes in version 0.1.0 (2024-11-08)

- **Fixed**: Resolved a NAMESPACE issue that was affecting documentation visibility for certain functions.

## Changes in version 0.1.0 (2024-11-12)

- **PruneGatingTree Function**: Corrected an issue where the pruning logic did not properly adhere to initial pruning criteria and failed to consistently retain parent nodes when any of their child nodes met retention conditions. This fix ensures that the function now correctly applies specified thresholds and respects the hierarchical structure of the gating tree.


## Changes in Version 0.1.1 (2024-11-17)

- **PlotDeltaEnrichment** and **PlotDeltaEnrichmentPrunedTree**: The option `significant` has been removed. Instead, a new function `PlotDeltaEnrichmentPrunedTree` is now available to assess the impacts of important markers.