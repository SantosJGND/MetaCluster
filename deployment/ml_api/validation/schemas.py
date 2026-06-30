from pydantic import BaseModel, Field
from typing import Dict, List, Optional


class RecallCutoffFromTableRow(BaseModel):
    taxid: int
    total_uniq_reads: float
    best_match_is_best: bool = False


class RecallCutoffFromTableRequest(BaseModel):
    model: str = "gp_clf"
    rows: List[RecallCutoffFromTableRow]
    tax_level: str = "order"
    target_recall: Optional[float] = None
    confidence: Optional[float] = None


class RecallCutoffResult(BaseModel):
    predicted_cutoff: float = Field(..., description="Predicted recall cutoff fraction (target_percentile)")
    predicted_recall_at_cutoff: float = Field(
        ..., description="Predicted recall value at the predicted cutoff fraction"
    )
    target_recall: Optional[float] = Field(
        None, description="Target recall τ that was applied"
    )
    confidence_score: float = Field(..., description="Confidence score of the prediction")
    model_version: Optional[str] = Field(None, description="Version of the model used")
    model_stage: Optional[str] = Field(None, description="MLflow stage of the model used")


CLUSTERING_TAXON_ORDER = [
    "Baculoviridae", "Coronaviridae", "Flaviviridae", "Herelleviridae",
    "Orthoherpesviridae", "Papillomaviridae", "Peduoviridae", "Picornaviridae",
    "Rhabdoviridae", "Schitoviridae", "Straboviridae", "unclassified",
]


class TelevirClusteringThresholdRequest(BaseModel):
    classifiers: List[str] = Field(..., description="List of classifiers to use for prediction")
    taxa: List[str] = Field(..., description="List of taxa to consider for prediction")
    taxon_hits: List[float] = Field(
        default_factory=lambda: [0.0] * len(CLUSTERING_TAXON_ORDER),
        description="Taxonomic proportion values in CLUSTERING_TAXON_ORDER order (12 values)",
    )
    taxonomic_diversity: float = Field(..., description="Taxonomic diversity of the sample")
    n_leaves: float = Field(..., description="Number of leaves in the cluster")
    proportion_shared_hits: float = Field(..., description="Proportion of hits shared with other clusters")
    proportion_unique_hits: float = Field(..., description="Proportion of hits unique to the cluster")

    def _feature_map(self):
        d = {}
        d["taxonomic_diversity"] = self.taxonomic_diversity
        d["n_leaves"] = self.n_leaves
        d["proportion_shared_hits"] = self.proportion_shared_hits
        d["proportion_unique_hits"] = self.proportion_unique_hits
        for name, val in zip(CLUSTERING_TAXON_ORDER, self.taxon_hits):
            d[name] = val
        return d

    def to_array(self, feature_names=None):
        if feature_names is None:
            return list(self._feature_map().values())
        fm = self._feature_map()
        return [fm[name] for name in feature_names]


class ClusteringThresholdResult(BaseModel):
    is_cluster: bool = Field(..., description="Indicates whether the cluster is valid or not")
    confidence_score: float = Field(..., description="Confidence score of the prediction")
    model_version: Optional[str] = Field(None, description="Version of the model used")
    model_stage: Optional[str] = Field(None, description="MLflow stage of the model used")


class CompositionStopTraversalRequest(BaseModel):
    features: Dict[str, float] = Field(
        ...,
        description="Feature dict — keys matching training column names, values are feature values at the node",
    )


class CompositionStopTraversalResult(BaseModel):
    stop_traversal: bool = Field(..., description="Whether traversal should stop at this node")
    probability: float = Field(..., description="Probability of the stop_traversal class")
    confidence_score: float = Field(..., description="Confidence score of the prediction")
    model_version: Optional[str] = Field(None, description="Version of the model used")
    model_stage: Optional[str] = Field(None, description="MLflow stage of the model used")