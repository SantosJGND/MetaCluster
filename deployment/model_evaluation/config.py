"""
Configuration classes using Pydantic v2.

Provides validated configuration for the evaluator module with
support for CLI args, YAML, JSON, and environment variables.
"""

from pathlib import Path
from typing import Optional, Any, Literal
from dataclasses import field

from pydantic import BaseModel, Field, field_validator, model_validator, ConfigDict
import argparse


ALLOWED_TAX_LEVELS = {'order', 'family', 'genus', 'species'}


class EvaluatorConfig(BaseModel):
    """
    Main configuration for evaluator.
    
    Uses Pydantic v2 for automatic validation and type coercion.
    """
    
    model_config = ConfigDict(
        extra='forbid',
        str_strip_whitespace=True,
        validate_assignment=True,
    )
    
    study_output_filepath: Path
    taxid_plan_filepath: Path
    analysis_output_filepath: Path
    
    data_set_divide: int = Field(default=5, ge=1, le=20)
    tax_level: str = Field(default='order')
    cross_hit_threshold: float = Field(default=0.9, ge=0.0, le=1.0)
    taxa_threshold: float = Field(default=0.02, ge=0.0, le=1.0)
    holdout_proportion: float = Field(default=0.3, ge=0.0, le=1.0)
    max_training: Optional[int] = Field(default=None, ge=1)
    threshold: float = Field(default=0.3, ge=0.0, le=1.0)
    output_db_dir: Optional[Path] = None
    
    verbose: bool = Field(default=False, description="Enable verbose console logging")
    generate_report: bool = Field(default=True, description="Generate HTML report")
    use_mlflow: bool = Field(default=True, description="Enable MLflow tracking")
    mlflow_uri: Optional[str] = Field(default=None, description="MLflow tracking URI")
    
    @field_validator('tax_level')
    @classmethod
    def validate_tax_level(cls, v: str) -> str:
        if v.lower() not in ALLOWED_TAX_LEVELS:
            raise ValueError(f'tax_level must be one of {ALLOWED_TAX_LEVELS}, got: {v}')
        return v.lower()
    
    @field_validator('cross_hit_threshold', 'taxa_threshold', 'holdout_proportion', 'threshold')
    @classmethod
    def validate_probability(cls, v: float) -> float:
        if not 0.0 <= v <= 1.0:
            raise ValueError(f'Value must be between 0.0 and 1.0, got: {v}')
        return v
    
    @model_validator(mode='after')
    def validate_paths(self):
        if not self.study_output_filepath.exists():
            raise ValueError(f'study_output_filepath does not exist: {self.study_output_filepath}')
        if not self.taxid_plan_filepath.exists():
            raise ValueError(f'taxid_plan_filepath does not exist: {self.taxid_plan_filepath}')
        self.analysis_output_filepath.mkdir(parents=True, exist_ok=True)
        return self
    
    @property
    def proportion_train(self) -> float:
        return 1.0 - self.holdout_proportion
    
    @property
    def models_dir(self) -> Path:
        return self.analysis_output_filepath / "models"
    
    @property
    def output_lineages(self) -> Path:
        return self.study_output_filepath / "lineages.tsv"
    
    @property
    def output_db(self) -> Path:
        if self.output_db_dir:
            return self.output_db_dir / "taxa.db"
        return self.study_output_filepath / "taxa.db"
    
    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'EvaluatorConfig':
        """Create config from argparse.Namespace."""
        return cls(
            study_output_filepath=Path(args.study_output_filepath),
            taxid_plan_filepath=Path(args.taxid_plan_filepath),
            analysis_output_filepath=Path(args.analysis_output_filepath),
            data_set_divide=args.data_set_divide,
            tax_level=args.tax_level_to_use,
            cross_hit_threshold=args.cross_hit_threshold,
            taxa_threshold=args.taxa_threshold,
            holdout_proportion=args.holdout_proportion,
            max_training=int(args.max_training) if args.max_training else None,
            threshold=args.threshold,
            output_db_dir=Path(args.output_db_dir) if args.output_db_dir else None,
        )
    
    @classmethod
    def from_yaml(cls, path: str | Path) -> 'EvaluatorConfig':
        """Load config from YAML file."""
        import yaml
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)
    
    @classmethod
    def from_json(cls, path: str | Path) -> 'EvaluatorConfig':
        """Load config from JSON file."""
        import json
        with open(path) as f:
            data = json.load(f)
        return cls(**data)
    
    def to_yaml(self, path: str | Path) -> None:
        """Save config to YAML file."""
        import yaml
        with open(path, 'w') as f:
            yaml.dump(self.model_dump(), f, default_flow_style=False)
    
    def to_json(self, path: str | Path) -> None:
        """Save config to JSON file."""
        import json
        with open(path, 'w') as f:
            json.dump(self.model_dump(), f, indent=2)


class ModelConfig(BaseModel):
    """
    Configuration for model training hyperparameters.
    """
    
    model_config = ConfigDict(extra='allow')
    
    recall_hyperparams: dict = field(default_factory=lambda: {
        'n_estimators': 100,
        'max_depth': 5,
        'learning_rate': 0.1,
        'objective': 'binary:logistic',
    })
    
    composition_hyperparams: dict = field(default_factory=lambda: {
        'n_estimators': 50,
        'max_depth': 3,
        'learning_rate': 0.1,
    })
    
    crosshit_hyperparams: dict = field(default_factory=lambda: {
        'probability_threshold': 0.9,
        'n_estimators': 50,
        'max_depth': 3,
    })
    
    k_folds: int = Field(default=5, ge=2, le=10)
    test_size: float = Field(default=0.2, ge=0.1, le=0.5)
    random_state: int = Field(default=42)


class VisualizationConfig(BaseModel):
    """
    Configuration for visualization output.
    """
    
    model_config = ConfigDict(extra='allow')
    
    style: str = Field(default='seaborn-v0_8-darkgrid')
    figure_dpi: int = Field(default=300, ge=72, le=600)
    figure_format: Literal['png', 'svg', 'pdf', 'jpg'] = Field(default='png')
    interactive: bool = Field(default=False)
    html_template: Optional[str] = None
    width: int = Field(default=12, ge=6, le=24)
    height: int = Field(default=8, ge=4, le=16)
    color_palette: str = Field(default='viridis')
    font_scale: float = Field(default=1.0, ge=0.5, le=2.0)


class MLflowConfig(BaseModel):
    """
    Configuration for MLflow tracking.
    """
    
    model_config = ConfigDict(extra='allow')
    
    enabled: bool = Field(default=False)
    tracking_uri: Optional[str] = None
    experiment_name: str = Field(default='metagenomics-evaluation')
    run_name: Optional[str] = None
    log_models: bool = Field(default=True)
    log_datasets: bool = Field(default=True)
    artifact_location: Optional[str] = None


class LoggingConfig(BaseModel):
    """
    Configuration for logging.
    """
    
    model_config = ConfigDict(extra='allow')
    
    level: Literal['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'] = Field(default='INFO')
    format: Literal['json', 'text'] = Field(default='text')
    log_file: Optional[Path] = None
    console_level: str = Field(default='INFO')
    file_level: str = Field(default='DEBUG')


class TrainedModels:
    """
    Container for trained model instances.
    """
    
    def __init__(
        self,
        recall_modeller: Any = None,
        composition_modeller: Any = None,
        crosshit_modeller: Any = None,
    ):
        self.recall_modeller = recall_modeller
        self.composition_modeller = composition_modeller
        self.crosshit_modeller = crosshit_modeller
    
    def __repr__(self) -> str:
        return (
            f"TrainedModels("
            f"recall={self.recall_modeller is not None}, "
            f"composition={self.composition_modeller is not None}, "
            f"crosshit={self.crosshit_modeller is not None})"
        )
    
    @property
    def is_complete(self) -> bool:
        return all([
            self.recall_modeller is not None,
            self.composition_modeller is not None,
            self.crosshit_modeller is not None,
        ])
