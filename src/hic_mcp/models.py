"""Typed response models - the stable output contract every tool returns.

server.py wraps the plain dicts from analysis.py into these, so the MCP
outputSchema clients see is explicit and versioned with the package.
"""

from pydantic import BaseModel, Field


class ResolutionInfo(BaseModel):
    resolution_bp: int
    bins: int
    nonzero_pixels: int
    total_contacts: int
    balanced: bool


class MatrixSummary(BaseModel):
    file: str
    is_bundled_demo: bool
    assembly: str
    chromosomes: dict[str, int] = Field(description="chromosome name -> length in bp")
    resolutions: list[ResolutionInfo]
    balanced: bool
    provenance: dict = Field(default_factory=dict)
    method: str


class ContactsAtLocus(BaseModel):
    region: str
    region2: str | None = None
    resolution_used: int
    shape_bins: list[int]
    raw_contacts_sum: int
    raw_contacts_max: int
    nonzero_fraction: float
    balanced: bool
    balanced_mean: float | None = None
    balanced_max: float | None = None
    balanced_matrix: list[list[float | None]] | None = None
    note: str | None = None
    method: str


class Boundary(BaseModel):
    locus: str
    strength: float | None
    log2_insulation: float | None
    windows_detected: list[int]


class InsulationTads(BaseModel):
    region: str | None = None
    resolution_used: int
    view: str = Field(description="The region model the score was normalised within")
    windows_bp: list[int]
    ranked_by: str
    boundary_counts_per_window: dict[str, int]
    top_boundaries: list[Boundary]
    balanced: bool
    scale_note: str | None = Field(
        default=None, description="Set when the windows are too large for TAD-scale calls"
    )
    method: str


class ArmEigenvalue(BaseModel):
    region: str
    eigval1: float | None
    variance_share: float | None = Field(
        default=None, description="Scale-free counterpart to eigval1; compare regions by this"
    )


class E1Point(BaseModel):
    start: int
    E1: float | None


class Compartments(BaseModel):
    resolution_used: int
    view: str
    sign_convention: str
    eigenvalues: list[ArmEigenvalue]
    eigenvalue_note: str | None = None
    balanced: bool
    region: str | None = None
    region_mean_E1: float | None = None
    region_call: str | None = None
    region_sign_consistency: float | None = Field(
        default=None, description="Null when too few bins for the figure to mean anything"
    )
    confidence_note: str | None = None
    transition_note: str | None = None
    E1_track: list[E1Point] | None = None
    genome_A_fraction: float | None = Field(
        default=None, description="Only set when the eigenvector is phased; null otherwise"
    )
    positive_E1_fraction: float | None = Field(
        default=None, description="Sign-neutral counterpart reported when unphased"
    )
    bins_used: int | None = None
    method: str


class ProfilePoint(BaseModel):
    start: int
    balanced: float | None


class Virtual4C(BaseModel):
    viewpoint: str
    resolution_used: int
    profile_points: list[ProfilePoint]
    profile_note: str
    distance_band_means: dict[str, float | None]
    distance_band_bins: dict[str, int] = Field(default_factory=dict)
    distance_bands_cover_bp: list[int] | None = None
    coverage_note: str | None = None
    balanced: bool
    method: str


class ExpectedCurvePoint(BaseModel):
    dist_bp: int
    expected: float | None


class ExpectedObserved(BaseModel):
    region: str
    resolution_used: int
    view: str = Field(description="The region whose expected curve this is")
    curve_scope: str = Field(description="States what the curve and slope describe")
    ps_slope: float | None = Field(
        default=None, description="log-log slope of P(s) over the range actually fitted"
    )
    ps_fit_range_bp: list[int] | None = Field(
        default=None, description="[min, max] separation the slope was fitted over"
    )
    ignored_diagonals: int = Field(
        description="Diagonals with no expected measurement; O/E is null there"
    )
    expected_curve_points: list[ExpectedCurvePoint]
    balanced: bool
    oe_matrix: list[list[float | None]] | None = None
    note: str | None = None
    method: str
