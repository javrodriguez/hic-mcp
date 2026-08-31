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
    windows_bp: list[int]
    ranked_by: str
    boundary_counts_per_window: dict[str, int]
    top_boundaries: list[Boundary]
    balanced: bool
    method: str


class ArmEigenvalue(BaseModel):
    region: str
    eigval1: float | None


class E1Point(BaseModel):
    start: int
    E1: float | None


class Compartments(BaseModel):
    resolution_used: int
    view: str
    sign_convention: str
    eigenvalues: list[ArmEigenvalue]
    balanced: bool
    region: str | None = None
    region_mean_E1: float | None = None
    region_call: str | None = None
    region_sign_consistency: float | None = None
    E1_track: list[E1Point] | None = None
    genome_A_fraction: float | None = None
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
    balanced: bool
    method: str


class ExpectedCurvePoint(BaseModel):
    dist_bp: int
    expected: float | None


class ExpectedObserved(BaseModel):
    region: str
    resolution_used: int
    ps_slope_100kb_10Mb: float | None
    expected_curve_points: list[ExpectedCurvePoint]
    balanced: bool
    oe_matrix: list[list[float | None]] | None = None
    note: str | None = None
    method: str
