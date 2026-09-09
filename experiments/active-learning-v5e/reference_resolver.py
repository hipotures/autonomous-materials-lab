"""V5e-2 automated reference identity and property resolution.

The resolver deliberately separates molecular identity evidence from thermophysical
reference evidence. PubChem is used for identity resolution, CoolProp is used as
an independent EOS/property reference when an exact identity match is available,
and NIST Chemistry WebBook is queried as an additional independent source.

A candidate is never marked ready for calibration from PubChem identity metadata
alone and the FeOS/Joback predictor under calibration is never used as reference.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from html.parser import HTMLParser
import json
import math
import re
import time
from typing import Any, Callable, Iterable
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


CAS_RE = re.compile(r"^(\d{2,7})-(\d{2})-(\d)$")
SPACE_RE = re.compile(r"\s+")


class ReferenceResolutionError(RuntimeError):
    """Base exception for reference-resolution failures."""


class SourceUnavailable(ReferenceResolutionError):
    """Raised when a remote source cannot be queried successfully."""


@dataclass(frozen=True)
class PubChemIdentity:
    cid: int
    title: str | None
    iupac_name: str | None
    canonical_smiles: str | None
    isomeric_smiles: str | None
    inchi_key: str
    molecular_formula: str | None
    molecular_weight_g_mol: float | None
    cas_numbers: tuple[str, ...]
    exact_inchi_key_match: bool
    source_url: str

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["cas_numbers"] = list(self.cas_numbers)
        return value


@dataclass(frozen=True)
class CoolPropReference:
    fluid_name: str
    cas_number: str | None
    inchi_key: str | None
    match_method: str
    version: str | None
    git_revision: str | None
    storage_temperature_k: float
    storage_pressure_pa: float
    saturation_pressure_pa: float | None
    storage_phase: str | None
    max_exit_temperature_k: float
    critical_temperature_k: float | None
    critical_pressure_pa: float | None
    acentric_factor: float | None
    normal_boiling_temperature_k: float | None
    storage_density_kg_m3: float | None
    storage_cp_j_kg_k: float | None
    latent_heat_vaporization_j_kg: float | None
    storage_enthalpy_j_kg: float | None
    exit_enthalpy_j_kg: float | None
    enthalpy_window_j_kg: float | None
    liquid_storage_feasible: bool
    required_property_complete: bool
    errors: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["errors"] = list(self.errors)
        return value


@dataclass(frozen=True)
class NistReference:
    cas_number: str
    compound_url: str
    compound_page_verified: bool
    fluid_url: str | None
    fluid_table_verified: bool
    title: str | None
    normal_boiling_temperature_k: float | None
    critical_temperature_k: float | None
    critical_pressure_pa: float | None
    latent_heat_vaporization_j_mol: float | None
    storage_temperature_k: float | None
    storage_pressure_pa: float | None
    storage_density_kg_m3: float | None
    storage_cp_j_kg_k: float | None
    storage_enthalpy_j_kg: float | None
    storage_phase: str | None
    errors: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value["errors"] = list(self.errors)
        return value


class HttpClient:
    """Small deterministic HTTP client with retries and source-friendly pacing."""

    def __init__(
        self,
        *,
        user_agent: str,
        timeout_s: float = 20.0,
        retries: int = 2,
        minimum_interval_s: float = 0.20,
        sleeper: Callable[[float], None] = time.sleep,
    ) -> None:
        self.user_agent = user_agent
        self.timeout_s = float(timeout_s)
        self.retries = int(retries)
        self.minimum_interval_s = float(minimum_interval_s)
        self._sleeper = sleeper
        self._last_request_monotonic: float | None = None

    def _pace(self) -> None:
        if self._last_request_monotonic is None:
            return
        elapsed = time.monotonic() - self._last_request_monotonic
        remaining = self.minimum_interval_s - elapsed
        if remaining > 0.0:
            self._sleeper(remaining)

    def get_bytes(self, url: str) -> bytes:
        last_error: Exception | None = None
        for attempt in range(self.retries + 1):
            self._pace()
            request = Request(
                url,
                headers={
                    "User-Agent": self.user_agent,
                    "Accept": (
                        "application/json,text/html,text/plain;q=0.9,*/*;q=0.8"
                    ),
                },
            )
            try:
                with urlopen(request, timeout=self.timeout_s) as response:
                    payload = response.read()
                self._last_request_monotonic = time.monotonic()
                return payload
            except HTTPError as exc:
                self._last_request_monotonic = time.monotonic()
                last_error = exc
                if exc.code == 404:
                    raise SourceUnavailable(f"HTTP 404 for {url}") from exc
                if 400 <= exc.code < 500 and exc.code not in (408, 429):
                    raise SourceUnavailable(
                        f"HTTP {exc.code} for {url}"
                    ) from exc
            except (URLError, TimeoutError, OSError) as exc:
                self._last_request_monotonic = time.monotonic()
                last_error = exc

            if attempt < self.retries:
                self._sleeper(min(2.0 ** attempt, 4.0))

        raise SourceUnavailable(f"request failed for {url}: {last_error}")

    def get_text(self, url: str) -> str:
        return self.get_bytes(url).decode("utf-8", errors="replace")

    def get_json(self, url: str) -> dict[str, Any]:
        try:
            return json.loads(self.get_text(url))
        except json.JSONDecodeError as exc:
            raise SourceUnavailable(f"invalid JSON from {url}") from exc


class CachedHttpClient(HttpClient):
    """HTTP client that preserves source responses for deterministic reruns."""

    def __init__(self, *, cache_dir: Any, **kwargs: Any) -> None:
        from pathlib import Path

        super().__init__(**kwargs)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def _cache_paths(self, url: str):
        digest = hashlib.sha256(url.encode("utf-8")).hexdigest()
        return (
            self.cache_dir / f"{digest}.bin",
            self.cache_dir / f"{digest}.json",
        )

    def get_bytes(self, url: str) -> bytes:
        payload_path, metadata_path = self._cache_paths(url)
        if payload_path.is_file():
            return payload_path.read_bytes()
        payload = super().get_bytes(url)
        payload_path.write_bytes(payload)
        metadata_path.write_text(
            json.dumps(
                {
                    "url": url,
                    "sha256": hashlib.sha256(payload).hexdigest(),
                    "size_bytes": len(payload),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        return payload


def _finite(value: Any) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _finite_positive(value: Any) -> float | None:
    value = _finite(value)
    if value is None or value <= 0.0:
        return None
    return value


def _safe_float(text: str | None) -> float | None:
    if text is None:
        return None
    cleaned = text.strip().replace("−", "-").replace("–", "-")
    match = re.search(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?",
        cleaned,
    )
    if not match:
        return None
    try:
        value = float(match.group(0))
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def valid_cas(value: str) -> bool:
    """Return True for syntactically valid CAS Registry Numbers."""
    match = CAS_RE.fullmatch(value.strip())
    if not match:
        return False
    digits = "".join(match.groups()[:2])
    check = int(match.group(3))
    total = sum(
        int(digit) * multiplier
        for multiplier, digit in enumerate(reversed(digits), start=1)
    )
    return total % 10 == check


def extract_cas_numbers(synonyms: Iterable[str]) -> tuple[str, ...]:
    result = {
        value.strip()
        for value in synonyms
        if isinstance(value, str) and valid_cas(value.strip())
    }
    return tuple(sorted(result))


def resolve_pubchem_identity(
    inchi_key: str,
    http: HttpClient,
) -> PubChemIdentity | None:
    """Resolve an exact PubChem compound record using InChIKey."""
    key = inchi_key.strip().upper()
    if not key:
        return None

    properties = (
        "Title,IUPACName,CanonicalSMILES,IsomericSMILES,"
        "InChIKey,MolecularFormula,MolecularWeight"
    )
    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    source_url = (
        f"{base}/inchikey/{quote(key, safe='')}/property/{properties}/JSON"
    )
    try:
        payload = http.get_json(source_url)
    except SourceUnavailable:
        return None

    rows = (
        payload.get("PropertyTable", {}).get("Properties", [])
        if isinstance(payload, dict)
        else []
    )
    if not rows:
        return None

    exact_rows = [
        row
        for row in rows
        if str(row.get("InChIKey", "")).upper() == key
    ]
    row = exact_rows[0] if exact_rows else rows[0]
    cid = row.get("CID")
    if not isinstance(cid, int):
        return None

    synonyms_url = f"{base}/cid/{cid}/synonyms/JSON"
    synonyms: list[str] = []
    try:
        synonym_payload = http.get_json(synonyms_url)
        infos = (
            synonym_payload.get("InformationList", {}).get("Information", [])
        )
        if infos:
            synonyms = list(infos[0].get("Synonym", []) or [])
    except SourceUnavailable:
        synonyms = []

    weight = _safe_float(str(row.get("MolecularWeight", "")))
    returned_key = str(row.get("InChIKey") or "").upper()
    return PubChemIdentity(
        cid=cid,
        title=(str(row["Title"]) if row.get("Title") else None),
        iupac_name=(
            str(row["IUPACName"]) if row.get("IUPACName") else None
        ),
        canonical_smiles=(
            str(row["CanonicalSMILES"])
            if row.get("CanonicalSMILES")
            else None
        ),
        isomeric_smiles=(
            str(row["IsomericSMILES"])
            if row.get("IsomericSMILES")
            else None
        ),
        inchi_key=returned_key,
        molecular_formula=(
            str(row["MolecularFormula"])
            if row.get("MolecularFormula")
            else None
        ),
        molecular_weight_g_mol=weight,
        cas_numbers=extract_cas_numbers(synonyms),
        exact_inchi_key_match=(returned_key == key),
        source_url=source_url,
    )


def _coolprop_fluid_records(cp: Any) -> list[dict[str, str | None]]:
    names = [
        item.strip()
        for item in cp.get_global_param_string("FluidsList").split(",")
        if item.strip()
    ]
    records: list[dict[str, str | None]] = []
    for name in names:

        def param(key: str) -> str | None:
            try:
                value = cp.get_fluid_param_string(name, key).strip()
            except Exception:
                return None
            return value or None

        records.append(
            {
                "name": name,
                "cas": param("CAS"),
                "inchi_key": param("INCHI_Key"),
            }
        )
    return records


def match_coolprop_fluid(
    *,
    inchi_key: str | None,
    cas_numbers: Iterable[str],
    cp: Any,
) -> tuple[str, str] | None:
    """Return (fluid name, match method) using exact persistent identifiers."""
    key = inchi_key.upper() if inchi_key else None
    cas_set = {item for item in cas_numbers if item}
    records = _coolprop_fluid_records(cp)

    if key:
        matches = [
            str(row["name"])
            for row in records
            if row.get("inchi_key")
            and str(row["inchi_key"]).upper() == key
        ]
        if len(matches) == 1:
            return matches[0], "inchi_key_exact"

    if cas_set:
        matches = [
            str(row["name"])
            for row in records
            if row.get("cas") in cas_set
        ]
        if len(matches) == 1:
            return matches[0], "cas_exact"

    return None


def _cp_scalar(cp: Any, output: str, fluid: str) -> float | None:
    try:
        return _finite(cp.PropsSI(output, fluid))
    except Exception:
        return None


def _cp_state(
    cp: Any,
    output: str,
    temperature_k: float,
    pressure_pa: float,
    fluid: str,
) -> float | None:
    try:
        return _finite(
            cp.PropsSI(
                output,
                "T",
                float(temperature_k),
                "P",
                float(pressure_pa),
                fluid,
            )
        )
    except Exception:
        return None


def _cp_saturation(
    cp: Any,
    output: str,
    *,
    temperature_k: float | None = None,
    pressure_pa: float | None = None,
    quality: int,
    fluid: str,
) -> float | None:
    try:
        if temperature_k is not None:
            value = cp.PropsSI(
                output,
                "T",
                float(temperature_k),
                "Q",
                quality,
                fluid,
            )
        elif pressure_pa is not None:
            value = cp.PropsSI(
                output,
                "P",
                float(pressure_pa),
                "Q",
                quality,
                fluid,
            )
        else:
            raise ValueError("temperature_k or pressure_pa required")
        return _finite(value)
    except Exception:
        return None


def resolve_coolprop_reference(
    *,
    inchi_key: str | None,
    cas_numbers: Iterable[str],
    cp: Any,
    storage_temperature_k: float,
    default_storage_pressure_pa: float,
    max_storage_pressure_pa: float,
    saturation_pressure_margin: float,
    max_exit_temperature_k: float,
) -> CoolPropReference | None:
    match = match_coolprop_fluid(
        inchi_key=inchi_key,
        cas_numbers=cas_numbers,
        cp=cp,
    )
    if match is None:
        return None
    fluid, match_method = match

    errors: list[str] = []
    try:
        cas = cp.get_fluid_param_string(fluid, "CAS").strip() or None
    except Exception:
        cas = None
    try:
        cp_inchi = (
            cp.get_fluid_param_string(fluid, "INCHI_Key").strip() or None
        )
    except Exception:
        cp_inchi = None
    try:
        cp_version = cp.get_global_param_string("version") or None
    except Exception:
        cp_version = None
    try:
        cp_revision = cp.get_global_param_string("gitrevision") or None
    except Exception:
        cp_revision = None

    tcrit = _cp_scalar(cp, "Tcrit", fluid)
    pcrit = _cp_scalar(cp, "Pcrit", fluid)
    acentric = _cp_scalar(cp, "acentric", fluid)

    saturation_pressure = None
    if tcrit is not None and storage_temperature_k < tcrit:
        saturation_pressure = _cp_saturation(
            cp,
            "P",
            temperature_k=storage_temperature_k,
            quality=0,
            fluid=fluid,
        )

    storage_pressure = float(default_storage_pressure_pa)
    if saturation_pressure is not None:
        storage_pressure = max(
            storage_pressure,
            saturation_pressure * float(saturation_pressure_margin),
        )

    liquid_storage_feasible = True
    if tcrit is not None and storage_temperature_k >= tcrit:
        liquid_storage_feasible = False
        errors.append("storage_temperature_at_or_above_critical")
    if storage_pressure > max_storage_pressure_pa:
        liquid_storage_feasible = False
        errors.append("required_storage_pressure_exceeds_limit")

    if liquid_storage_feasible:
        try:
            storage_phase = cp.PhaseSI(
                "T",
                float(storage_temperature_k),
                "P",
                float(storage_pressure),
                fluid,
            )
        except Exception:
            storage_phase = None
        if (
            storage_phase is not None
            and storage_phase.lower()
            not in {"liquid", "supercritical_liquid"}
        ):
            liquid_storage_feasible = False
            errors.append(f"storage_phase_not_liquid:{storage_phase}")
    else:
        storage_phase = None

    density = (
        _cp_state(
            cp,
            "Dmass",
            storage_temperature_k,
            storage_pressure,
            fluid,
        )
        if liquid_storage_feasible
        else None
    )
    cp_mass = (
        _cp_state(
            cp,
            "Cpmass",
            storage_temperature_k,
            storage_pressure,
            fluid,
        )
        if liquid_storage_feasible
        else None
    )
    h_storage = (
        _cp_state(
            cp,
            "Hmass",
            storage_temperature_k,
            storage_pressure,
            fluid,
        )
        if liquid_storage_feasible
        else None
    )
    h_exit = (
        _cp_state(
            cp,
            "Hmass",
            max_exit_temperature_k,
            storage_pressure,
            fluid,
        )
        if liquid_storage_feasible
        else None
    )
    enthalpy_window = (
        h_exit - h_storage
        if h_exit is not None and h_storage is not None
        else None
    )

    normal_boiling = _cp_saturation(
        cp,
        "T",
        pressure_pa=101325.0,
        quality=0,
        fluid=fluid,
    )
    h_vap_liq = _cp_saturation(
        cp,
        "Hmass",
        pressure_pa=101325.0,
        quality=0,
        fluid=fluid,
    )
    h_vap_vap = _cp_saturation(
        cp,
        "Hmass",
        pressure_pa=101325.0,
        quality=1,
        fluid=fluid,
    )
    latent = (
        h_vap_vap - h_vap_liq
        if h_vap_vap is not None and h_vap_liq is not None
        else None
    )

    required = {
        "critical_temperature_k": tcrit,
        "critical_pressure_pa": pcrit,
        "storage_density_kg_m3": density,
        "storage_cp_j_kg_k": cp_mass,
        "storage_enthalpy_j_kg": h_storage,
        "exit_enthalpy_j_kg": h_exit,
        "enthalpy_window_j_kg": enthalpy_window,
    }
    for name, value in required.items():
        if value is None:
            errors.append(f"missing_{name}")
    required_complete = liquid_storage_feasible and all(
        value is not None for value in required.values()
    )

    return CoolPropReference(
        fluid_name=fluid,
        cas_number=cas,
        inchi_key=cp_inchi,
        match_method=match_method,
        version=cp_version,
        git_revision=cp_revision,
        storage_temperature_k=float(storage_temperature_k),
        storage_pressure_pa=float(storage_pressure),
        saturation_pressure_pa=saturation_pressure,
        storage_phase=storage_phase,
        max_exit_temperature_k=float(max_exit_temperature_k),
        critical_temperature_k=tcrit,
        critical_pressure_pa=pcrit,
        acentric_factor=acentric,
        normal_boiling_temperature_k=normal_boiling,
        storage_density_kg_m3=density,
        storage_cp_j_kg_k=cp_mass,
        latent_heat_vaporization_j_kg=latent,
        storage_enthalpy_j_kg=h_storage,
        exit_enthalpy_j_kg=h_exit,
        enthalpy_window_j_kg=enthalpy_window,
        liquid_storage_feasible=liquid_storage_feasible,
        required_property_complete=required_complete,
        errors=tuple(errors),
    )


class _TableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self._in_cell = False
        self._cell_parts: list[str] = []
        self._row: list[str] = []
        self.rows: list[list[str]] = []
        self.title_parts: list[str] = []
        self._in_title = False

    def handle_starttag(
        self,
        tag: str,
        attrs: list[tuple[str, str | None]],
    ) -> None:
        tag = tag.lower()
        if tag in ("td", "th"):
            self._in_cell = True
            self._cell_parts = []
        elif tag == "tr":
            self._row = []
        elif tag == "title":
            self._in_title = True

    def handle_endtag(self, tag: str) -> None:
        tag = tag.lower()
        if tag in ("td", "th") and self._in_cell:
            text = SPACE_RE.sub(
                " ",
                " ".join(self._cell_parts),
            ).strip()
            self._row.append(text)
            self._cell_parts = []
            self._in_cell = False
        elif tag == "tr":
            if self._row:
                self.rows.append(self._row)
            self._row = []
        elif tag == "title":
            self._in_title = False

    def handle_data(self, data: str) -> None:
        if self._in_cell:
            self._cell_parts.append(data)
        if self._in_title:
            self.title_parts.append(data)

    @property
    def title(self) -> str | None:
        value = SPACE_RE.sub(
            " ",
            " ".join(self.title_parts),
        ).strip()
        return value or None


def _parse_nist_phase_rows(
    html: str,
) -> dict[str, tuple[float, str | None]]:
    parser = _TableParser()
    parser.feed(html)
    found: dict[str, tuple[float, str | None]] = {}
    label_map = {
        "tboil": "normal_boiling_temperature_k",
        "boiling point": "normal_boiling_temperature_k",
        "tc": "critical_temperature_k",
        "critical temperature": "critical_temperature_k",
        "pc": "critical_pressure",
        "critical pressure": "critical_pressure",
        "δvaph": "latent_heat_vaporization",
        "enthalpy of vaporization": "latent_heat_vaporization",
    }
    for row in parser.rows:
        if len(row) < 2:
            continue
        label = SPACE_RE.sub(" ", row[0]).strip().lower()
        mapped = None
        for candidate, canonical in label_map.items():
            if candidate.lower() == label or candidate.lower() in label:
                mapped = canonical
                break
        if mapped is None:
            continue
        value = _safe_float(row[1])
        if value is None:
            continue
        unit = row[2] if len(row) >= 3 else None
        found.setdefault(mapped, (value, unit))
    return found


def _convert_pressure_to_pa(
    value: float | None,
    unit: str | None,
) -> float | None:
    if value is None:
        return None
    normalized = (unit or "").strip().lower().replace(" ", "")
    if normalized == "pa":
        return value
    if normalized == "kpa":
        return value * 1.0e3
    if normalized == "mpa":
        return value * 1.0e6
    if normalized == "bar":
        return value * 1.0e5
    if normalized == "atm":
        return value * 101325.0
    return None


def _convert_vap_h_to_j_mol(
    value: float | None,
    unit: str | None,
) -> float | None:
    if value is None:
        return None
    normalized = (unit or "").strip().lower().replace(" ", "")
    if normalized in ("j/mol", "jmol-1", "jmol−1"):
        return value
    if normalized in ("kj/mol", "kjmol-1", "kjmol−1"):
        return value * 1.0e3
    return None


def _parse_nist_fluid_table(html: str) -> dict[str, Any]:
    parser = _TableParser()
    parser.feed(html)
    rows = parser.rows
    header_index = None
    header: list[str] = []
    for index, row in enumerate(rows):
        lowered = [cell.lower() for cell in row]
        if any("temperature" in cell for cell in lowered) and any(
            "density" in cell for cell in lowered
        ):
            header_index = index
            header = row
            break
    if header_index is None:
        return {}

    data_row = None
    for row in rows[header_index + 1 :]:
        if (
            len(row) >= min(4, len(header))
            and _safe_float(row[0]) is not None
        ):
            data_row = row
            break
    if data_row is None:
        return {}

    output: dict[str, Any] = {}
    for index, label in enumerate(header):
        if index >= len(data_row):
            continue
        key = SPACE_RE.sub(" ", label).strip().lower()
        value_text = data_row[index]
        value = _safe_float(value_text)
        if "temperature" in key and value is not None:
            output["storage_temperature_k"] = value
        elif "pressure" in key and value is not None:
            output["storage_pressure_pa"] = value * 1.0e6
        elif "density" in key and value is not None:
            output["storage_density_kg_m3"] = value
        elif (
            key.startswith("cp")
            or "c p" in key
            or "heat capacity at constant pressure" in key
        ) and value is not None:
            output["storage_cp_j_kg_k"] = value * 1.0e3
        elif "enthalpy" in key and value is not None:
            output["storage_enthalpy_j_kg"] = value * 1.0e3
        elif "phase" in key:
            output["storage_phase"] = value_text.strip() or None
    return output


def _nist_id(cas_number: str) -> str:
    return "C" + cas_number.replace("-", "")


def resolve_nist_reference(
    *,
    cas_number: str,
    http: HttpClient,
    storage_temperature_k: float,
    storage_pressure_pa: float,
) -> NistReference | None:
    if not valid_cas(cas_number):
        return None
    nist_id = _nist_id(cas_number)
    compound_url = (
        "https://webbook.nist.gov/cgi/cbook.cgi?"
        + urlencode(
            {
                "ID": nist_id,
                "Units": "SI",
                "Mask": "4",
            }
        )
    )
    errors: list[str] = []
    try:
        html = http.get_text(compound_url)
    except SourceUnavailable:
        return None

    lowered = html.lower()
    compound_verified = (
        "nist chemistry webbook" in lowered
        and "name not found" not in lowered
        and cas_number in html
    )
    parser = _TableParser()
    parser.feed(html)
    phase = _parse_nist_phase_rows(html)

    def phase_value(name: str) -> tuple[float | None, str | None]:
        value = phase.get(name)
        if value is None:
            return None, None
        return value

    tboil, _ = phase_value("normal_boiling_temperature_k")
    tcrit, _ = phase_value("critical_temperature_k")
    pc_value, pc_unit = phase_value("critical_pressure")
    hvap_value, hvap_unit = phase_value("latent_heat_vaporization")
    pcrit = _convert_pressure_to_pa(pc_value, pc_unit)
    hvap = _convert_vap_h_to_j_mol(hvap_value, hvap_unit)

    fluid_params = {
        "P": f"{storage_pressure_pa / 1.0e6:.10g}",
        "TLow": f"{storage_temperature_k:.10g}",
        "THigh": f"{storage_temperature_k + 0.01:.10g}",
        "TInc": "1",
        "Digits": "8",
        "ID": nist_id,
        "Action": "Load",
        "Type": "IsoBar",
        "TUnit": "K",
        "PUnit": "MPa",
        "DUnit": "kg/m3",
        "HUnit": "kJ/kg",
        "WUnit": "m/s",
        "VisUnit": "Pa*s",
        "STUnit": "N/m",
        "RefState": "DEF",
    }
    fluid_url = (
        "https://webbook.nist.gov/cgi/fluid.cgi?"
        + urlencode(fluid_params)
    )
    fluid_values: dict[str, Any] = {}
    fluid_verified = False
    try:
        fluid_html = http.get_text(fluid_url)
        fluid_values = _parse_nist_fluid_table(fluid_html)
        fluid_verified = bool(
            fluid_values.get("storage_density_kg_m3")
        )
        if not fluid_verified:
            errors.append("nist_fluid_table_not_available")
    except SourceUnavailable:
        errors.append("nist_fluid_endpoint_unavailable")

    return NistReference(
        cas_number=cas_number,
        compound_url=compound_url,
        compound_page_verified=compound_verified,
        fluid_url=fluid_url,
        fluid_table_verified=fluid_verified,
        title=parser.title,
        normal_boiling_temperature_k=tboil,
        critical_temperature_k=tcrit,
        critical_pressure_pa=pcrit,
        latent_heat_vaporization_j_mol=hvap,
        storage_temperature_k=fluid_values.get(
            "storage_temperature_k"
        ),
        storage_pressure_pa=fluid_values.get(
            "storage_pressure_pa"
        ),
        storage_density_kg_m3=fluid_values.get(
            "storage_density_kg_m3"
        ),
        storage_cp_j_kg_k=fluid_values.get(
            "storage_cp_j_kg_k"
        ),
        storage_enthalpy_j_kg=fluid_values.get(
            "storage_enthalpy_j_kg"
        ),
        storage_phase=fluid_values.get("storage_phase"),
        errors=tuple(errors),
    )


def identity_consistency(
    *,
    expected_inchi_key: str | None,
    expected_formula: str | None,
    expected_molecular_weight_g_mol: float | None,
    pubchem: PubChemIdentity | None,
    max_molecular_weight_relative_error: float,
) -> dict[str, Any]:
    if pubchem is None:
        return {
            "identity_verified": False,
            "checks": {
                "pubchem_record_found": False,
                "inchi_key_exact": False,
                "formula_exact": False,
                "molecular_weight_consistent": False,
            },
            "reasons": ["pubchem_identity_not_resolved"],
        }

    expected_key = (
        expected_inchi_key.upper()
        if expected_inchi_key
        else None
    )
    inchi_exact = bool(
        expected_key
        and pubchem.inchi_key.upper() == expected_key
        and pubchem.exact_inchi_key_match
    )
    formula_exact = bool(
        expected_formula
        and pubchem.molecular_formula
        and expected_formula == pubchem.molecular_formula
    )
    weight_consistent = False
    relative_error = None
    if (
        expected_molecular_weight_g_mol is not None
        and pubchem.molecular_weight_g_mol is not None
        and expected_molecular_weight_g_mol > 0.0
    ):
        relative_error = abs(
            pubchem.molecular_weight_g_mol
            - expected_molecular_weight_g_mol
        ) / expected_molecular_weight_g_mol
        weight_consistent = (
            relative_error
            <= max_molecular_weight_relative_error
        )

    reasons: list[str] = []
    if not inchi_exact:
        reasons.append("inchi_key_mismatch")
    if expected_formula and not formula_exact:
        reasons.append("molecular_formula_mismatch")
    if (
        expected_molecular_weight_g_mol is not None
        and not weight_consistent
    ):
        reasons.append("molecular_weight_mismatch")

    verified = inchi_exact and (
        (not expected_formula or formula_exact)
        and (
            expected_molecular_weight_g_mol is None
            or weight_consistent
        )
    )
    return {
        "identity_verified": verified,
        "checks": {
            "pubchem_record_found": True,
            "inchi_key_exact": inchi_exact,
            "formula_exact": formula_exact,
            "molecular_weight_consistent": weight_consistent,
        },
        "molecular_weight_relative_error": relative_error,
        "reasons": reasons,
    }


def classify_resolution(
    *,
    identity_verified: bool,
    coolprop: CoolPropReference | None,
    nist: NistReference | None,
) -> tuple[str, str, bool]:
    """Return (status, quality, usable_for_v5b_calibration)."""
    if not identity_verified:
        return "unresolved_identity", "insufficient", False

    if coolprop is not None and coolprop.required_property_complete:
        return "resolved", "high", True

    if nist is not None and nist.fluid_table_verified:
        return "resolved_partial", "medium", False
    if nist is not None and nist.compound_page_verified:
        return "resolved_identity_only", "low", False
    if coolprop is not None:
        return "resolved_property_incomplete", "medium", False
    return "resolved_identity_only", "low", False


def disagreement_metrics(
    coolprop: CoolPropReference | None,
    nist: NistReference | None,
) -> dict[str, Any]:
    if coolprop is None or nist is None:
        return {}

    pairs = {
        "critical_temperature": (
            coolprop.critical_temperature_k,
            nist.critical_temperature_k,
        ),
        "critical_pressure": (
            coolprop.critical_pressure_pa,
            nist.critical_pressure_pa,
        ),
        "normal_boiling_temperature": (
            coolprop.normal_boiling_temperature_k,
            nist.normal_boiling_temperature_k,
        ),
        "storage_density": (
            coolprop.storage_density_kg_m3,
            nist.storage_density_kg_m3,
        ),
        "storage_cp": (
            coolprop.storage_cp_j_kg_k,
            nist.storage_cp_j_kg_k,
        ),
    }
    output: dict[str, Any] = {}
    for name, (left, right) in pairs.items():
        if left is None or right is None:
            continue
        denom = max(abs(left), abs(right), 1.0e-30)
        output[name] = {
            "coolprop": left,
            "nist": right,
            "relative_difference": abs(left - right) / denom,
        }
    return output


def resolve_target(
    target: dict[str, Any],
    *,
    http: HttpClient,
    cp: Any,
    config: dict[str, Any],
) -> dict[str, Any]:
    identity_cfg = config["identity"]
    properties_cfg = config["properties"]

    expected_key = target.get("inchi_key")
    pubchem = (
        resolve_pubchem_identity(str(expected_key), http)
        if expected_key
        else None
    )
    identity = identity_consistency(
        expected_inchi_key=(
            str(expected_key)
            if expected_key
            else None
        ),
        expected_formula=(
            str(target["molecular_formula"])
            if target.get("molecular_formula")
            else None
        ),
        expected_molecular_weight_g_mol=_finite_positive(
            target.get("molecular_weight_g_mol")
        ),
        pubchem=pubchem,
        max_molecular_weight_relative_error=float(
            identity_cfg[
                "max_molecular_weight_relative_error"
            ]
        ),
    )

    cas_numbers = pubchem.cas_numbers if pubchem else tuple()
    coolprop = None
    if identity["identity_verified"]:
        coolprop = resolve_coolprop_reference(
            inchi_key=(
                pubchem.inchi_key
                if pubchem
                else expected_key
            ),
            cas_numbers=cas_numbers,
            cp=cp,
            storage_temperature_k=float(
                properties_cfg["storage_temperature_k"]
            ),
            default_storage_pressure_pa=float(
                properties_cfg[
                    "default_storage_pressure_pa"
                ]
            ),
            max_storage_pressure_pa=float(
                properties_cfg["max_storage_pressure_pa"]
            ),
            saturation_pressure_margin=float(
                properties_cfg[
                    "saturation_pressure_margin"
                ]
            ),
            max_exit_temperature_k=float(
                properties_cfg["max_exit_temperature_k"]
            ),
        )

    nist = None
    if identity["identity_verified"] and cas_numbers:
        ordered = list(cas_numbers)
        if (
            coolprop
            and coolprop.cas_number in ordered
        ):
            ordered.remove(coolprop.cas_number)
            ordered.insert(0, coolprop.cas_number)
        for cas in ordered:
            nist = resolve_nist_reference(
                cas_number=cas,
                http=http,
                storage_temperature_k=(
                    coolprop.storage_temperature_k
                    if coolprop
                    else float(
                        properties_cfg[
                            "storage_temperature_k"
                        ]
                    )
                ),
                storage_pressure_pa=(
                    coolprop.storage_pressure_pa
                    if coolprop
                    else float(
                        properties_cfg[
                            "default_storage_pressure_pa"
                        ]
                    )
                ),
            )
            if (
                nist is not None
                and nist.compound_page_verified
            ):
                break

    status, quality, usable = classify_resolution(
        identity_verified=bool(identity["identity_verified"]),
        coolprop=coolprop,
        nist=nist,
    )
    disagreements = disagreement_metrics(coolprop, nist)
    max_disagreement = float(
        config["quality"][
            "maximum_cross_source_relative_difference"
        ]
    )
    conflict_fields = [
        key
        for key, value in disagreements.items()
        if value["relative_difference"] > max_disagreement
    ]
    if conflict_fields and usable:
        status = "resolved_conflict_review"
        quality = "medium"
        usable = False

    preferred_name = None
    if pubchem is not None:
        preferred_name = (
            pubchem.title
            or pubchem.iupac_name
        )

    return {
        "acquisition_rank": target.get("acquisition_rank"),
        "candidate_id": target.get("candidate_id"),
        "smiles": target.get("smiles"),
        "canonical_smiles": target.get(
            "canonical_smiles"
        ),
        "expected_inchi_key": expected_key,
        "expected_molecular_formula": target.get(
            "molecular_formula"
        ),
        "expected_molecular_weight_g_mol": target.get(
            "molecular_weight_g_mol"
        ),
        "family": target.get("family"),
        "resolution_status": status,
        "reference_quality": quality,
        "usable_for_v5b_calibration": usable,
        "preferred_name": preferred_name,
        "identity": identity,
        "pubchem": (
            pubchem.to_dict()
            if pubchem
            else None
        ),
        "coolprop": (
            coolprop.to_dict()
            if coolprop
            else None
        ),
        "nist": (
            nist.to_dict()
            if nist
            else None
        ),
        "cross_source_disagreement": disagreements,
        "conflict_fields": conflict_fields,
        "reference_independence": {
            "predictor_under_calibration": (
                "FeOS GC-PC-SAFT + Joback"
            ),
            "accepted_calibration_backend": (
                "CoolProp"
                if usable and coolprop is not None
                else None
            ),
            "predictor_reused_as_reference": False,
        },
    }
