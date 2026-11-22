"""The data sources Stargazer knows about, how to cite them, and a probe of which ones this machine can reach.

`SOURCES` lists every host the project has considered (the list began as spikes/first-light/sources.txt). Only
some have connectors in this package; those carry full citation or acknowledgement text, which every
`Provenance` built from them repeats.
"""
from __future__ import annotations

import concurrent.futures as futures
import logging
import urllib.error
import urllib.request
from dataclasses import dataclass

log = logging.getLogger(__name__)

USER_AGENT = "stargazer/0.1"

ROUTES = ("cloud", "api", "code")
STATUSES = ("reachable", "blocked", "unreachable")

ZTF_CITATION = (
    "Zwicky Transient Facility Data Release 24 light curves, NASA/IPAC Infrared Science Archive, "
    "doi:10.26131/IRSA598. Based on observations obtained with the Samuel Oschin Telescope 48-inch and the "
    "60-inch Telescope at the Palomar Observatory as part of the Zwicky Transient Facility project "
    "(Bellm et al. 2019, PASP 131, 018002; Masci et al. 2019, PASP 131, 018003)."
)
GAIA_CITATION = (
    "This work has made use of data from the European Space Agency (ESA) mission Gaia "
    "(https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, "
    "https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national "
    "institutions, in particular the institutions participating in the Gaia Multilateral Agreement. "
    "Gaia Collaboration, Vallenari et al. 2023, A&A 674, A1 (Gaia DR3). Read from the HATS copy at "
    "s3://stpubdata/gaia/gaia_dr3/public/hats (STScI/MAST)."
)
KEPLER_CITATION = (
    "This paper includes data collected by the Kepler mission and obtained from the MAST data archive at the "
    "Space Telescope Science Institute (STScI). Funding for the Kepler mission is provided by the NASA Science "
    "Mission Directorate. STScI is operated by the Association of Universities for Research in Astronomy, Inc., "
    "under NASA contract NAS 5-26555. Kepler long-cadence light curves, doi:10.17909/T9488N."
)
FERMI_CITATION = (
    "Fermi Large Area Telescope Light Curve Repository (Fermi-LAT Collaboration, Abdollahi et al. 2023, "
    "ApJS 265, 31, arXiv:2301.01607), served by the Fermi Science Support Center, NASA GSFC."
)


@dataclass(frozen=True)
class SourceInfo:
    """One data source: which host serves it, what it holds, and how to acknowledge it.

    route is "cloud" (public object storage read directly), "api" (an archive web service) or "code"
    (software and package hosts).
    """

    id: str
    title: str
    host: str
    serves: str
    route: str
    citation: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if self.route not in ROUTES:
            raise ValueError(f"route must be one of {ROUTES}, not {self.route!r}")


@dataclass
class Reachability:
    """Result of probing one source's host: status is "reachable", "blocked" or "unreachable"."""

    id: str
    host: str
    status: str
    detail: str


def _info(id: str, title: str, host: str, serves: str, route: str, citation: str = "", notes: str = "") -> SourceInfo:
    return SourceInfo(id, title, host, serves, route, citation, notes)


_ENTRIES = [
    # Archive web services and APIs
    _info("mast", "MAST", "mast.stsci.edu", "MAST portal and API (Kepler, K2, TESS, Pan-STARRS, HST)", "api",
          KEPLER_CITATION, "Mast.Caom.Cone via /api/v0/invoke works anonymously."),
    _info("mast_archive", "MAST legacy archive", "archive.stsci.edu", "MAST legacy archive and TESS bulk downloads",
          "api"),
    _info("mast_catalogs", "MAST catalogue API", "catalogs.mast.stsci.edu", "MAST catalogue API (Pan-STARRS, TIC)",
          "api", notes="Refused by the sandbox proxy on 21 Nov 2026."),
    _info("ztf_irsa", "IRSA ZTF light-curve service", "irsa.ipac.caltech.edu",
          "IRSA (ZTF light-curve API, WISE/NEOWISE, TAP)", "api", ZTF_CITATION,
          "Used by stargazer.sources.ztf with route='irsa'."),
    _info("ned", "NASA/IPAC Extragalactic Database", "ned.ipac.caltech.edu", "NASA/IPAC Extragalactic Database", "api"),
    _info("gaia_archive", "ESA Gaia archive", "gea.esac.esa.int", "ESA Gaia archive (TAP)", "api", GAIA_CITATION),
    _info("simbad", "SIMBAD", "simbad.cds.unistra.fr", "SIMBAD object database", "api"),
    _info("vizier", "VizieR", "vizier.cds.unistra.fr", "VizieR catalogue service", "api"),
    _info("heasarc", "HEASARC", "heasarc.gsfc.nasa.gov", "HEASARC (Swift, NICER, MAXI, Fermi GBM)", "api"),
    _info("fermi_lcr", "Fermi-LAT Light Curve Repository", "fermi.gsfc.nasa.gov",
          "Fermi Science Support Center (LAT Light Curve Repository)", "api", FERMI_CITATION,
          "Undocumented JSON endpoint queryDB.php; used by stargazer.sources.fermi."),
    _info("sdss_skyserver", "SDSS SkyServer", "skyserver.sdss.org", "SDSS SkyServer SQL", "api"),
    _info("sdss_data", "SDSS science archive", "data.sdss.org", "SDSS science archive files", "api"),
    _info("desi", "DESI public data", "data.desi.lbl.gov", "DESI public data releases", "api"),
    _info("alerce", "ALeRCE", "api.alerce.online", "ALeRCE alert broker API", "api"),
    _info("fink", "Fink", "api.fink-portal.org", "Fink alert broker API", "api"),
    _info("lasair", "Lasair", "lasair-ztf.lsst.ac.uk", "Lasair alert broker", "api"),
    _info("antares", "ANTARES", "antares.noirlab.edu", "ANTARES alert broker", "api"),
    _info("datalab", "NOIRLab Astro Data Lab", "datalab.noirlab.edu", "NOIRLab Astro Data Lab", "api"),
    _info("rubin", "Rubin Science Platform", "data.lsst.cloud", "Rubin Science Platform", "api"),
    _info("bl_portal", "Breakthrough Listen open data portal", "seti.berkeley.edu",
          "Breakthrough Listen open data portal", "api"),
    _info("bl_search", "Breakthrough Listen open data search", "breakthroughinitiatives.org",
          "Breakthrough Listen open data search", "api"),
    _info("bl_data", "Breakthrough Listen data server", "bldata.berkeley.edu", "Breakthrough Listen data server",
          "api"),
    _info("atnf", "ATNF Pulsar Catalogue", "www.atnf.csiro.au", "ATNF Pulsar Catalogue", "api"),
    _info("csiro_dap", "CSIRO Data Access Portal", "data.csiro.au", "CSIRO Data Access Portal (Parkes pulsar data)",
          "api"),
    _info("chime_frb", "CHIME/FRB", "www.chime-frb.ca", "CHIME/FRB catalogues", "api"),
    _info("tns", "Transient Name Server", "www.wis-tns.org", "Transient Name Server", "api"),
    _info("gwosc", "GWOSC", "gwosc.org", "Gravitational Wave Open Science Center", "api"),
    _info("icecube", "IceCube", "icecube.wisc.edu", "IceCube public data releases", "api"),
    _info("asassn", "ASAS-SN Sky Patrol", "asas-sn.ifa.hawaii.edu", "ASAS-SN Sky Patrol", "api"),
    _info("atlas", "ATLAS forced photometry", "fallingstar-data.com", "ATLAS forced-photometry server", "api"),
    _info("ogle", "OGLE", "ogle.astrouw.edu.pl", "OGLE photometry archive", "api"),
    _info("nanograv", "NANOGrav", "nanograv.org", "NANOGrav data releases", "api"),
    _info("ads", "NASA ADS", "api.adsabs.harvard.edu", "NASA ADS literature API", "api"),
    _info("arxiv", "arXiv", "arxiv.org", "arXiv preprints", "api"),
    _info("zenodo", "Zenodo", "zenodo.org", "Zenodo data deposits", "api"),
    _info("lsdb", "LSDB HATS catalogues", "data.lsdb.io", "LSDB-hosted HATS catalogues", "api"),
    # Public cloud copies of the same data
    _info("stpubdata", "STScI public bucket", "stpubdata.s3.amazonaws.com",
          "STScI public bucket (TESS, Kepler, K2, Gaia DR3, Pan-STARRS, SDSS DR7 spectra, HST, JWST)", "cloud",
          KEPLER_CITATION, "Kepler light curves and the Gaia DR3 HATS catalogue are read from here."),
    _info("ztf_dr24_lc", "ZTF DR24 light curves (HATS)", "ipac-irsa-ztf.s3.amazonaws.com",
          "IRSA ZTF bucket (DR23 and DR24 light curves and objects, HATS format)", "cloud", ZTF_CITATION,
          "ztf/enhanced/dr24/lc/hats/ztf_dr24_lc-hats; used by stargazer.sources.ztf (route='cloud')."),
    _info("wise_s3", "IRSA WISE bucket", "nasa-irsa-wise.s3.amazonaws.com",
          "IRSA WISE bucket (AllWISE, NEOWISE single exposures, unWISE)", "cloud"),
    _info("spherex_s3", "IRSA SPHEREx bucket", "nasa-irsa-spherex.s3.amazonaws.com", "IRSA SPHEREx bucket", "cloud"),
    _info("heasarc_s3", "HEASARC bucket", "nasa-heasarc.s3.amazonaws.com",
          "HEASARC bucket (Fermi, Swift, NICER, MAXI, RXTE, XMM, Chandra and more)", "cloud"),
    _info("gcs", "Google Cloud Storage", "storage.googleapis.com",
          "Google Cloud Storage (Breakthrough Listen sample files)", "cloud"),
    # Code and packages
    _info("github", "GitHub", "github.com", "Source repositories", "code"),
    _info("github_raw", "GitHub raw files", "raw.githubusercontent.com", "Raw files from repositories", "code"),
    _info("pypi", "PyPI", "pypi.org", "Python packages", "code"),
    _info("conda", "Anaconda", "conda.anaconda.org", "Conda packages", "code"),
]

SOURCES: dict[str, SourceInfo] = {entry.id: entry for entry in _ENTRIES}

# Connector sources that are not hosts of their own but datasets on a host above.
SOURCES["gaia_dr3_hats"] = _info(
    "gaia_dr3_hats", "Gaia DR3 source table (HATS)", "stpubdata.s3.amazonaws.com",
    "Gaia DR3 gaia_source as HATS Parquet at gaia/gaia_dr3/public/hats/gaia", "cloud", GAIA_CITATION,
    "Used by stargazer.sources.gaia.")
SOURCES["kepler_lc"] = _info(
    "kepler_lc", "Kepler long-cadence light curves", "stpubdata.s3.amazonaws.com",
    "Kepler light-curve FITS files at kepler/public/lightcurves", "cloud", KEPLER_CITATION,
    "Used by stargazer.sources.kepler.")


def probe(host: str, timeout: float = 20.0) -> tuple[str, str]:
    """GET https://<host>/ once and classify the outcome as (status, detail).

    Any HTTP status from the server means "reachable"; a refused proxy tunnel (or HTTP 407) means "blocked";
    a timeout, DNS or TLS failure means "unreachable".
    """
    request = urllib.request.Request(f"https://{host}/", method="GET", headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            return "reachable", f"HTTP {response.status}"
    except urllib.error.HTTPError as error:  # the server itself answered
        return "reachable", f"HTTP {error.code}"
    except urllib.error.URLError as error:
        reason = str(error.reason)
        if "Tunnel connection failed" in reason or "407" in reason:
            return "blocked", reason
        return "unreachable", reason
    except (OSError, ValueError) as error:  # timeouts, TLS failures, malformed responses
        return "unreachable", f"{type(error).__name__}: {error}"


def check(ids: list[str] | None = None, timeout: float = 20.0, workers: int = 4) -> list[Reachability]:
    """Probe the hosts of the given source ids (all sources when None), a few at a time.

    Returns one Reachability per id, in the order given. Each distinct host is probed only once.
    """
    chosen = list(SOURCES) if ids is None else list(ids)
    unknown = [i for i in chosen if i not in SOURCES]
    if unknown:
        raise KeyError(f"unknown source ids: {unknown}")
    hosts = sorted({SOURCES[i].host for i in chosen})
    with futures.ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        outcome = dict(zip(hosts, pool.map(lambda h: probe(h, timeout), hosts)))
    results = [Reachability(i, SOURCES[i].host, *outcome[SOURCES[i].host]) for i in chosen]
    for r in results:
        log.debug("%s %s %s %s", r.id, r.host, r.status, r.detail)
    return results
