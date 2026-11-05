"""Download the 10x Genomics human tonsil dataset used by the app (~70 MB).

Dataset: "CytAssist FFPE Protein Expression – Human Tonsil" (Space Ranger 2.1.0), public from 10x Genomics.
Only the two processed output files the app needs are downloaded:
    filtered_feature_bc_matrix.h5   RNA + protein counts per spot
    spatial/                        H&E images, spot positions and scale factors

Run:  python download_data.py
"""
import tarfile
import urllib.request
from pathlib import Path

SAMPLE = "CytAssist_FFPE_Protein_Expression_Human_Tonsil"
BASE_URL = f"https://cf.10xgenomics.com/samples/spatial-exp/2.1.0/{SAMPLE}/{SAMPLE}"
OUT_DIR = Path(__file__).parent / "data" / "Human_Tonsil_FFPE"


def download(url, dest):
    """Download url to dest. Writes to a .part file first, so a failed download leaves no broken file."""
    print(f"Downloading {url.rsplit('/', 1)[-1]} …")
    # The 10x CDN rejects Python's default User-Agent (HTTP 403), so send a generic browser one.
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req) as resp, open(tmp, "wb") as f:
        while chunk := resp.read(1 << 20):
            f.write(chunk)
    tmp.rename(dest)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    h5 = OUT_DIR / "filtered_feature_bc_matrix.h5"
    if not h5.exists():
        download(f"{BASE_URL}_filtered_feature_bc_matrix.h5", h5)

    if not (OUT_DIR / "spatial" / "scalefactors_json.json").exists():
        archive = OUT_DIR / "spatial.tar.gz"
        download(f"{BASE_URL}_spatial.tar.gz", archive)
        with tarfile.open(archive) as tar:
            # Only extract the expected spatial/ files (guards against unsafe paths in the archive)
            members = [m for m in tar.getmembers() if m.name.startswith("spatial/") and ".." not in m.name]
            tar.extractall(OUT_DIR, members=members, filter="data")
        archive.unlink()

    print(f"Done. Data is in {OUT_DIR}")


if __name__ == "__main__":
    main()
