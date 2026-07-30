"""Archive a tagged release to Zenodo as a new version of the concept record.

Creates a new version under the existing concept DOI rather than starting a
fresh lineage, so ``10.5281/zenodo.20357011`` keeps resolving to the newest
release and the DOI cited by CITATION.cff, README and the paper stays valid.

Metadata is read from CITATION.cff so it is maintained in exactly one place.

Configuration comes from the environment, since this runs in CI:

* ``ZENODO_TOKEN``    — personal access token, scopes ``deposit:write`` and
  ``deposit:actions``.  Required.
* ``RELEASE_TAG``     — tag to archive, e.g. ``v0.4.0``.  Defaults to
  ``GITHUB_REF_NAME``.
* ``ZENODO_PUBLISH``  — ``"false"`` leaves the deposit as an unpublished
  draft for inspection.  Anything else publishes, minting the DOI.
* ``ZENODO_DRY_RUN``  — ``"true"`` resolves the record and builds the
  metadata, then stops without creating anything.

Publishing mints a DOI and cannot be undone, so exercise a dry run first when
changing this script.
"""

import json
import os
import pathlib
import sys
import urllib.error
import urllib.request

import yaml

# --- configuration -------------------------------------------------------
ZENODO_API = "https://zenodo.org/api"
CONCEPT_RECID = "20357011"  # parent of the opf-potpourri version series
GITHUB_REPO = "RWTH-IAEW/opf-potpourri"
CITATION_FILE = "CITATION.cff"
UPLOAD_TYPE = "software"
# -------------------------------------------------------------------------

TOKEN = os.environ.get("ZENODO_TOKEN", "")
TAG = os.environ.get("RELEASE_TAG") or os.environ.get("GITHUB_REF_NAME", "")
PUBLISH = os.environ.get("ZENODO_PUBLISH", "true").lower() != "false"
DRY_RUN = os.environ.get("ZENODO_DRY_RUN", "false").lower() == "true"


def log(msg):
    print(msg, flush=True)


def api(method, url, *, data=None, headers=None, raw=None):
    """Call the Zenodo API and return the decoded JSON body (or None)."""
    hdrs = {"Authorization": f"Bearer {TOKEN}"}
    if headers:
        hdrs.update(headers)
    body = raw
    if data is not None:
        body = json.dumps(data).encode()
        hdrs["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=body, headers=hdrs, method=method)
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            payload = r.read()
            return json.loads(payload) if payload else None
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode(errors="replace")[:800]
        raise SystemExit(
            f"Zenodo API {method} {url} failed: HTTP {exc.code}\n{detail}"
        ) from exc


def citation_metadata(version):
    """Build Zenodo metadata from CITATION.cff."""
    cff = yaml.safe_load(pathlib.Path(CITATION_FILE).read_text())

    creators = []
    for a in cff.get("authors", []):
        name = " ".join(
            p for p in (a.get("family-names"), a.get("given-names")) if p
        )
        if a.get("family-names") and a.get("given-names"):
            name = f"{a['family-names']}, {a['given-names']}"
        entry = {"name": name or a.get("name", "Unknown")}
        if a.get("affiliation"):
            entry["affiliation"] = " ".join(a["affiliation"].split())
        if a.get("orcid"):
            entry["orcid"] = str(a["orcid"]).rsplit("/", 1)[-1]
        creators.append(entry)

    meta = {
        "title": cff.get("title", GITHUB_REPO),
        "upload_type": UPLOAD_TYPE,
        "description": " ".join(str(cff.get("abstract", "")).split()),
        "creators": creators,
        "version": version,
        "publication_date": str(cff.get("date-released", "")),
        "access_right": "open",
        "related_identifiers": [
            {
                "identifier": (f"https://github.com/{GITHUB_REPO}/tree/{TAG}"),
                "relation": "isSupplementTo",
                "scheme": "url",
            }
        ],
    }
    lic = cff.get("license")
    if lic:
        meta["license"] = str(lic).lower()
    kw = cff.get("keywords")
    if kw:
        meta["keywords"] = [str(k) for k in kw]
    return meta


def latest_version_id():
    """Resolve the newest published record in the concept series."""
    # Resolve via search rather than GET /records/{conceptrecid}: that
    # endpoint returns HTTP 410 "record has been deleted" for a concept id,
    # apparently because the concept is a placeholder rather than a real
    # record. Searching on conceptrecid is answered consistently.
    hits = api(
        "GET",
        f"{ZENODO_API}/records?q=conceptrecid:{CONCEPT_RECID}"
        f"&all_versions=true&size=100",
    )
    records = ((hits or {}).get("hits") or {}).get("hits") or []
    if not records:
        raise SystemExit(
            f"no published records found for concept {CONCEPT_RECID}"
        )

    # Any member of the series knows the current head; ask it rather than
    # guessing from ids or dates.
    latest_link = (records[0].get("links") or {}).get("latest")
    if latest_link:
        newest = api("GET", latest_link)
    else:
        newest = max(records, key=lambda r: int(r["id"]))
    return str(newest["id"]), (newest.get("metadata") or {}).get("version")


def main():
    if not TOKEN:
        raise SystemExit("ZENODO_TOKEN is not set")
    if not TAG:
        raise SystemExit("RELEASE_TAG / GITHUB_REF_NAME is not set")

    # Zenodo stores the existing series as "v0.3.1", so keep the tag
    # form for the version field and use the bare number in filenames.
    version = TAG
    file_version = TAG.lstrip("v")
    log(f"tag={TAG} version={version} concept={CONCEPT_RECID}")

    latest_id, latest_version = latest_version_id()
    log(f"latest published version: {latest_version} (record {latest_id})")
    if latest_version == version:
        raise SystemExit(
            f"version {version} is already archived on Zenodo; refusing to "
            f"create a duplicate"
        )

    meta = citation_metadata(version)
    log(
        f"metadata: {len(meta['creators'])} creators, "
        f"license={meta.get('license')}, date={meta['publication_date']}"
    )

    archive_url = (
        f"https://github.com/{GITHUB_REPO}/archive/refs/tags/{TAG}.zip"
    )
    filename = f"opf-potpourri-{file_version}.zip"

    if DRY_RUN:
        log("dry run: resolved everything, creating nothing")
        log(json.dumps(meta, indent=2)[:1500])
        return

    log("creating new version draft ...")
    newv = api(
        "POST",
        f"{ZENODO_API}/deposit/depositions/{latest_id}/actions/newversion",
    )
    draft_url = newv["links"]["latest_draft"]
    draft = api("GET", draft_url)
    draft_id = draft["id"]
    log(f"draft deposition {draft_id}")

    # A new version inherits the previous version's files; drop them so the
    # deposit contains only this release's archive.
    for f in draft.get("files", []):
        api(
            "DELETE",
            f"{ZENODO_API}/deposit/depositions/{draft_id}/files/{f['id']}",
        )
        log(f"  removed inherited file {f.get('filename')}")

    log(f"downloading {archive_url}")
    with urllib.request.urlopen(archive_url, timeout=300) as r:
        blob = r.read()
    log(f"  {len(blob)} bytes")

    bucket = draft["links"]["bucket"]
    log(f"uploading {filename}")
    api(
        "PUT",
        f"{bucket}/{filename}",
        raw=blob,
        headers={"Content-Type": "application/octet-stream"},
    )

    log("setting metadata")
    api(
        "PUT",
        f"{ZENODO_API}/deposit/depositions/{draft_id}",
        data={"metadata": meta},
    )

    if not PUBLISH:
        log(f"ZENODO_PUBLISH=false — draft {draft_id} left unpublished")
        log(f"  review: https://zenodo.org/deposit/{draft_id}")
        return

    log("publishing (mints the DOI)")
    published = api(
        "POST",
        f"{ZENODO_API}/deposit/depositions/{draft_id}/actions/publish",
    )
    doi = published.get("doi") or published.get("metadata", {}).get("doi")
    log(f"published: {doi}")
    log(f"  record: https://zenodo.org/records/{published['id']}")
    log(
        "  concept DOI still resolves to latest: "
        f"10.5281/zenodo.{CONCEPT_RECID}"
    )


if __name__ == "__main__":
    sys.exit(main())
