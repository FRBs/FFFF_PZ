# Fix mw_ebv

## Goals

Fix what might be a new API for the mw_ebv calculation for new FRBs

## Context

When we ingest a new FRB using the ingest_frbs API, it creates a new FRBTransient object.
At the end of its creation, it calls the execute_after_save() signal handler.
That executes a call to the IRSA dust API to get the mw_ebv value:

```python
c = SkyCoord(ra=instance.ra, dec=instance.dec, unit='deg')
instance.mw_ebv = irsa_dust.IrsaDust.get_query_table(c)['ext SandF mean'][0]
```

This call is currently failing.  We need to explore why and then 
implement a fix.

See the following code:

- YSE_App.models.frbtransient_models.py 
- YSE_App.models.frbtransient_models.execute_after_save()
- data_utils.ingest_frbs()

## Prompts

1.  Let's create a simple test case of the call to recover the mw_ebv value
for the IrsaDust API.  Create a method in a new module named testing/test_ingest.py
that does the following:

    - Creates a SkyCoord object
    - Calls the IrsaDust.get_query_table() method
    - Recovers the mw_ebv value from the result
    - Returns the mw_ebv value

    Then run the test and see if it passes.
    If you have any questions, put them in the Q&A/Testing section below.
    Use Opus 5.5.  Log your work

2.  I have put the error message from the server logs in the Q&A/Testing section below.
I believe the container is running astroquery 0.4.3.  Can you explore 
that code base to confirm/refute your hypothesis that it points to an old IRSA DUST endpoint that
IRSA has since retired/moved?  
If you have any questions, put them in the Q&A/Testing section below.
Use Opus 5.5.  Log your work

3. I am confident you have the right certificate and Python version.
Please provide the fixes for this, both short-term (hot-fix) and long-term (proper fix).
Note that we have not rebuilt the Docker image in awhile.
If you have any questions, put them in the Q&A/Solution section below.
Use Opus 5.5.  Log your work

## Solution

Root cause: the web image pins `certifi==2020.12.5`. IRSA's TLS cert (renewed
2026-11-12) chains to `SSL.com TLS RSA Root CA 2022`, which first appears in
certifi 2023.7.22. `requests` (and therefore astroquery) uses certifi's bundle,
so verification fails. Fix: get certifi >= 2023.7.22 into the container.
Current certifi (2026.11.20) requires Python >= 3.7, so it installs on the
image's Python 3.8. `requests==2.25.1` only needs `certifi>=2017.4.17`, so there
is no conflict.

### Short-term

Hot-fix on the running container (no rebuild, no code change; the app code is
volume-mounted at `/app`):

```bash
# 1. upgrade the CA bundle inside the running web container
docker exec ysepz_web_container pip install -U "certifi>=2023.7.22"

# 2. restart so the Django/gunicorn process picks up the new bundle
#    (runserver's autoreloader does not watch site-packages)
docker restart ysepz_web_container

# 3. verify (after pulling this branch on the server so testing/ exists)
docker exec ysepz_web_container python3 testing/test_ingest.py
#   expect: python 3.8.x, astroquery 0.4.3, certifi 2026.x ... / mw_ebv = 0.0155
```

Caveat: this change lives only in the container's writable layer. It survives
`docker restart` and `docker compose stop/start`, but it is **lost** on
`docker compose down`/`up`, `--force-recreate`, or a new image pull. Hence the
long-term fix.

Then re-ingest the FRBs that failed since ~2026-11-12. The handler deleted them
on failure, so re-running `ingest_frbs` for them is sufficient.

### Long-term

Bake the newer bundle into the derived image (this avoids rebuilding the old
2021 base image):

- **Done (uncommitted):** `docker/requirements.web.dev` now includes
  `certifi>=2023.7.22`. `Dockerfile.web.dev` (FROM
  `ghcr.io/davecoulter/yse_pz:latest`) `pip install`s that file, so pip will
  upgrade the base's 2020.12.5. A lower bound, not an exact pin, means every
  rebuild picks up the latest root certificates.
- Rebuild and redeploy (on the server, in `docker/`):
  ```bash
  docker compose -f dev-docker-compose.yml build web
  docker compose -f dev-docker-compose.yml up -d web     # recreates the container
  docker exec ysepz_web_container python3 testing/test_ingest.py
  # optionally: docker push ghcr.io/frbs/ffff_pz:latest
  ```
- Recommended follow-ups (not done):
  - Rebuild periodically (e.g. each semester). Any frozen CA bundle will break
    again the next time a remote host rotates to a new root.
  - Do **not** switch the import in `frbtransient_models.py` to
    `astroquery.ipac.irsa.irsa_dust` while the image has astroquery 0.4.3: that
    sub-package does not exist in 0.4.3 (verified in the sdist). Change it only
    together with an astroquery upgrade.
  - Optionally bump `certifi` in
    `docker/web_image_history/base/requirements/requirements_2.txt`, so that a
    future base rebuild does not reintroduce the old bundle.

## Q&A

### Solution

- Q: Is production started from `docker/dev-docker-compose.yml` (image
  `ghcr.io/frbs/ffff_pz:latest`, container `ysepz_web_container`)? I assumed so
  for the commands above. If you use `docker-compose.yml` instead, it pulls
  `ghcr.io/davecoulter/yse_pz:latest` directly and has **no** derived build
  step. In that case the long-term fix needs that compose file to switch to
  the `Dockerfile.web.dev` build.
- Q: Should a failed IRSA query still be fatal (delete the FRB and raise
  `FRBIngestionError`)? An alternative is to log a warning, leave
  `mw_ebv=None`, and continue, which would have kept ingestion running during
  this outage. The trade-off is that tagging runs after the mw_ebv calculation
  and may depend on it. I have left the current behaviour unchanged.

### Testing

- Q: The test passes locally (astroquery 0.4.10), but the docker requirements
  (`docker/web_image_history/base/requirements/requirements_{2,3}.txt`) pin
  `astroquery==0.4.3`. Which astroquery version is actually installed in the
  production web container? (e.g. `docker exec <web> python -c "import astroquery; print(astroquery.__version__)"`)
  My working hypothesis is that 0.4.3 points at an old IRSA DUST endpoint that
  IRSA has since retired/moved, which would explain the failure on the server
  but not locally.
- Q: Do you have the actual error message from the server logs
  (`Failed to query IRSA for mw_ebv: ...`)? That would confirm or rule out the
  hypothesis above.

>A. Here is the error message from the server logs:
```
Failed to query IRSA for mw_ebv: HTTPSConnectionPool(host='irsa.ipac.caltech.edu', port=443): Max retries exceeded with url: /cgi-bin/DUST/nph-dust (Caused by SSLError(SSLCertVerificationError(1, '[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: unable to get local issuer certificate (_ssl.c:1131)')))
```
- Q (prompt 2): Can you confirm the CA bundle in the running web container?
  `docker exec <web> python -c "import certifi, sys; print(sys.version, certifi.__version__, certifi.where())"`
  and `docker exec <web> env | grep -i -E "REQUESTS_CA_BUNDLE|SSL_CERT"`.
  I expect certifi 2020.12.5 (pinned in `requirements_2.txt`), Python 3.8,
  and no CA-bundle override.
- Q (prompt 2): For the fix, which do you prefer? (a) bump the certifi pin in the
  image (`certifi>=2023.7.22`, ideally latest) and rebuild; (b) a quick hot-fix
  `docker exec <web> pip install -U certifi` + restart, then do (a) properly;
  (c) also add `pip install -U certifi` to `docker/Dockerfile.web.dev` so the
  derived image refreshes it without rebuilding the base. My recommendation is
  (b) now, then (a)/(c). Upgrading astroquery is **not** needed for this bug.
- Note: 0.4.10 emits a DeprecationWarning: `irsa_dust` has moved to
  `astroquery.ipac.irsa.irsa_dust`. Not the cause of the failure, but worth
  updating the import in `frbtransient_models.py` as part of the fix.
  **Correction (prompt 3):** do not change the import yet. `astroquery.ipac`
  does not exist in 0.4.3, so changing it would break the container. Change it
  only together with an astroquery upgrade.

## Logs

### 2026-11-23 — Prompt 1 (Opus 5.5)

- Created `testing/test_ingest.py` with:
  - `get_mw_ebv(ra=150., dec=30.)` — builds a SkyCoord, calls
    `irsa_dust.IrsaDust.get_query_table(c)`, returns `['ext SandF mean'][0]`
    (same call as `execute_after_save()`).
  - `test_mw_ebv()` — pytest check that 0 <= mw_ebv < 1.
  - `__main__` block prints the value.
- Ran locally (conda env `astro`, astroquery 0.4.10):
  `cd testing && python -m pytest test_ingest.py` → **1 passed** (1.5 s).
  `python testing/test_ingest.py` → `mw_ebv = 0.0155`.
- Result: the IRSA dust API itself works with a current astroquery. The
  failure is most likely environment-specific to the server (see Q&A/Testing:
  pinned astroquery 0.4.3). Next step: run the same test inside the web
  container.

### 2026-11-23 — Prompt 2 (Opus 5.5)

**Verdict: hypothesis REFUTED.** The endpoint has not moved. The failure is an
outdated CA certificate bundle (certifi) in the container, triggered by IRSA
renewing its TLS certificate on 2026-11-12.

- astroquery 0.4.3 source (sdist from PyPI, unpacked in the scratchpad):
  `astroquery/irsa_dust/__init__.py` sets the server to
  `https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust`. That is **identical**
  to 0.4.10 (`astroquery/ipac/irsa/irsa_dust/__init__.py`). Diffing `core.py`
  shows only keyword-only-arg and docstring changes, with no change to the
  request or SSL behaviour. The server error also shows this same URL.
- The error itself is `SSLCertVerificationError ... unable to get local issuer
  certificate`. That is a trust-store problem, not a 404 or retired service.
- IRSA's live cert chain (`openssl s_client`):
  - leaf `irsaproxy1.ipac.caltech.edu`, **NotBefore Nov 12 2026**
  - ← `SSL.com TLS Issuing RSA CA R1` ← root **`SSL.com TLS RSA Root CA 2022`**
- Container image (`docker/web_image_history/base/`): Ubuntu focal, Python 3.8
  (matching `_ssl.c:1131` in the error), and `requirements_2.txt` pins
  `certifi==2020.12.5`, `requests==2.25.1`.
- `SSL.com TLS RSA Root CA 2022` in each certifi release's `cacert.pem`:
  2020.12.5 ✗, 2022.12.7 ✗, 2023.5.7 ✗, **2023.7.22 ✓** (first), 2023.11.17 ✓,
  2024.2.2 ✓.
- Reproduction: `requests.get(<nph-dust>, verify=<cacert.pem>)` with certifi
  2020.12.5's bundle reproduces the exact server error (`SSLError ...
  CERTIFICATE_VERIFY_FAILED`). With 2023.11.17's bundle it returns HTTP 200.
  This also explains why it worked before ~Nov 12 (the old cert presumably
  chained to a root that is in the 2020 bundle) and why it works locally
  (certifi 2025.04.26).
- Fix (not yet applied; see Q&A/Testing): upgrade certifi in the container to
  >= 2023.7.22. No astroquery change is required. The other `requests`-based
  HTTPS calls from the container to SSL.com-certified hosts will be fixed by
  the same upgrade.
### 2026-11-23 — Prompt 3 (Opus 5.5)

- Read the Docker setup. The base image (`web_image_history/base`, Ubuntu focal
  + Python 3.8) pins `certifi==2020.12.5`. `Dockerfile.web.dev` builds
  `ghcr.io/frbs/ffff_pz:latest` FROM `ghcr.io/davecoulter/yse_pz:latest` and
  `pip install`s `docker/requirements.web.dev`. No later `web_image_history/0.x`
  requirements file touches certifi or requests. The app code is
  volume-mounted (`${VOL}:/app`), so only the package change needs an image
  rebuild.
- Checked compatibility: certifi 2026.11.20 requires Python >= 3.7 (OK on 3.8).
  requests 2.25.1 requires certifi >= 2017.4.17 (OK).
- Checked the astroquery 0.4.3 sdist: there is no `astroquery/ipac/`. I
  corrected my earlier note so the import is not changed yet.
- **Edited** `docker/requirements.web.dev`: added `certifi>=2023.7.22` with a
  comment (long-term fix; takes effect on the next derived-image rebuild).
- **Edited** `testing/test_ingest.py`: the `__main__` block now also prints the
  python/astroquery/certifi versions, so
  `docker exec ysepz_web_container python3 testing/test_ingest.py` doubles as
  the hot-fix/rebuild verification. Re-ran locally: `python 3.13.2,
  astroquery 0.4.10, certifi 2025.04.26` / `mw_ebv = 0.0155`.
- Wrote the Short-term (hot-fix) and Long-term sections under `## Solution`.
  The hot-fix commands must be run on the server; I have not run them.
- Open questions are in Q&A/Solution (which compose file prod uses; whether an
  IRSA failure should stay fatal). Nothing committed.
