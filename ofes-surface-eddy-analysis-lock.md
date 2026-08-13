# OFES surface-eddy analysis lock

Status: locked before any 56-event PET association input was read.
Lock date: 2026-08-14 (Asia/Shanghai)
Scope: OFES NP30 2003 daily model-SSH eddies detected by PET.

## 1. Environment and implementation identity

- Worktree: `/mnt/w2/scratch/user3/Oxygen-ofes-surface`
- Branch: `feat/ofes-surface-eddies`
- PET environment: `ofes-pet`
- Python: 3.10.20
- `pyeddy_tracker`: 3.6.1
- NumPy/SciPy/xarray/netCDF4/pandas/pyarrow/Matplotlib/PyYAML:
  1.21.6 / 1.11.2 / 2023.7.0 / 1.5.8 / 2.0.3 / 17.0.0 /
  3.8.4 / 6.0.3
- PET module path:
  `/home/user3/.conda/envs/ofes-pet/lib/python3.10/site-packages/py_eddy_tracker/__init__.py`
- PET module SHA-256:
  `157c5f8b2b2c40477f2b94aefaf2242e2cab8910985a2b880608ed9bd8fe123c`
- Audited PET source SHA-256:
  `dataset/grid.py` = `2f6f22dbb8200c4e337852de4fb8ce974d09267995119a9a5014b684fde39c7f`
  `observations/observation.py` = `134be8a9d82165acf143508b386e284f0bea832f4a52b1925ac729cd21aa434b`
  `poly.py` = `a667a404c322b4820592fb82228aff802f82ed095779abb943399a239708e1f2`
  `tracking.py` = `5fbcb6b7714a5b9a2e8666a00c04588c5ce8abbc906118a5c394370b4583a95c`
- `h5py` is absent from the PET environment. The existing `track.py` is not
  imported by the detector or runner.

## 2. Input and scientific scope

- Input files: `./data/OFES_NP30/eta/eta.MM.DD.2003.nc`.
- Date window: 2003-01-01 through 2003-12-31, exactly 365 dates.
- Delivered domain: approximately 25--45°N, 140--170°E; no global wrapping
  or outside-domain interpolation.
- ETA is read as `(time, lat, lon)` and converted from raw cm to m by `0.01`.
  Longitude and latitude must be one-dimensional, strictly increasing, and the
  data are transposed to PET's `(lon, lat)` order. The source land mask is
  retained as a NumPy masked array.
- No DO, temperature, salinity, deep-event label, deep Ro, McCoy result, or
  expected hit rate is read by the identity detector. OFES native u/v are not
  identity inputs.
- No daily or annual spatial mean is removed. The input remains model SSH,
  not SLA.

## 3. High-pass and geostrophic speed

The available worktree/HPC search did not contain a usable META3.1 code tag
`v3.3.1` production filter call. Therefore this is explicitly a
**META-family PET implementation**, not a byte-for-byte META3.1/META3.2
reproduction:

```python
grid.bessel_high_filter(
    "eta", 700.0, order=1, lat_max=85.0, extend=False
)
```

The filter is the PET 3.6.1 Bessel kernel. It uses masked ocean values and
does not fill the clipped domain with outside data. The PET official order is
then applied to the filtered field:

```python
grid.add_uv("eta", "u", "v", stencil_halfwidth=4)
```

The resulting u/v are used only for PET speed-contour construction; detection
identity remains the 700-km high-pass ETA contour.

`filter_valid_mask` is defined per latitude using the actual PET
`kernel_bessel(lat, 700, order=1)` footprint: the complete kernel footprint
must lie inside the delivered lon/lat array and the PET high-pass result must
not be masked. A clipped-domain point is not filter-valid. The PET ocean/land
mask remains a separate validity condition. A contour touching the delivered
outer boundary or any filter-invalid point is marked `boundary_censored`.
An ocean contour near a real coastline is not discarded solely for being
coastal; its `coastal` flag is recorded using PET's land-mask semantics.

## 4. Detection parameters

The fixed PET call is:

```python
grid.eddy_identification(
    "eta", "u", "v", date,
    step=0.002,
    shape_error=70.0,
    pixel_limit=(5, 2000),
    presampling_multiplier=10,
    sampling=50,
    sampling_method="visvalingam",
    mle=1,
    nb_step_min=2,
    nb_step_to_be_mle=2,
    force_height_unit="m",
    force_speed_unit="m/s",
)
```

PET 3.6.1 source documents `amplitude_threshold = step * nb_step_min`,
so `0.002 * 2 = 0.004 m`. `mle=1` allows one accepted SSH extremum;
`nb_step_to_be_mle=2` is PET's secondary-extremum separation threshold. The
“three closed contour” wording is represented by the PET level sequence of
the outer contour plus two 0.002-m intervals to the extremum; no extra
post-hoc contour detector is introduced. The source/API meaning, rather than
an object-count calibration, is the frozen rule.

The output retains PET effective contour, speed contour, effective/speed
radius, amplitude, polarity, SSH extremum, center, and mask/boundary flags.

## 5. Tracking parameters

- Cyclonic and anticyclonic objects are tracked separately.
- Links use effective-contour intersection-over-union, threshold 0.05.
- For multiple candidates, the largest overlap is selected; ties are resolved
  deterministically by object id.
- Candidate observations may be searched 1--5 days ahead.
- At most four consecutive virtual observations are inserted. Virtual rows
  are explicitly flagged and are never counted as actual daily PET detections.
- No minimum track duration is imposed. In-window lifetime categories are:
  `untracked` (1 day), `short` (2--9), `long` (>=10), and `persistent`
  (>=30, sensitivity label).
- Year-start/year-end, delivered-domain, and filter-boundary censor flags are
  retained. A one-year track lifetime is a lower bound, not a complete-life
  claim.

## 6. Validation and execution gates

Before formal output, the runner must pass the official PET demo, source/API
checks, unit/transpose/mask tests, one fixed non-event OFES day, and four
pre-registered engineering dates: 2003-01-15, 2003-04-15, 2003-07-15,
2003-10-15. Object counts cannot be used to tune the locked parameters.

Only a manifest with all 365 daily statuses `complete` may enter event
association. The runner uses explicit `run_id`/manifest paths and never
selects a run by `sorted(glob)[-1]`. Scientific signatures include code,
protocol, config, PET source, and ETA inventory hashes; worker count is not in
the signature.

The authoritative DO50 population/catalog/diagnostic inputs will be supplied
explicitly to the postprocessor and hashed at that later gate. The 56-event
association, nulls, rotation cross-tab, and McCoy cross-tab are intentionally
not part of this lock commit and have not been read.
