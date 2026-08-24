# Three parameters that stopped doing anything

`marginCutoffCoreNetwork`, `reducedNetworkRepresentation` and
`allowCysCoreSeeding` were accepted by the GUI, the CLI and `single_script.py`,
and assigned onto the pHinder instance — but **nothing in the package read them**.
The algorithms that once consumed them were rewritten for 7.0 and the parameters
were left behind.

They have been removed from the GUI. The CLI still accepts the flags so existing
scripts keep running, but warns and ignores them.

## What they used to do

Recovered from the archived sources under
`IsomLabPrivate/pythonScripts/pHinder/z_archives/` (2015–2017), so that
restoring any of them is a matter of decision, not rediscovery.

### `allowCysCoreSeeding`  (default 0)

`pHinder_2_4_triangulation`, ~line 3298. With the flag off **and**
`residueSet == "ionizableSet"`, cysteines were not allowed to nucleate cores
blindly — a CYS seeded a core only if it had one or more positively charged
neighbours, the signature of a pH-sensing cysteine. With the flag on, cysteines
seeded cores freely.

> `if not allowCysCoreSeeding and residueSet == "ionizableSet":`

The 7.0 core-seeding path contains no cysteine special case at all.

### `reducedNetworkRepresentation`  (default 1)

`minimizeNetworks`, ~line 363 of the same file. It chose the distance threshold
for keeping loop edges:

- `0` — keep every edge ≤ **6.0 Å**, described in the source as an arbitrary
  choice that "makes the resultant networks a bit easier to visualize"
- `1` — keep only edges ≤ **4.0 Å**, the author's "operational definition for
  salt bridges and anti-pairs"

The current `geometry/minimizeNetworks.py` is a complete rewrite with no
loop-edge stage, so there is no threshold left to switch.

### `marginCutoffCoreNetwork`  (default −2.0)

`pHinder_2_5_triangulation`, ~line 3582. Margin sidechains are classified in
`classifySidechains()`; this value decided which of them were then admitted into
**core networks**, by depth:

| Value  | Behaviour          |
|--------|--------------------|
| −100   | no margin admitted |
| −2.0   | restricted margin (the documented default) |
| 0.0    | liberal margin     |

> `if depth < marginCutoffCoreNetwork: networkBigs.update(...)`

Admitting margin nodes let separate sub-networks connect through margin as well
as core. 7.0 keeps `core_networks` attributes but not this admission step.

## Restoring one

Each needs its algorithm re-implemented against the 7.0 data structures, not
just a variable re-wired — and each changes scientific output, so it is a
deliberate choice rather than a cleanup. The archived implementations above are
the reference.
