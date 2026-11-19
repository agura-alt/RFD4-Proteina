"""Keep the seed hotspots that SampleNearbyHotspots would otherwise replace.

`rfproteina.transforms.hotspot_sampling.SampleNearbyHotspots` REPLACES the C_HOT
annotation: it centres a patch on the seed atoms' centroid, excludes the seed
residues from the candidate pool (hotspot_sampling.py:87), and writes a fresh mask
over the result. The seeds therefore steer where the patch sits and are then
dropped.

For a campaign whose epitope was chosen deliberately -- CSS picked A/204 LEU CD2,
A/206 LEU CD2 and A/246 TYR CE1 for solvent exposure, after the originally
requested CD1/CD1/CZ atoms measured 0.00/3.52/0.50 A^2 and would have made the
hotspot gate unreachable -- that is the wrong behaviour. This wrapper runs the
sampler and then ORs the seeds back in, so every design is conditioned on the
chosen atoms PLUS a varying draw from around them.

Consequences worth stating, because they change what the numbers mean:

  * The conditioned hotspot count is seeds + draw, not draw. With three seeds and
    an inner `num_hotspots: [2, 3]` that is 5 or 6 atoms per design.
  * `hotspot_obedience` is a FRACTION of conditioned hotspots contacted, so a gate
    of `min: 1.0` now demands all 5 or 6 within max_hotspot_distance (4.5 A), not
    all 3. That is a materially harder gate than the fixed arm's, in the opposite
    direction from the seed-dropping arm's easier 2-of-2 draws. Neither pass rate
    is comparable with the other; record the per-design count alongside the gate.
  * The inner sampler excludes the seed RESIDUES, so the draw can never duplicate
    a seed atom. The counts add exactly.

This module lives under `src/rfproteina/transforms/sampling/`, which
configs/inference_experiment/TEMPLATE.yaml designates for user sampling transforms.
That path is git-ignored (.gitignore:9), so this file was force-added to keep it with
the CSS config that names it. The directory stays ignored: a new module beside this one
will not be tracked unless it is force-added too.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from atomworks.ml.transforms.base import Transform

from rfproteina.conditions import C_HOT


class KeepSeedHotspots(Transform):
    """Run an inner hotspot sampler, then restore the seed hotspots it dropped.

    Args:
        sampler: The hotspot transform to wrap, instantiated by Hydra from a nested
            `_target_`. It must leave atom order and identity unchanged; anything
            that reindexes the structure is rejected rather than silently remapped.
        min_sampled: Minimum number of NEW hotspots the sampler must contribute on
            top of the seeds. 1 by default, so a sampler that quietly returned the
            seeds unchanged -- the failure this wrapper is most likely to mask --
            raises instead of looking like a successful draw.
    """

    def __init__(self, sampler: Transform, min_sampled: int = 1) -> None:
        if not isinstance(sampler, Transform):
            raise TypeError(f"sampler must be a Transform, got {type(sampler).__name__}.")
        if not isinstance(min_sampled, int) or min_sampled < 0:
            raise ValueError("min_sampled must be a non-negative integer.")
        self.sampler = sampler
        self.min_sampled = min_sampled

    @staticmethod
    def _identity(array) -> tuple:
        return (
            array.array_length(),
            tuple(map(str, array.chain_id)),
            tuple(map(int, array.res_id)),
            tuple(map(str, array.atom_name)),
        )

    def forward(self, data: dict[str, Any]) -> dict[str, Any]:
        array = data["atom_array"]
        seeds = C_HOT.mask(array, default="generate")
        if not seeds.any():
            raise ValueError(
                f"{data.get('example_id', 'Input')}: no C_HOT seed to keep. "
                "KeepSeedHotspots is only meaningful on inputs that carry seed hotspots."
            )
        before = self._identity(array)

        result = self.sampler.forward(data)
        sampled_array = result["atom_array"]
        if self._identity(sampled_array) != before:
            raise ValueError(
                f"{data.get('example_id', 'Input')}: {type(self.sampler).__name__} changed atom "
                "identity or order, so the seed mask cannot be restored positionally."
            )

        output = sampled_array.copy() if sampled_array is array else sampled_array
        sampled = C_HOT.mask(output, default="generate")
        new = sampled & ~seeds
        if int(new.sum()) < self.min_sampled:
            raise ValueError(
                f"{data.get('example_id', 'Input')}: {type(self.sampler).__name__} contributed "
                f"{int(new.sum())} hotspot(s) outside the seeds, below min_sampled={self.min_sampled}."
            )
        C_HOT.set_annotation(output, np.asarray(seeds | sampled, dtype=bool))
        return {**result, "atom_array": output}
