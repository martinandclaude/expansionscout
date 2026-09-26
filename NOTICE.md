# Third-party content

`expansionscout` itself is MIT licensed; see `LICENSE`. This file covers content
redistributed with it under different terms.

## STRchive locus catalogue

This repository redistributes a copy of the STRchive disease-locus catalogue
at `data/strchive/STRchive-loci.json`. It is the source of locus
coordinates, gene strand, repeat motifs in both orientations, documented
interruptions, and clinical size thresholds.

The file is compiled into the `expansionscout` binary, so every binary built
from this repository, including the release binaries, redistributes it too,
under the same terms.

- **Upstream**: <https://github.com/dashnowlab/STRchive> and <https://strchive.org>
- **Pinned commit and fetch date**: see `data/strchive/VERSION`
- **Data licence**: Creative Commons Attribution 4.0 International (CC BY 4.0),
  <https://creativecommons.org/licenses/by/4.0/>
- **Upstream code licence** (MIT), reproduced verbatim:
  `data/strchive/LICENSE.STRchive`

CC BY 4.0 requires attribution. If you use this tool, cite:

> Hiatt L, Weisburd B, Dolzhenko E, Rubinetti V, Avvaru AK, VanNoy GE,
> Kurtas NE, Rehm HL, Quinlan AR, Dashnow H. STRchive: a dynamic resource
> detailing population-level and locus-specific insights at tandem repeat
> disease loci. *Genome Medicine* 2025;17(1):29.
> doi:10.1186/s13073-025-01454-4

The catalogue file is redistributed unmodified. Locus interpretation applied
on top of it (methylation windows, which interruptions are scored, clinical
size-band labels) is this project's own and is not endorsed by the STRchive
authors.
