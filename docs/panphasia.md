# PANPHASIA descriptors

```ini
[random]
generator  = PANPHASIA
descriptor = [Panph6,L18,(56034,71400,250000),S1,KK1025,CH3774755196,Flamingo_Gpc1]
```

- `generator = PANPHASIA`: use [PANPHASIA](https://arxiv.org/abs/1306.5968)
  white noise. Needs `ENABLE_PANPHASIA`.
- `descriptor`: the noise to use. Both original `Panph1` and high-order
  `Panph6` descriptors (such as FLAMINGO's) work; spaces are ignored. Copy it
  exactly: the `CH` check digit must match.

## `Panph6`

`Panph6` noise is monofonIC's field sampled at cell centres, where MUSIC
places particles (monofonIC uses cell corners), on every level. In zooms each
refined patch is computed on its padded region, and its long-wavelength modes
are replaced by the coarser level's.

Requirements:

- the base size `S` must be a power of two;
- with a `KK` term, which fixes the amplitudes of all modes with |k|² ≤ KK
  (box units), `levelmin` must satisfy 2^(levelmin-1) > √KK (e.g. 7 for
  KK1025), and `fix_mode_amplitude` must not be set;
- the noise needs about 200 bytes per cell of the largest padded level.

The high-order routines are a serial copy of monofonIC's, in
`ext/panphasia_ho` (see its README for provenance).
