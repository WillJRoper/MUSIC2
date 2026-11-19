# PANPHASIA descriptors

Select [PANPHASIA](https://arxiv.org/abs/1306.5968) white noise with

```ini
[random]
generator  = PANPHASIA
descriptor = [Panph6,L18,(56034,71400,250000),S1,KK1025,CH3774755196,Flamingo_Gpc1]
```

Both original `Panph1` descriptors and high-order `Panph6` descriptors (such
as FLAMINGO's) are supported; spaces inside the descriptor are ignored.

## `Panph6`

`Panph6` white noise is monofonIC's field sampled at cell centres, where MUSIC
places its particles (monofonIC places them at cell corners), on every
refinement level. In zooms each refinement patch is computed on its padded
region, and its long-wavelength modes are replaced by those of the coarser
level.

Requirements:

- the base size `S` must be a power of two;
- a `KK` limit fixes the amplitudes of all modes with |k|² ≤ KK (in box
  units). This is done on the periodic top level and carried into every
  refinement level, so `levelmin` must satisfy 2^(levelmin-1) > √KK, and
  `fix_mode_amplitude` cannot be used as well;
- generating the noise needs roughly 200 bytes per cell of the largest padded
  refinement level.

The high-order routines are a serial copy of monofonIC's, in
`ext/panphasia_ho` (see its README for provenance).
