# PLX.jl

PLX.jl reads Plexon PLX files in Julia. It is several times faster than the MATLAB SDK provided by Plexon Inc.

## Quick Start

Because PLX.jl loads the entire contents of a given Plexon file into memory, you will want at least as much RAM as your largest Plexon file.

To read a Plexon file:

```julia
load("PLX")
using PLX
plx = PLXFile("/path/to/plexon/file.plx")
```

To read a Plexon file without LFPs:

```julia
load("PLX")
using PLX
plx = PLXFile("/path/to/plexon/file.plx", lfps=false)
```

To read a Plexon file including spike waveforms:

```julia
load("PLX")
using PLX
plx = PLXFile("/path/to/plexon/file.plx", waveforms=true)
```

To access spike times:

```julia
plx.spike_channels[n].units[m].spike_times
```

To access encodes and encode times:

```julia
plx.event_channels[257].times
plx.event_channels[257].codes
```

To find samples around given time points in a continuous channel, use:

```julia
channel = plx.continuous_channels[n]
channel.data[sample_index(channel.times, index_or_indices)]
```

For further documentation of the PLXFile type, read the source or use `idump(PLXFile)`.

## Implementation Notes

PLX.jl relies heavily on the functionality provided by `mmap_array` to read files. I'm not actually sure if this works on Windows. However, it provides a large (~2X) performance boost on Linux.

PLX.jl implements its own object (`SampleTimes`) to handle the timestamps on continuous channels, both to save memory and to optimize searching for sample indices corresponding to specific time points.

Documentation of the PLX format is available on Plexon's [website](http://www.plexon.com/downloads.html), or in PDF form [here](http://hardcarve.com/wikipic/PlexonDataFileStructureDocumentation.pdf).
