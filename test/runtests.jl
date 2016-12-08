# PLX.jl
# Copyright (C) 2012   Simon Kornblith

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

using PLX, MAT, BinDeps, Base.Test

# Unfortunately, these test files come from the SDK provided by Plexon Inc. and
# are assumed not to be redistributable
testdir = joinpath(dirname(@__FILE__), "Matlab Offline Files SDK", "mexPlex", "tests")
if !isdir(testdir)
	parentdir = dirname(@__FILE__)
	bundlezip = joinpath(parentdir, "Plexon Offline SDKs.zip")
	run(download_cmd("http://www.plexon.com/sites/default/files/downloads/OmniPlex%20and%20MAP%20Offline%20SDK%20Bundle_0.zip",
		             bundlezip))
	cd(parentdir) do
		run(`unzip $bundlezip`)
	end
	cd(parentdir) do
		run(`unzip $(joinpath(parentdir, "Plexon Offline SDKs", "Matlab Offline Files SDK.zip"))`)
	end
end

to_vec(x::Vector) = x
to_vec(x::Array) = vec(x)
to_vec(x) = [x]

m = matread(joinpath(testdir, "mexPlexData1.dat"))
for plx in m["data"]["plxs"]
	p = PLXFile(joinpath(testdir, plx["FileName"]), waveforms=true)
	@test p.header.Version == plx["Version"]
	@test p.header.ADFrequency == plx["Freq"]
	@test p.header.Comment == plx["Comment"]
	@test p.header.Trodalness == plx["Trodalness"]
	@test p.header.NumPointsWave == plx["NPW"]
	@test p.header.NumPointsPreThr == plx["PreThresh"]
	@test p.header.SpikeMaxMagnitudeMV == plx["SpikePeakV"]
	@test p.header.BitsPerSpikeSample == plx["SpikeADResBits"]
	@test p.header.SlowMaxMagnitudeMV == plx["SlowPeakV"]
	@test p.header.BitsPerSlowSample == plx["SlowADResBits"]
	@test p.header.LastTimestamp/p.header.ADFrequency == plx["Duration"]

	# Spikes
	n = plx["ts"]["n"]
	ts = plx["ts"]["ts"]
	wf = plx["wf"]["wf"]
	wf_v = plx["wf_v"]["wf"]
	for iunit = 1:30, ich = -1:130
		myts = ts[iunit+2, ich+2]
		if myts == -1
			@test !haskey(p.spike_channels, ich) ||
			      iunit > length(p.spike_channels[ich].units) ||
			      isempty(p.spike_channels[ich].units[iunit].spike_times)
		else
			ch = p.spike_channels[ich]
			unit = ch.units[iunit]
			@test vec(myts) == unit.spike_times
			@test wf[iunit+2, ich+2]' == unit.spike_waveforms
			@test_approx_eq wf_v[iunit+2, ich+2]' unit.spike_waveforms/3276.8
		end
	end

	# Continuous channels
	freq = plx["ad"]["freq"]
	ts = plx["ad"]["ts"]
	val = plx["ad"]["val"]
	val_v = plx["ad_v"]["val"]
	for ich = -1:400
		myts = ts[ich+2]
		if myts == -1
			@test !haskey(p.continuous_channels, ich) ||
			      isempty(p.continuous_channels[ich].samples)
		else
			ch = p.continuous_channels[ich]
			@test to_vec(myts) == ch.times.timestamps[1:end-1]/ch.times.timestamp_frequency
			@test frequency(ch.times) == freq[ich+2]
			@test ch.samples == vec(val[ich+2])
			@test ch.samples*ch.voltage_multiplier == vec(val_v[ich+2])
		end
	end

	# Events
	tsevs = plx["tsevs"]
	for ich = -1:310
		myts = tsevs[ich+2]
		if myts == -1
			@test !haskey(p.event_channels, ich) ||
			      isempty(p.event_channels[ich].times)
		else
			ch = p.event_channels[ich]
			@test_approx_eq ch.times to_vec(myts)
			if ich == 257
				v = to_vec(plx["strobed"]["v"])
				@test v == ch.codes
			end
		end
	end
end
