# PLX.jl
# Methods for reading Plexon PLX files in Julia

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

module PLX
using Compat

import Base.length, Base.read, Base.searchsortedlast, Base.searchsortedfirst

export PL_FileHeader, PL_ChanHeader, PL_EventHeader, PL_SlowChannelHeader, SampleTimes, PLXUnit,
	PLXSpikeChannel, PLXEventChannel, PLXContinuousChannel, PLXFile, sample_index, dt, frequency

function struct_string(bytes::Vector{UInt8})
    last_byte = findfirst(bytes, 0)
    String(last_byte == 0 ? bytes : bytes[1:last_byte-1])
end

macro struct(typename, contents)
	if contents.head != :block
		error("Invalid struct declaration")
	end

	blk = Expr(:call, typename)

	for typedecl in contents.args
		if typedecl.head == :line
			continue
		elseif typedecl.head != :(::)
			error("Invalid struct declaration")
		end
		fieldname = typedecl.args[1]
		fieldtype = typedecl.args[2]
		if isa(fieldtype, Expr)
			if fieldtype.head != :call
				error("Invalid struct declaration")
			end
			# Type has size parameters

			if fieldtype.args[1] in (:String, :AbstractString)
				typedecl.args[2] = fieldtype.args[1]
				push!(blk.args, :(struct_string(read(ios, UInt8, ($(fieldtype.args[2:end]...))))))
			else
				typedecl.args[2] = :(Array{$(fieldtype.args[1]), $(length(fieldtype.args)-1)})
				push!(blk.args, :(read(ios, $(fieldtype.args[1]), ($(transpose(fieldtype.args[2:end])...)))))
			end
		else
			push!(blk.args, :(read(ios, $fieldtype)))
		end
	end
	
	quote
		global read
		type $typename
			$contents
		end
		read(ios::IOStream, ::Type{$typename}) = $blk
	end
end

@struct PL_FileHeader begin
	MagicNumber::UInt32
	Version::Int32
	Comment::String(128)
	ADFrequency::Int32
	NumDSPChannels::Int32
	NumEventChannels::Int32
	NumSlowChannels::Int32
	NumPointsWave::Int32
	NumPointsPreThr::Int32

	Year::Int32
	Month::Int32
	Day::Int32
	Hour::Int32
	Minute::Int32
	Second::Int32

	FastRead::Int32
	WaveformFreq::Int32
	LastTimestamp::Float64

	Trodalness::UInt8
	DataTrodalness::UInt8
	BitsPerSpikeSample::UInt8
	BitsPerSlowSample::UInt8
	SpikeMaxMagnitudeMV::UInt16
	SlowMaxMagnitudeMV::UInt16
	SpikePreAmpGain::UInt16

	Padding::UInt8(46)

	TSCounts::Int32(5, 130)
	WFCounts::Int32(5, 130)

	EVCounts::Int32(512)
end

@struct PL_ChanHeader begin
	Name::String(32)
	SIGName::String(32)
	Channel::Int32
	WFRate::Int32
	SIG::Int32
	Ref::Int32
	Gain::Int32
	Filter::Int32
	Threshold::Int32
	Method::Int32
	NUnits::Int32
	Template::Int16(64, 5)
	Fit::Int32(5)
	SortWidth::Int32
	Boxes::Int16(4, 2, 5)
	SortBeg::Int32
	Comment::String(128)
	Padding::Int32(11)
end

@struct PL_EventHeader begin
	Name::String(32)
	Channel::Int32
	Comment::String(128)
	Padding::Int32(33)
end

@struct PL_SlowChannelHeader begin
	Name::String(32)
	Channel::Int32
	ADFreq::Int32
	Gain::Int32
	Enabled::Int32
	PreAmpGain::Int32

	SpikeChannel::Int32

	Comment::String(128)
	Padding::Int32(28)
end

abstract SpikeChannel
abstract SpikeFile

type SampleTimes
	# Timestamps as integers
	timestamps::Vector{Int}
	# Indices of samples corresponding to each timestamp
	timestamp_indices::Vector{Int}
	# timestamps/timestamp_frequency gives times in seconds
	timestamp_frequency::UInt64
	# Distance between samples in units of 1/timestamp_frequency
	# Sample frequency is timestamp_frequency/sample_dt
	sample_dt::UInt64

	SampleTimes(n::Int, timestamp_frequency::Integer, sample_dt::Integer) =
		new(Array(Int, n), Array(Int, n), timestamp_frequency, sample_dt)
	SampleTimes(timestamps::Vector{Int}, timestamp_indices::Vector{Int},
		timestamp_frequency::Integer, sample_dt::Integer) =
		optimize(new(timestamps, timestamp_indices, timestamp_frequency, sample_dt))
end

type PLXUnit
	id::UInt64
	spike_times::Vector{Float64}
	spike_waveforms::Union{Array{Int16, 2}, Void}
	voltage_multiplier::Float64

	PLXUnit(unit_number::Int, channel::SpikeChannel, n_spikes::Int, waveforms::Bool=false) =
		new(hash((channel.id, unit_number)), Array(Float64, n_spikes),
			waveforms ? Array(Int16, (convert(Int, channel.points_per_waveform), n_spikes)) : nothing,
			channel.voltage_multiplier)
end

type PLXSpikeChannel <: SpikeChannel
	id::UInt64
	header::PL_ChanHeader
	units::Vector{PLXUnit}
	unclustered::PLXUnit
	voltage_multiplier::Float64
	points_per_waveform::Int32

	function PLXSpikeChannel(plxfile::SpikeFile, header::PL_ChanHeader)
		x = new(hash((plxfile.id, header.Channel)), header)
		x.units = Array(PLXUnit, header.NUnits)
		x.voltage_multiplier = if plxfile.header.Version >= 105
				plxfile.header.SpikeMaxMagnitudeMV/
					((1 << (plxfile.header.BitsPerSpikeSample-1))*header.Gain*plxfile.header.SpikePreAmpGain*1000)
			elseif plxfile.header.Version >= 103
				plxfile.header.SpikeMaxMagnitudeMV/
					((1 << (plxfile.header.BitsPerSpikeSample-1))*header.Gain*1000)
			else
				3000/(2048*header.Gain*1000)
			end
		x.points_per_waveform = plxfile.header.NumPointsWave
		return x
	end
end

type PLXEventChannel
	header::PL_EventHeader
	times::Vector{Float64}
	codes::Union{Vector{Int16}, Void}

	PLXEventChannel(plxfile::SpikeFile, header::PL_EventHeader) = new(header)
end

type PLXContinuousChannel
	id::UInt64
	header::PL_SlowChannelHeader
	samples::Vector{Int16}
	times::SampleTimes
	voltage_multiplier::Float64

	function PLXContinuousChannel(plxfile::SpikeFile, header::PL_SlowChannelHeader)
		x = new(hash((plxfile.id, header.Channel+1)), header)

		x.voltage_multiplier = if plxfile.header.Version >= 103
				plxfile.header.SlowMaxMagnitudeMV/
					((1 << (plxfile.header.BitsPerSpikeSample-1))*header.Gain*header.PreAmpGain)
			elseif plxfile.header.Version == 102
				5000/(2048*header.Gain*header.PreAmpGain)
			else
				5000/(2048*1000*header.Gain)
			end

		return x
	end
end

type PLXFile <: SpikeFile
	id::UInt64
	header::PL_FileHeader
	spike_channels::Dict{Int, PLXSpikeChannel}
	event_channels::Dict{Int, PLXEventChannel}
	continuous_channels::Dict{Int, PLXContinuousChannel}

	function PLXFile(file_name::AbstractString; lfps::Bool=true, waveforms::Bool=false)
		const maxChannels = 32*1024

		x = new()
		ios = open(file_name, "r")

		x.header = read(ios, PL_FileHeader)
		if x.header.MagicNumber != 0x58454c50
			error("$file_name does not appear to be a PLX file")
		elseif x.header.NumDSPChannels < 0 || x.header.NumDSPChannels > maxChannels ||
				x.header.NumEventChannels < 0 || x.header.NumEventChannels > maxChannels ||
				x.header.NumSlowChannels < 0 || x.header.NumSlowChannels > maxChannels
			error("PLX file header specifies an invalid number of channels")
		elseif x.header.ADFrequency <= 0
			error("PLX file header specifies ADFrequency <= 0")
		elseif x.header.NumPointsWave <= 0
			error("PLX file header specifies no points in waveform")
		end

		x.id = hash((x.header.Year, x.header.Month, x.header.Day, x.header.Hour, x.header.Minute, x.header.Second))

		max_unit = 0
		x.spike_channels = Dict{Int, PLXSpikeChannel}()
		sizehint!(x.spike_channels, x.header.NumDSPChannels)
		for i=1:x.header.NumDSPChannels
			header = read(ios, PL_ChanHeader)
			ch = convert(Int, header.Channel)
			x.spike_channels[ch] = PLXSpikeChannel(x, header)
			if header.NUnits > max_unit
				max_unit = header.NUnits
			end
		end

		x.event_channels = Dict{Int, PLXEventChannel}()
		sizehint!(x.spike_channels, x.header.NumEventChannels)
		for i=1:x.header.NumEventChannels
			header = read(ios, PL_EventHeader)
			ch = convert(Int, header.Channel)
			x.event_channels[ch] = PLXEventChannel(x, header)
		end

		x.continuous_channels = Dict{Int, PLXContinuousChannel}()
		sizehint!(x.spike_channels, x.header.NumSlowChannels)
		for i=1:x.header.NumSlowChannels
			header = read(ios, PL_SlowChannelHeader)
			ch = convert(Int, header.Channel)
			x.continuous_channels[ch] = PLXContinuousChannel(x, header)
		end

		data_offset = position(ios)

		# Read through the file once to determine how much memory to allocate
		n_spikes = zeros(Int, (convert(Int, maximum(keys(x.spike_channels))), max_unit+1))
		n_events = zeros(Int, maximum(keys(x.event_channels)))
		n_timestamps = zeros(Int, maximum(keys(x.continuous_channels))+1)
		n_samples = zeros(Int, size(n_timestamps))

		n_blocks = 0
		seekend(ios)
		fsize = position(ios)
		max_offset = div(fsize-data_offset, 2)
		contents = Mmap.mmap(ios, Matrix{Int16}, (1, max_offset), data_offset)
		cur_offset = 1
		while cur_offset < max_offset
			block_type = contents[cur_offset]
			ch = contents[cur_offset+4]

			if block_type == 1			# spike
				n_spikes[ch, contents[cur_offset+5]+1] += 1
			elseif block_type == 4		# event
				n_events[ch] += 1
			elseif block_type == 5		# continuous
				n_timestamps[ch+1] += 1
				n_samples[ch+1] += contents[cur_offset+7]
			else
				error(strcat("Invalid data block type ", t))
			end

			cur_offset += contents[cur_offset+7]+8
			n_blocks += 1
		end

		# Allocate
		for (i, channel)=x.spike_channels
			channel.unclustered = PLXUnit(0, channel, n_spikes[i, 1], waveforms)
			for j=1:channel.header.NUnits
				channel.units[j] = PLXUnit(j, channel, n_spikes[i, j+1], waveforms)
			end
		end
		for (i, channel)=x.event_channels
			channel.times = Array(Float64, n_events[i])
			if i == 257
				channel.codes = Array(Int16, n_events[i])
			else
				channel.codes = nothing
			end
		end
		for (i, channel)=x.continuous_channels
			sample_dt = x.header.ADFrequency/channel.header.ADFreq
			if sample_dt % 1 != 0
				error("Channel $i frequency $(channel.header.ADFreq) is non-integer multiple of AD frequency $(x.header.ADFrequency)")
			end
			channel.samples = Array(Int16, lfps ? n_samples[i+1] : 0)
			channel.times = SampleTimes(lfps ? n_timestamps[i+1] : 0, x.header.ADFrequency, @compat Int64(sample_dt))
		end

		# Read through the file again
		cur_spike = zeros(Int, size(n_spikes))
		cur_event = zeros(Int, size(n_events))
		cur_timestamp = zeros(Int, size(n_timestamps))
		cur_sample = zeros(Int, size(n_samples))
		cur_offset = 1
		while cur_offset < max_offset
			block_type = contents[cur_offset]
			timestamp = (convert(Int64, reinterpret(UInt16, contents[cur_offset+1])) << 32) +
				reinterpret(UInt16, contents[cur_offset+2]) +
				(convert(Int64, reinterpret(UInt16, contents[cur_offset+3])) << 16)
			ch = convert(Int, contents[cur_offset+4])

			if block_type == 1			# spike
				unit_number = contents[cur_offset+5]
				if unit_number == 0
					unit = x.spike_channels[ch].unclustered
				else
					unit = x.spike_channels[ch].units[unit_number]
				end
				c = (cur_spike[ch, unit_number+1] += 1)

				unit.spike_times[c] = timestamp/x.header.ADFrequency

				block_waveforms = contents[cur_offset+7]
				if waveforms
					# Unrolling doesn't seem to help here
					unit.spike_waveforms[:, c] = block_waveforms == 0 ? -32768 : contents[cur_offset+8:cur_offset+8+block_waveforms-1]
				end
				cur_offset += block_waveforms
			elseif block_type == 4		# event
				channel = x.event_channels[ch]
				c = (cur_event[ch] += 1)
				
				channel.times[c] = timestamp/x.header.ADFrequency
				if ch == 257
					channel.codes[c] = contents[cur_offset+5]
				end
			elseif block_type == 5		# continuous
				block_samples = contents[cur_offset+7]
				if block_samples > 0
					t = (cur_timestamp[ch+1] += 1)
					if lfps
						channel = x.continuous_channels[ch]
						c = cur_sample[ch+1]

						channel.times.timestamps[t] = timestamp
						channel.times.timestamp_indices[t] = c+1

						# Unrolling seems to help here
						start_offset = cur_offset + 7
						samples = channel.samples
						for i = 1:block_samples
							samples[c+i] = contents[start_offset+i]
						end
					end
					cur_sample[ch+1] += block_samples
					cur_offset += block_samples
				end
			end
			cur_offset += 8
		end

		@assert cur_spike == n_spikes
		@assert cur_event == n_events
		@assert cur_timestamp == n_timestamps
		@assert cur_sample == n_samples

		if lfps
			for (i, channel)=x.continuous_channels
				times = channel.times
				if !isempty(times.timestamp_indices)
					times.timestamp_indices[end] = cur_sample[i+1]
					times.timestamps[end] = times.timestamps[end-1]+
						times.sample_dt*(cur_sample[i+1]-times.timestamp_indices[end-1])
					optimize_times(times)
				end
			end
		end

		close(ios)

		return x
	end
end

Base.:(==)(A::SampleTimes, B::SampleTimes) =
	A.timestamps == B.timestamps && A.timestamp_indices == B.timestamp_indices &&
	A.timestamp_frequency == B.timestamp_frequency && A.sample_dt == B.sample_dt

# Find the timestamp at a given point (not dividing by timestamp_frequency)
function _int_ref(x::SampleTimes, y::Int)
	index = searchsortedlast(x.timestamp_indices, y)
	if index == 0
		throw(BoundsError())
	end
	return x.timestamps[index]+(y-x.timestamp_indices[index])*x.sample_dt
end

# Find the timestamps in a given range (not dividing by timestamp_frequency)
function _int_ref(x::SampleTimes, y::Range{Int})
	n = length(y)
	m = length(x.timestamps)

	out = similar(x.timestamps, n)
	index = _binary_search_leq(x.timestamp_indices, y[1])
	if index == 0
		throw(BoundsError())
	end

	y_start = 1
	next_timestamp_index = x.timestamp_indices[index]

	while index <= m
		timestamp_index = next_timestamp_index
		if index == m
			next_timestamp_index = Inf
		else
			next_timestamp_index = x.timestamp_indices[index+1]
		end
		if timestamp_index <= y[y_start]
			d = next_timestamp_index-y[y_start]
			if isa(y, Range1)
				y_end = y_start+d-1
			else
				y_end = y_start+div(d, y.step)
				if d % y.step == 0
					y_end -= 1
				end
			end

			if y_end >= n
				out[y_start:n] = x.timestamps[index]+((y[y_start:end])-timestamp_index)*x.sample_dt
				return out
			else
				out[y_start:y_end] = x.timestamps[index]+((y[y_start:y_end])-timestamp_index)*x.sample_dt
				y_start = y_end+1
			end
		end
		index += 1
	end
	throw(BoundsError())
end

# Get the sample index corresponding to a given timepoint
function _sample_index(x::SampleTimes, timepoint::Real)
	timestamp = float64(timepoint)*x.timestamp_frequency
	index = findfirst(y -> y >= timestamp, x.timestamps)

	if index == 0
		return 0
	elseif x.timestamps[index] == timepoint
		return index
	elseif index == 1
		return 0
	end

	sample = x.timestamp_indices[index-1] + (timestamp - x.timestamps[index-1])/x.sample_dt
	if sample > x.timestamp_indices[index]
		return 0
	end
	return sample
end

# Optimize SampleTimes by eliminating redundant timestamps
function optimize_times(x::SampleTimes)
	if isempty(x.timestamps)
		return x
	end

	ts_diff = diff(x.timestamps)
	index_diff = diff(x.timestamp_indices)
	keep = find(ts_diff .!= index_diff*x.sample_dt).+1
	last = length(x.timestamps)
	if isempty(keep) || keep[end] != last
		keep = [1; keep; last]
	else
		keep = [1; keep]
	end
	x.timestamps = x.timestamps[keep]
	x.timestamp_indices = x.timestamp_indices[keep]
	return x
end

# Allow indexing with integer indices (yields times) and float indices (yields integer indices)
Base.getindex(x::SampleTimes, y::Union{Range{Int}, Int}) = _int_ref(x, y)/x.timestamp_frequency
Base.getindex(x::SampleTimes, y::AbstractFloat) = sample_index(x, y)

length(x::SampleTimes) = x.timestamp_indices[end]

start(x::SampleTimes) = 1
next(x::SampleTimes, state) = (x[state], state+1)
done(x::SampleTimes, state) = state > length(x)

dt(x::SampleTimes) = x.sample_dt/x.timestamp_frequency
frequency(x::SampleTimes) = x.timestamp_frequency/x.sample_dt

# Integer sample indices
searchsortedfirst(x::SampleTimes, y::Real) = ifloor(_sample_index(x, y))
searchsortedlast(x::SampleTimes, y::Real) = iceil(_sample_index(x, y))
sample_index(x::SampleTimes, y::Real) = iround(_sample_index(x, y))
sample_index(x::SampleTimes, y1::Real, y2::Real) = iround(_sample_index(x, y1)):iround(_sample_index(x, y2))
sample_index{T <: Real}(x::SampleTimes, y::Vector{T}) = [sample_index(x, z) for z in y]
end
