using MRIFieldmaps: b0map, b0model, b0init, b0scale

"""
CalculateSensitivitityMap(RawData::RawAcquisitionData,params::Dict{Symbol,Any})
MtxSize = Tuple (nx,ny,nz) = Final size of Map
It expects raw data to be 3D
Calculates Sensitivity map
"""
function CalculateSensitivityMap(recon,MtxSize::Tuple;calib_size::Int=12)

    SensitivityMap = Array{ComplexF32}(undef,MtxSize)
    calibration = Array{ComplexF32}(undef,size(recon)[1:4])

    numChan = size(recon)[4]

    for i_ch = 1:numChan
        calibration[:,:,:,i_ch] = fftshift(fft(recon[:,:,:,i_ch,1]),3)
    end

    # # Back to k-space
    # calibration = recon[:,:,:,:,1]
    # calibration = fft(calibration,1)
    # calibration = fft(calibration,2)
    # calibration = fftshift(fft(calibration,3),3)

    calibration = calibration[Int(end/2+1-calib_size):Int(end/2+calib_size),Int(end/2+1-calib_size):Int(end/2+calib_size),Int(end/2+1-6):Int(end/2+6),:]
    SensitivityMap = espirit(calibration,MtxSize)
    SensitivityMap = fftshift(fftshift(SensitivityMap,1),2)

    # ToDo: Sometimes I need to do fftshift in 3 dim also
    SensitivityMap = fftshift(SensitivityMap,3)

    SensitivityMap = dropdims(SensitivityMap; dims = 5)

    return SensitivityMap
    
end

"""
CalculateOffresonanceMap(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Calculates B0 map
"""
function CalculateOffresonanceMap(recon_b0,SensitivityMap,EchoTimes::Vector{Float64})

    mask = SensitivityMap[:,:,:,1]
    mask[abs.(mask) .> 0] .= 1
    mask = imresize(mask,(size(recon_b0)[1],size(recon_b0)[2],size(recon_b0)[3]))
    mask = isone.(mask)

    # Lowering resolution of sensitivities
    smap = imresize(SensitivityMap,size(recon_b0)[1:4])

    # Taking only the first 3 echos
    ydata = recon_b0[:,:,:,:,1:3] .*1e3
    echotime = EchoTimes[1:3].*1e-3 *1s   # Original

    # ToDo: Do I need to flip diemension?
    # ydata = reverse(ydata, dims=1)

    yik_sos = sum(conj(smap) .* ydata; dims=4) # coil combine
    yik_sos = yik_sos[:,:,:,1,:] # (dims..., ne)

    (yik_sos_scaled, scale) = b0scale(yik_sos, echotime) # todo
    yik_scale = ydata / scale # original

    finit = b0init(ydata, echotime; SensitivityMap)

    @info ("Stop... B0 map...")
    @infiltrate

    # Original "Low smoothing":
    fmap_run = (niter, precon, track; kwargs...) ->
        b0map(yik_scale, echotime; finit, smap, mask,
        order=1, l2b=0.002, gamma_type=:PR, niter, precon, track, kwargs...)

    # # A lot of smoothing:
    # fmap_run = (niter, precon, track; kwargs...) ->
    #     b0map(yik_scale, echotime; finit, smap, mask,
    #     order=1, l2b=4, gamma_type=:PR, niter, precon, track, kwargs...)


    function runner(niter, precon; kwargs...)
        (fmap, times, out) = fmap_run(niter, precon, true; kwargs...)
        return (fmap, out.fhats, out.costs, times)
    end;
    if !@isdefined(fmap_cg_d)
        niter_cg_d = 200  # (400)
        (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :diag)
        # (fmap_cg_d, fhat_cg_d, cost_cg_d, time_cg_d) = runner(niter_cg_d, :ichol)
    end

    b0 = fmap_cg_d

    b0 = b0.*2π.*-1    # Original 
    # b0 = b0.*-1

    # Temp: Trying to rescale
    # ToDo: Are the maps in MRIFieldmaps.jl off by 10? 
    # B0 correction is optimal after rescaling... (and not converting to rada/s) Also values are as expected ~100Hz 
    # This seems to be the case only for sv_01
    # b0 = b0.*10

    b0 = ustrip.(b0)

    b0 = imresize(b0,size(SensitivityMap)[1:3],method=BSpline(Cubic()))
    b0 = 1im.*b0
    b0 = convert(Array{ComplexF32,3},b0)

    @info ("Stop... B0 map...")
    @infiltrate

    return b0
end