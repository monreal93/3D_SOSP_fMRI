
"""
CorrectFFTShift(RawData::RawAcquisitionData,k0::AbstractVector,params::Dict{Symbol,Any})

Performs FFT shift
"""

function CorrectFFTShift(tmp,fftShift::Tuple,ksTraj::Matrix{Float32},params::Dict{Symbol,Any})
    
    ksTraj = ksTraj[:,1:params[:numRead]]'

    tmp = tmp.*exp.(1im*2*π.*ksTraj[:,1].*(fftShift[2]./2))  # kx shift
    tmp = tmp.*exp.(1im*2*π.*ksTraj[:,2].*(fftShift[1]./2))  # ky shift
    # tmp = tmp.*exp.(1im*2*π.*ksTraj[3,:].*fftShift[3])  # kz shift
    
    return tmp
end


"""
CorrectPartitionDORK(RawData::RawAcquisitionData,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function CorrectPartitionDORK(tmp,nav_ref,params::Dict{Symbol,Any})
    
    if params[:kz_enc] == 0  # Linear
        nav = tmp[6:30,Int(params_pulseq["gen"]["n_ov"][3]/2)+1,:,1]
        nav = mean(nav, dims=2)
    elseif params[:kz_enc] == 1 # Center-out
        nav = tmp[6:30,1,:,1]
        nav = mean(nav, dims=2)
    end

    del_phi = mean(angle.(nav)-angle.(nav_ref))/params[:TE]
    
    tmp = tmp.*exp.(-1im.*del_phi.*params[:acq_times]')

    return tmp
end


"""
Correctk0(RawData::RawAcquisitionData,k0::AbstractVector,params::Dict{Symbol,Any})

Performs partition DORK correction
"""

function Correctk0(tmp,k0_sk::AbstractVector{Float64},k0_sim::AbstractMatrix{Float64},params::Dict{Symbol,Any})
    
    @infiltrate
    
    # Un-do Siemens ECC
    tmp = tmp./exp.(1im.*k0_sim)

    # Apply Skope K0
    tmp = tmp.*exp.(1im.*k0_sk.*-2π)

    return tmp
end

