function [gx_spoil,gy_spoil,gz_spoil] = prepare_spoilers(params,mag_spoil,sr_spoil)
    lims = params.gen.lims;
    
    gx_spoil=mr.makeTrapezoid('x','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);
    gy_spoil=mr.makeTrapezoid('y','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);
    gz_spoil=mr.makeTrapezoid('z','maxGrad',mag_spoil*lims.gamma,'maxSlew',sr_spoil*lims.gamma,'Area',params.gen.del_k(1)*params.gen.n(1)*1.5);

end