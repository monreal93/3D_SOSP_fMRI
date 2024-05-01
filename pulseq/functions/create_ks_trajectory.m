function ks_traj = create_ks_trajectory(seq,adc,params)

    [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing, slicepos, t_slicepos] = seq.calculateKspacePP();
    j = 1;
    if params.gen.ro_type == 's'
        plane_samples = adc.numSamples;
    elseif params.gen.ro_type == 'c'
        plane_samples = adc.numSamples*params.epi.n_lines;
        % Discarding the EPI navigator samples, here I have 3
        j = j+(adc.numSamples*3);
        tmp = ktraj_adc(1,j:j+adc.numSamples*3);
        tmp_mx = max(tmp(:));
    end
    for i=1:params.gen.n(3)
        l = 1;
        for m=1:params.spi.interl
                ks_traj.kx(l:l+plane_samples-1,i) = ktraj_adc(1,j:j+plane_samples-1);
                ks_traj.ky(l:l+plane_samples-1,i) = ktraj_adc(2,j:j+plane_samples-1);
                ks_traj.kz(l:l+plane_samples-1,i) = ktraj_adc(3,j:j+plane_samples-1);
                j = j+(plane_samples*params.gen.echos);
                l = l+plane_samples;
                if params.gen.dork; j = j+adc_post.numSamples;  end
                if params.gen.ro_type == 'c'
                    j = j+(adc.numSamples*3);
                end
        end
    end
    % Splitting trajectories for ech1 and ech2 of IN-OUT
    % AMM ToDo : Fix this part so it works for more echos..
    if params.spi.type == 3 && params.gen.ro_type == 's'
        % Echo 1 (IN) trajectory
        ks_traj.e1.kx = ks_traj.kx([1:(plane_samples/2)-1,plane_samples+1:end-(plane_samples/2)-1],:);
        ks_traj.e1.ky = ks_traj.ky([1:(plane_samples/2)-1,plane_samples+1:end-(plane_samples/2)-1],:);
        ks_traj.e1.kz = ks_traj.kz([1:(plane_samples/2)-1,plane_samples+1:end-(plane_samples/2)-1],:);
        ks_traj.e1.kx = [padarray(ks_traj.e1.kx(1:end/2,:),[1 0],'post');padarray(ks_traj.e1.kx((end/2)+1:end,:),[1 0],'post')];
        ks_traj.e1.ky = [padarray(ks_traj.e1.ky(1:end/2,:),[1 0],'post');padarray(ks_traj.e1.ky((end/2)+1:end,:),[1 0],'post')];
        ks_traj.e1.kz = [padarray(ks_traj.e1.kz(1:end/2,:),[1 0],'post');padarray(ks_traj.e1.kz((end/2)+1:end,:),[1 0],'post')];
         % Echo 2 (OUT) trajectory
        ks_traj.e2.kx = ks_traj.kx([(plane_samples/2):plane_samples-2,plane_samples+(plane_samples/2):end-2],:);
        ks_traj.e2.ky = ks_traj.ky([(plane_samples/2):plane_samples-2,plane_samples+(plane_samples/2):end-2],:);
        ks_traj.e2.kz = ks_traj.kz([(plane_samples/2):plane_samples-2,plane_samples+(plane_samples/2):end-2],:);
        ks_traj.e2.kx = [padarray(ks_traj.e2.kx(1:end/2,:),[1 0],'pre');padarray(ks_traj.e2.kx((end/2)+1:end,:),[1 0],'pre')];
        ks_traj.e2.ky = [padarray(ks_traj.e2.ky(1:end/2,:),[1 0],'pre');padarray(ks_traj.e2.ky((end/2)+1:end,:),[1 0],'pre')];
        ks_traj.e2.kz = [padarray(ks_traj.e2.kz(1:end/2,:),[1 0],'pre');padarray(ks_traj.e2.kz((end/2)+1:end,:),[1 0],'pre')];
    end
    
    % Plotting traj partition by partition
    figure(17);
    plot3(ks_traj.kx(:),ks_traj.ky(:),ks_traj.kz(:),'DisplayName',sprintf('Kz = %i',i)); title('3D K-space'); xlabel('Kx (1/m)'); ylabel('Ky (1/m)'); zlabel('Kz (1/m)');
    hold on
    view(2)

%     % Plotting trajectory by segment
%     figure(15);
%     hold on
%     adcSamples = params.gen.adc_segments*params.gen.adc_split;
%     for i = 1:params.gen.n(3)
%         idx=1;
%         for j=1:params.spi.interl                                    
%             plot3(ks_traj.kx(idx:idx+adcSamples-1,i),ks_traj.ky(idx:idx+adcSamples-1,i),ks_traj.kz(idx:idx+adcSamples-1,i)); title('3D K-space');
%             idx = idx+adcSamples;
%         end
%     end
%     view(2)

end