function gz_blips = prepare_gz_blips(params)
lims = params.gen.lims;    
% Kz Blips:
    tmp = params.gen.n(3);
    if mod(tmp,2) == 1; tmp = tmp + 1; end
    for i=1:tmp
        if mod(params.gen.n(3),2) == 1
            area = -(params.gen.del_k(3)*(params.gen.n(3)/2))+(params.gen.del_k(3)*(i-1))+(params.gen.del_k(3)/2);
        else
            area = -(params.gen.del_k(3)*(params.gen.n(3)/2))+(params.gen.del_k(3)*(i-1));
        end
        dur = ceil(2*sqrt(area/lims.maxSlew)/10e-6)*10e-6;
        if area ~= 0
            if i == 1
                gz_blips(i) = mr.makeTrapezoid('z',lims,'Area',area);
            else
                gz_blips(i) = mr.makeTrapezoid('z',lims,'Area',area,'Duration',mr.calcDuration(gz_blips(1)));
            end
        end
    end
    if mod(floor(params.gen.n(3)),2) == 1
        gz_blips(round(params.gen.n(3)/2)) = [];    % Removing the empty blip...
    else
        gz_blips(round(params.gen.n(3)/2)+1) = [];    % Removing the empty blip...
    end
    
    % Reshuffling blips if center-out
    if params.gen.kz_enc == 1
        tmp = [];
        j = floor(params.gen.n(3)/params.gen.kz)/2;
        for i=1:length(gz_blips)
            tmp = [tmp gz_blips(j)];
            if mod(i,2) == 0
                j = j-i;
            else
                j = j+i;
            end
        end
        gz_blips = tmp;
    end

    % Repeting for interlaves
    gz_blips = repmat(gz_blips',[1,params.gen.seg]);

end