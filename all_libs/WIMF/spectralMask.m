function D_masked_hat = spectralMask(D_hat,Trans,V_x,V_y,V_t)
    
    D_masked_hat = zeros(size(D_hat));
    for i = 1:size(D_hat,2)
        D_masked_hat(:,i) = spectral_wave_transform(Trans.*inv_spectral_wave_transform(D_hat(:,i),V_x,V_y,V_t),V_x,V_y,V_t);
    end

end