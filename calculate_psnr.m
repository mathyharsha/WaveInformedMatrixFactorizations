function psnr_value = calculate_psnr(original, reconstructed, peak_value)
    % Calculate PSNR between original and reconstructed signals/images
    %
    % Inputs:
    %   original - reference signal/image
    %   reconstructed - reconstructed signal/image
    %   peak_value - maximum possible pixel value (optional)
    %                defaults to max(original(:)) if not provided
    %
    % Output:
    %   psnr_value - PSNR in dB
    
    if nargin < 3
        peak_value = max(original(:));
    end
    
    % Calculate Mean Squared Error (MSE)
    mse = mean((original(:) - reconstructed(:)).^2);
    
    % Handle perfect reconstruction case
    if mse == 0
        psnr_value = Inf;
        return;
    end
    
    % Calculate PSNR in decibels
    psnr_value = 10 * log10(peak_value^2 / mse);
end