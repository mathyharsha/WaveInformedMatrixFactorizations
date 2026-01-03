function Y = inv_spectral_wave_transform(yhat,V_x,V_y,V_t)

    % yhat is vector
    % Y is tensor
    
    data_temp_trans = reshape(yhat,size(V_y,1),size(V_x,1),size(V_t,1));
    
    data_temp_trans = reshape(data_temp_trans,size(V_y,1)*size(V_x,1),size(V_t,1));
    
    data_temp = data_temp_trans*transpose(V_t);
    
    data_hat = reshape(data_temp,size(V_y,1),size(V_x,1),size(V_t,1));
    
    Y = data_hat;
    
    for i = 1:size(V_t,1)
        Y(:,:,i)= V_y*Y(:,:,i)*transpose(V_x);
    end
    
    
end