function yhat = spectral_wave_transform(new_data,V_x,V_y,V_t)
    
    % Y is a tensor
    % yhat is a vector
    
    data_hat = new_data;

    for i = 1:size(new_data,3)
        data_hat(:,:,i)= transpose(V_y)*data_hat(:,:,i)*V_x;
    end
    
    data_temp = reshape(data_hat,size(data_hat,1)*size(data_hat,2),size(data_hat,3));
    
    data_temp_trans = data_temp*V_t;
    
    data_temp_trans = reshape(data_temp_trans,size(data_hat,1),size(data_hat,2),size(data_hat,3));
    
    yhat = reshape(data_temp_trans,numel(data_temp_trans),1);

end