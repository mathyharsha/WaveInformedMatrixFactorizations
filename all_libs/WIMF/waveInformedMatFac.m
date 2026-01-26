function [Dc,x,C] = waveInformedMatFac(data,Trans,V_x,V_y,V_t,varargin)

    % if all(varargin{1}=='count')
    %     thres = varargin{2}; 
    % end
    
    algo_parameters = containers.Map({'count','gradient_descent','tolerance','threshold'},{5,false,0.1,0.1});


    wave_data = data;

    for i=1:2:length(varargin)
        if isKey(algo_parameters,varargin{i})==true
            algo_parameters(varargin{i}) = varargin{i+1};
        end
    end

    thres = algo_parameters('count');

    Lx = 4*del2(eye(size(wave_data,2)));
    Ly = 4*del2(eye(size(wave_data,1)));
    Lt = 4*del2(eye(size(wave_data,3)));

    [V_x,Dx] = eig(Lx);
    [V_y,Dy] = eig(Ly);
    [V_t,Dt] = eig(Lt);

    A1 = kron(kron(speye(size(Dt,1)),Dx),speye(size(Dy,1)));
    A1_ = kron(kron(speye(size(Dt,1)),Dx),speye(size(Dy,1)));
    A2 = kron(speye(size(Dt,1)),kron(speye(size(Dx,1)),Dy));
    A2_ = kron(speye(size(Dt,1)),kron(speye(size(Dx,1)),Dy));
    A3 = kron(Dt,kron(speye(size(Dx,1)),speye(size(Dy,1))));
    A3_ = kron(Dt,kron(speye(size(Dx,1)),speye(size(Dy,1))));
    A = speye(size(Dx,1)*size(Dy,1)*size(Dt,1));

    wave_data_ = reshape(wave_data,numel(wave_data),1);
    data_hat = spectral_wave_transform(wave_data,V_x,V_y,V_t);

    vv = diag(A1+A2)./diag(A3);
    vv = full(vv);

    vv = diag(A1+A2)./diag(A3);
    vv = full(vv);

    a = min(abs(vv));
    b = max(abs(vv));

    % gamma = 0.2*(size(wave_data,1)/pi)^2;
    % lambda = 0.1*norm(data_hat);

    gamma = 0.2*(size(wave_data,1)/pi)^2;
    lambda = 0.1*norm(data_hat);

    Dc = [];
    x = [];
    C = [];

    threshold = algo_parameters('threshold');

    polar = 100;

    count = 0;
    alpha = 0.001; 
    gradient_steps = 0;

while polar>1+threshold
    alpha = 5; 
    count = count+1;
    % if isempty(Dc)
    %     z = data_hat;
    % else
    %     z = data_hat-Dc*x;
    % end
    
    %tol = 0.1;
    %tol = 10;
    tol = algo_parameters('tolerance');

    if isempty(Dc)
        z = data_hat;
    else
        z = spectral_wave_transform(Trans.*inv_spectral_wave_transform(data_hat-Dc*x,V_x,V_y,V_t),V_x,V_y,V_t);
    end

    
    ObjectiveFunction = @(c) -forPolar_actual(c,gamma,A,A1,A2,A3,z);
    tic
    [cbar,funcval] = simulannealbnd(ObjectiveFunction,(1/sqrt(a)+1/sqrt(b))/3,1/sqrt(b),1/sqrt(a));
    % Define bounds
    % lb = 1/sqrt(b);  % Lower bound
    % ub = 1/sqrt(a);  % Upper bound
    % 
    % % Set options (optional)
    % options = optimoptions('surrogateopt', ...
    %     'Display', 'iter', ...
    %     'MaxFunctionEvaluations', 200, ...
    %     'PlotFcn', 'surrogateoptplot');
    % 
    % % Run optimization
    % [cbar, funcval] = surrogateopt(ObjectiveFunction, lb, ub, options);
    toc
    disp(cbar);
    C(count,count) = cbar;
    
    d_temp = (diag((speye(size(Dx,1)*size(Dy,1)*size(Dt,1)) + (gamma)*((A1+A2) - (1/cbar^2)*(A3)).^2)).^(-0.5)).*z;
    d_temp = full((diag((speye(size(Dx,1)*size(Dy,1)*size(Dt,1)) + (gamma)*((A1+A2) - (1/cbar^2)*(A3)).^2)).^(-0.5)).*d_temp)/norm(d_temp);
    
    % if isempty(Dc)
    %     polar = (1/lambda)*d_temp'*data_hat;
    % else
    %     polar = (1/lambda)*d_temp'*(data_hat-Dc*x);
    % end

    polar = (1/lambda)*d_temp'*z;
   
    
    tau = sqrt( d_temp'*z - lambda )/(norm(d_temp));
    
    x = [x;tau];
    
    Dc = [Dc, d_temp*tau];
    
    %polar = (1/lambda)*forPolar(cbar,gamma,Dx,Dt,data_hat-Dc*x);
    
    disp('Iteration:')
    disp(count)
    disp('Value of Polar:')
    disp(polar)
    
    %tic
    % Gradient descent
    grad_desc = algo_parameters('gradient_descent');
    if all(varargin{3}=='gradient_descent')
        grad_desc = varargin{4}; 
    end
    if grad_desc==true 
    % alpha = 5; 
    Dnew = zeros(size(Dc));
    xnew = zeros(size(x));
    %cnew = zeros(size(C,1));
    Cnew = zeros(size(C));
    i = 1;
    %alpha = 1e-3;
    % fig = figure('Name', 'Gradient Descent Progress', 'Position', [100, 100, 800, 600]);
    % loss_history = [];
    % iteration_history = [];
    % (norm(Dnew-Dc)/norm(Dc)) + (norm(xnew-x)/norm(xnew))

    while norm(Dnew-Dc,'fro')/(norm(Dc,'fro')) + norm(xnew-x)/(norm(x)) > 2*size(Dc,2)*(tol)
       if i > 1   
           Dc = Dnew;
           x = xnew;
           C = Cnew;
       end

       while 1

       summer = 0;

        for j = 1:size(Dc,2)
           summer = summer + norm(  full( ( gamma*((A1+A2) - (1/C(j,j)^2)*A3 )  )*Dc(:,j) )  )^2; 
        end

        objall = 0.5*norm( data_hat - spectralMask(Dc*x,Trans,V_x,V_y,V_t)  )^2 + 0.5*lambda*(summer + norm(Dc,'fro')^2 + norm(x)^2 );
        
        % % Store initial loss for this outer iteration
        % loss_history = [loss_history, objall];
        % iteration_history = [iteration_history, i];

           for j = 1:size(Dc,2)

               % Dnew(:,j) = Dc(:,j) - alpha* ( 2*x(j)*x(j)*Dc(:,j) - x(j)*Dc*x + lambda* full( ( gamma*((A1+A2) - (1/C(j,j)^2)*A3 )  )*Dc(:,j) ) + lambda*Dc(:,j) );
               % cnew(j,j) = ((1/C(j,j))^2 - alpha*lambda*gamma*( ((1/C(j,j))^2)*Dc(:,j)'*full( (A1+A2).^2 *Dc(:,j)) -Dc(:,j)'* full( (A1+A2)*A3*Dc(:,j)  )  ))^(-0.5); 
               % C(j,j) = 1/sqrt(cnew(j,j));
               Dnew(:,j) = Dc(:,j) - alpha*( -x(j)*data_hat + x(j)*spectral_wave_transform(Trans.*inv_spectral_wave_transform(Dc*x,V_x,V_y,V_t),V_x,V_y,V_t) + lambda*full( (A + ( gamma*((A1+A2) - (1/C(j,j)^2)*A3 ).^2)  )*Dc(:,j)) );
               %c = (Dnew(:,j)'* ( (A1+A2).*A3 )*Dnew(:,j)) / (Dnew(:,j)'* ( A3.^2) * Dnew(:,j));
               c = ((1/C(j,j))^2 - alpha*lambda*gamma*( ((1/C(j,j))^2)*Dc(:,j)'*full( (A1+A2).^2 *Dc(:,j)) -Dc(:,j)'* full( (A1+A2)*A3*Dc(:,j)  )  ))^(-0.5); 
               Cnew(j,j) = 1/sqrt(c);
               gradient_steps = gradient_steps + 1;
           end
        y_masked = spectralMask(Dnew*x,Trans,V_x,V_y,V_t);
        xnew = x - alpha*( Dnew'*y_masked - Dnew'*data_hat + lambda*x   );

        summer = 0;
        for j = 1:size(Dc,2)
            summer = summer + norm(  full( ( gamma*((A1+A2) - (1/Cnew(j,j)^2)*(A3) )  )*Dnew(:,j) )  )^2; 
        end

        objnew = 0.5*norm( data_hat - spectralMask(Dnew*xnew,Trans,V_x,V_y,V_t)  )^2 + 0.5*lambda*(summer + norm(Dnew,'fro')^2 + norm(xnew)^2 );
        % disp('---')
        % disp('Gradient steps taken')
        % disp(gradient_steps)
        % disp('Learning rate ' )
        % disp(alpha)
        % disp('How bad is the reduce: ')
        % disp((norm(Dnew - Dc)/(norm(Dc))) + (norm(xnew-x)/(norm(x))))
        % disp('Present Loss Function Value')
        % disp(objall)
        % disp('Increase caused by higher learning rate')
        % disp(objnew)
        %disp(objnew<=objall)

        % ((norm(Dnew - Dc)/norm(Dc))+((norm(xnew-x)/norm(x))))
        if (objnew<objall) || ( (norm(Dnew-Dc,'fro')/(norm(Dc,'fro')) + norm(xnew-x)/(norm(x))) < 2*size(Dc,2)*(tol))      
            break; 
        end
        % disp('Decreasing learning rate...')
        % 
        % disp('---')
        alpha = alpha/1.3;

       end


       i = i+1; 
    end

    % Dc = Dnew;
    % x = xnew;
    
    fprintf("Gradient Steps: %d",gradient_steps-int16(log10(0.001/alpha))/log10(1.3))

    end
    
    % % Gradient descent - Modified using Claude
    % grad_desc = false;
    % if all(varargin{3}=='gradient_descent')
    %     grad_desc = varargin{4}; 
    % end
    % 
    % if grad_desc == true 
    %     alpha = 5; 
    %     Dnew = zeros(size(Dc));
    %     xnew = zeros(size(x));
    %     % Dnew = Dc;
    %     % xnew = x;
    %     Cnew = C;  % Fixed: work with C directly, not c²
    %     i = 1;
    %     gradient_steps = 0;
    % 
    %     while (norm(Dnew-Dc,'fro')/norm(Dc,'fro')) + (norm(xnew-x)/norm(x)) > 2*size(Dc,2)*tol
    %         if i > 1   
    %             Dc = Dnew;
    %             x = xnew;
    %             C = Cnew;
    %         end
    % 
    %         % Compute current objective value
    %         y_current = spectral_wave_transform(Trans.*inv_spectral_wave_transform(Dc*x,V_x,V_y,V_t),V_x,V_y,V_t);
    %         residual = data_hat - y_current;
    % 
    %         summer = 0;
    %         for j = 1:size(Dc,2)
    %             summer = summer + norm(full(gamma*((A1+A2) - (1/Cnew(j,j)^2)*A3)*Dc(:,j)))^2; 
    %         end
    %         objall = 0.5*norm(residual)^2 + 0.5*lambda*(gamma*summer + norm(Dc,'fro')^2 + norm(x)^2);
    % 
    %         % Backtracking line search
    %         while true
    %             % GAUSS-SEIDEL: Update each (D_j, c_j) pair sequentially
    %             for j = 1:size(Dc,2)
    %                 % Recompute residual with CURRENT state (using already-updated columns)
    %                 y_current = spectral_wave_transform(Trans.*inv_spectral_wave_transform(Dnew*x,V_x,V_y,V_t),V_x,V_y,V_t);
    %                 residual_current = data_hat - y_current;
    % 
    %                 % Gradient w.r.t. D(:,j) using current residual
    %                 grad_dj = -x(j)*residual_current + lambda*full((A + gamma*((A1+A2) - (1/Cnew(j,j)^2)*A3).^2)*Dnew(:,j));
    % 
    %                 % Update D(:,j)
    %                 Dnew(:,j) = Dc(:,j) - alpha*grad_dj;
    % 
    %                 % Immediately update c_j analytically based on new D(:,j)
    %                 numerator = Dnew(:,j)'*((A1+A2).*A3)*Dnew(:,j);
    %                 denominator = Dnew(:,j)'*(A3.^2)*Dnew(:,j);
    % 
    %                 if abs(denominator) > 1e-12
    %                     Cnew(j,j) = sqrt(numerator / denominator);
    %                 end
    % 
    %                 gradient_steps = gradient_steps + 1;
    %             end
    % 
    %             % Update x using the fully updated D and C
    %             y_new = spectral_wave_transform(Trans.*inv_spectral_wave_transform(Dnew*x,V_x,V_y,V_t),V_x,V_y,V_t);
    %             grad_x = Dnew'*(y_new - data_hat) + lambda*x;
    %             xnew = x - alpha*grad_x;
    % 
    %             % Compute new objective value
    %             y_final = spectral_wave_transform(Trans.*inv_spectral_wave_transform(Dnew*xnew,V_x,V_y,V_t),V_x,V_y,V_t);
    %             residual_new = data_hat - y_final;
    % 
    %             summer_new = 0;
    %             for j = 1:size(Dnew,2)
    %                 summer_new = summer_new + norm(full(gamma*((A1+A2) - (1/Cnew(j,j)^2)*A3)*Dnew(:,j)))^2; 
    %             end
    %             objnew = 0.5*norm(residual_new)^2 + 0.5*lambda*(gamma*summer_new + norm(Dnew,'fro')^2 + norm(xnew)^2);
    % 
    %             % Compute relative change
    %             rel_change_D = norm(Dnew - Dc,'fro') / (norm(Dc,'fro') + 1e-12);
    %             rel_change_x = norm(xnew - x) / (norm(x) + 1e-12);
    %             total_change = rel_change_D + rel_change_x;
    % 
    %             % Display progress
    %             disp('---')
    %             disp(['Gradient steps taken: ' num2str(gradient_steps)])
    %             disp(['Learning rate: ' num2str(alpha)])
    %             disp(['Relative change (D + x): ' num2str(total_change)])
    %             disp(['Previous objective: ' num2str(objall)])
    %             disp(['New objective: ' num2str(objnew)])
    %             disp(['Objective decrease: ' num2str(objall - objnew)])
    % 
    %             % Accept step if objective decreased OR change is small enough
    %             if (objnew <= objall) || (total_change < 2*size(Dnew,2)*tol)
    %                 disp('Step accepted')
    %                 disp('---')
    %                 break; 
    %             end
    % 
    %             % Reduce learning rate and retry
    %             alpha = alpha / 1.3;
    %             disp('Decreasing learning rate...')
    %             disp('---')
    % 
    %             % Safety check: prevent infinite loop
    %             if alpha < 1e-10
    %                 disp('Warning: Learning rate too small, accepting step anyway')
    %                 break;
    %             end
    %         end
    % 
    %         i = i + 1;
    % 
    %         % Safety check: prevent excessive iterations
    %         if i > 1000
    %             disp('Warning: Gradient descent exceeded 1000 iterations')
    %             break;
    %         end
    %     end
    % 
    % end  % end if grad_desc


    if count>thres-1
        break
    end
    
end
    
end