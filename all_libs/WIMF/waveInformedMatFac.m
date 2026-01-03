function [Dc,x,C] = waveInformedMatFac(data,Trans,V_x,V_y,V_t,varargin)

    if all(varargin{1}=='count')
        thres = varargin{2}; 
    end
    
    wave_data = data;

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

    threshold = 0.1;

    polar = 100;

    count = 0;
    alpha = 3; 

while polar>1+threshold
    
    count = count+1;
    % if isempty(Dc)
    %     z = data_hat;
    % else
    %     z = data_hat-Dc*x;
    % end
    
    tol = 0.1;
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
    grad_desc = false;
    if all(varargin{3}=='gradient_descent')
        grad_desc = varargin{4}; 
    end
    if grad_desc==true 
    alpha = 3; 
    Dnew = zeros(size(Dc));
    xnew = zeros(size(x));
    %cnew = zeros(size(C,1));
    cnew = 1./sqrt(C);
    i = 1;
    gradient_steps = 0;

    while (norm(Dnew-Dc)/norm(Dc)) + (norm(xnew-x)/norm(xnew)) > 2*size(Dc,2)*tol
       if i > 1   
           Dc = Dnew;
           x = xnew;
           C = 1./sqrt(cnew);
       end

       while 1

       summer = 0;

        for j = 1:size(Dc,2)
           summer = summer + norm(  full( ( gamma*((A1+A2) - (1/C(j,j)^2)*A3 )  )*Dc(:,j) )  )^2; 
        end

        objall = 0.5*norm( data_hat - spectralMask(Dc*x,Trans,V_x,V_y,V_t)  )^2 + 0.5*lambda*(gamma*summer + norm(Dc,'fro')^2 + norm(x)^2 );

           for j = 1:size(Dc,2)

               % Dnew(:,j) = Dc(:,j) - alpha* ( 2*x(j)*x(j)*Dc(:,j) - x(j)*Dc*x + lambda* full( ( gamma*((A1+A2) - (1/C(j,j)^2)*A3 )  )*Dc(:,j) ) + lambda*Dc(:,j) );
               % cnew(j,j) = ((1/C(j,j))^2 - alpha*lambda*gamma*( ((1/C(j,j))^2)*Dc(:,j)'*full( (A1+A2).^2 *Dc(:,j)) -Dc(:,j)'* full( (A1+A2)*A3*Dc(:,j)  )  ))^(-0.5); 
               % C(j,j) = 1/sqrt(cnew(j,j));
               Dnew(:,j) = Dc(:,j) - alpha*( -x(j)*data_hat + x(j)*spectral_wave_transform(Trans.*inv_spectral_wave_transform(Dc*x,V_x,V_y,V_t),V_x,V_y,V_t) + lambda*full( (A + ( gamma*((A1+A2) - (1/C(j,j)^2)*A3 ).^2)  )*Dc(:,j)) );
               cnew(j,j) = (Dnew(:,j)'* ( (A1+A2)*A3 )*Dnew(:,j)) / (Dnew(:,j)'* ( A3.^2) * Dnew(:,j));
               C(j,j) = 1/sqrt(cnew(j,j));
               gradient_steps = gradient_steps + 1;
           end
        y_masked = spectralMask(Dnew*x,Trans,V_x,V_y,V_t);
        xnew = x - alpha*( Dnew'*y_masked - Dnew'*data_hat + lambda*x   );

        summer = 0;
        for j = 1:size(Dc,2)
            summer = summer + norm(  full( ( gamma*((A1+A2) - (1/C(j,j)^2)*(A3) )  )*Dnew(:,j) )  )^2; 
        end

        objnew = 0.5*norm( data_hat - spectralMask(Dnew*xnew,Trans,V_x,V_y,V_t)  )^2 + 0.5*lambda*(gamma*summer + norm(Dnew,'fro')^2 + norm(xnew)^2 );
        disp('---')
        disp('Gradient steps taken')
        disp(gradient_steps)
        disp('Learning rate ' )
        disp(alpha)
        disp('How bad is the reduce: ')
        disp((norm(Dnew - Dc)/norm(Dc)) + (norm(xnew-x)/norm(x)))
        disp('Present Loss Function Value')
        disp(objall)
        disp('Increase caused by higher learning rate')
        disp(objnew)
        %disp(objnew<=objall)
        disp('Decreasing learning rate...')
        
        disp('---')
        if (objnew<=objall) || (((norm(Dnew - Dc)/norm(Dc))+((norm(xnew-x)/norm(x)))) < 2*size(Dnew,2)*tol)       
         break; 
        end

        alpha = alpha/1.3;


       end


       i = i+1; 
    end

    % Dc = Dnew;
    % x = xnew;
    end
    if count>thres-1
        break
    end
    
end
    
end