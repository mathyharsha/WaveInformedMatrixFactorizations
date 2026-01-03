function d = forPolar_actual(c,gamma,A,A1,A2,A3,z)

%     A1 = kron(kron(speye(size(Dt,1)),Dx),speye(size(Dy,1)));
%     A2 = kron(speye(size(Dt,1)),kron(speye(size(Dx,1)),Dy));
%     A3 = kron(Dt,kron(speye(size(Dx,1)),speye(size(Dy,1))));
    B = A + gamma*((A1+A2) - (1/c^2)*A3)^2;
    d = z'*(B\z);
end