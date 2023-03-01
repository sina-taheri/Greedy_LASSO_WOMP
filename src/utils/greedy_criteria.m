function quantity = greedy_criteria(A, r, w, x, lambda, S, T, type)
switch type
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'wsrlassoGrad'
        quantity = abs(A'*r/norm(r,2) + lambda * w.*sign(x));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'lassol0'
        quantity(T) = min([-abs(A(:,T)'*r).^2 + lambda*w(T).^2, zeros(size(x(T)))],[],2);       
        quantity(S) = (1-(x(S)==0)) .* min([abs(x(S)).^2 - lambda * w(S).^2, zeros(size(x(S)))],[],2);
        quantity = - quantity; % switch sign to maximize it
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'srlassol0'
        quantity(T) = min([sqrt(norm(r)^2 - abs(A(:,T)'*r).^2) - norm(r) ...
            + lambda * w(T).^2, zeros(size(x(T)))],[],2);       
        quantity(S) = (1-(x(S)==0)) .* min([sqrt(norm(r)^2 + abs(x(S)).^2) ...
            - norm(r) - lambda * w(S).^2, zeros(size(x(S)))],[],2);
        quantity = - quantity; % switch sign to maximize it
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'ladlassol0'
        len_T = length(T);
        t = zeros(len_T, 1);
        q = zeros(len_T, 1);
        for k = 1:len_T
            t(k) = lad_solver(A(:, T(k)), r);
            q(k) = norm(r - t(k)*A(:, T(k)), 1) - norm(r, 1) + lambda*w(T(k))^2;
        end
        quantity(T) = min([q, zeros(size(x(T)))],[],2);

        len_S = length(S);
        t = zeros(1, len_S);
        q = zeros(1, len_S);
        for k = 1:len_S
            t(k) = lad_solver(A(:, S(k)), r);
            q(k) = ((x(S(k)) == 0))*min([norm(r - t(k)*A(:, S(k)), 1) - norm(r, 1) + lambda*w(S(k))^2, 0]) + ...
                ((x(S(k)) ~= 0))*min([norm(r - t(k)*A(:, S(k)), 1) - norm(r, 1), ... 
                norm(r + x(S(k))*A(:, S(k)), 1) - norm(r, 1) - lambda*w(S(k))^2]);
        end
        quantity(S) = q;
        quantity = - quantity;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'lassol1'
%         quantity(T) = min([-(abs(A(:, T)'*r) - lambda/2*w(T)).^2, zeros(size(x(T)))],[],2);
%         quantity(S) = (1-(x(S)==0)).* min([(lambda*w(S)/2).^2 + lambda*w(S).*abs(abs(x(S)) - lambda*w(S)/2) - lambda*w(S).*abs(x(S)), ...
%             abs(x(S)).^2 - lambda*w(S).*abs(x(S)), zeros(size(x(S)))],[],2);
%         quantity = - quantity; % switch sign to maximize it
        quantity(T) = min([-(abs(A(:, T)'*r) - lambda*w(T)/2), zeros(size(x(T)))],[],2).^2;
        quantity(S) = -min([-abs(x(S)).*(lambda*w(S) - abs(x(S))), -lambda*w(S).*(abs(x(S)) - lambda*w(S)/4 - abs(abs(x(S)) - lambda*w(S)/2)), ...
            zeros(size(x(S)))],[],2);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'srlassol1'
%         quantity(T) = min([lambda*w(T).*abs(A(:, T)'*r) + sqrt((1-(lambda*w(T)).^2).*(norm(r,2)^2 - abs(A(:, T)'*r).^2)) - norm(r, 2), ...
%             zeros(size(x(T)))],[],2);
%         quantity(S) = (1-(x(S)==0)).* min([zeros(size(x(S))), ...
%             sqrt(abs(x(S)).^2 + norm(r, 2)^2) - norm(r, 2) - lambda*w(S).*abs(x(S)), ... 
%             norm(r, 2)./sqrt(1 - (lambda*w(S)).^2) + lambda*w(S).*abs(abs(x(S)) - lambda*w(S)*norm(r, 2)./sqrt(1 - (lambda*w(S)).^2)) - norm(r, 2) - lambda*w(S).*abs(x(S))],[],2);           
%         quantity = - quantity; % switch sign to maximize it
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'srlassol1' % based on a simpler formulation: ran a test with synthetic data and obtained the exact same result as srlassol1
        quantity(T) = min([lambda*w(T).*abs(A(:, T)'*r) + sqrt((1-(lambda*w(T)).^2).*(norm(r,2)^2 - abs(A(:, T)'*r).^2)) - norm(r, 2), ...
            zeros(size(x(T)))],[],2);
        quantity(S) = min([sqrt(abs(x(S)).^2 + norm(r, 2)^2) - norm(r, 2) - lambda*w(S).*abs(x(S)), ...
            (lambda*w(S) <= 1).*(norm(r, 2)./sqrt(1 - (lambda*w(S)).^2) + lambda*w(S).*abs(abs(x(S)) - lambda*w(S)*norm(r, 2)./sqrt(1 - (lambda*w(S)).^2)) - norm(r, 2) - lambda*w(S).*abs(x(S)))],[],2);           
        quantity = - quantity; % switch sign to maximize it
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'ladlassol1'
        len_T = length(T);
        t = zeros(len_T, 1);
        q = zeros(len_T, 1);
        for k = 1:len_T
            a1 = [A(:, T(k)); lambda*w(T(k))];
            r1 = [r; 0];
            t(k) = lad_solver(a1, r1);
            q(k) = norm(r1 - t(k)*a1, 1) - norm(r, 1);
        end
        quantity(T) = min([q, zeros(size(x(T)))],[],2);

        len_S = length(S);
        t = zeros(1, len_S);
        q = zeros(1, len_S);
        for k = 1:len_S
            a2 = [A(:, S(k)); lambda*w(S(k))];
            r2 = [r; -lambda*w(S(k))*x(S(k))];
            t(k) = lad_solver(a2, r2);
            q(k) = ((x(S(k)) == 0))*min([norm(r2 - t(k)*a2, 1) - norm(r), 0]) + ...
                ((x(S(k)) ~= 0))*min([norm(r2 - t(k)*a2, 1) - norm(r, 1) - lambda*w(S(k))*abs(x(S(k))), ...
                norm(r + x(S(k))*A(:, S(k)), 1) - norm(r, 1) - lambda*w(S(k))*abs(x(S(k)))]);
        end
        quantity(S) = q;
        quantity = - quantity;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'scaling-based'
        quantity(T) = -abs(A(:,T)'*r./w(T));       
        quantity(S) = zeros(size(x(S)));
        quantity = - quantity; % switch sign to maximize it
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
        error('type not supported')
end 