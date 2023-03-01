function [x_vals, res, norms, stat] = womp(A, y, w, lambda, n_iter, type, cvx_opt)

N = size(A,2);
for j = 1:N
    if abs(norm(A(:,j),2) - 1) > 1e-12
        error('The columns of A must have unit length')
    end
end

if ~exist('cvx_opt','var')
    cvx_opt = [];
end
cvx_opt = set_options(cvx_opt);

w = w(:); % weights must form a column vectors

% Initialize data structures
x_vals = zeros(N,n_iter+1);
res = zeros(1,n_iter+1); 
res(1) = norm(y); % first residual norm is ||y||_2 since x_0 = 0
norms = zeros(1,n_iter+1);
stat = cell(1,n_iter+1);
stat{1} = 'Solved';
support = [];

switch type
    case 'wlasso'
        cvx_begin quiet
            cvx_precision high  % apparently setting cvx precisoin to best with the solver set to mosek creates a conflict, see CVX forums (it is fine with CVX alone though).
            variable x(N)
            minimize(sum_square_abs(A*x - y) + lambda*norm(w.*x, 1))
        cvx_end
%         [~, ii] = sort(abs(x), 'descend');
%         x(setdiff(1:N, ii(1:n_iter))) = 0;
        x_vals = repmat(x, [1, n_iter + 1]);
        
    case 'wsrlasso'
        cvx_begin quiet
            cvx_precision high
            variable x(N)
            minimize(norm(A*x - y, 2) + lambda*norm(w.*x, 1))
        cvx_end
%         [~, ii] = sort(abs(x), 'descend');
%         x(setdiff(1:N, ii(1:n_iter))) = 0;
        x_vals = repmat(x, [1, n_iter + 1]);
        
    case 'wladlasso'
        cvx_begin quiet
            cvx_precision high
            variable x(N)
            minimize(norm(A*x - y, 1) + lambda*norm(w.*x, 1))
        cvx_end
%         [~, ii] = sort(abs(x), 'descend');
%         x(setdiff(1:N, ii(1:n_iter))) = 0;
        x_vals = repmat(x, [1, n_iter + 1]);
        
    otherwise
        for iter = 2 : n_iter+1   
            r = y - A * x_vals(:,iter-1);  % current residual
            search_set = 1:N; 

            % define quantity to be maximized for the new index search
            S = support;      
            T = setdiff(search_set, S); % complement of S
            x = x_vals(:,iter-1);     % previous approximation (row vector)

            quantity = greedy_criteria(A, r, w, x, lambda, S, T, type);

            [~, sel_index] = max(quantity(search_set)); % find new index
            
            % if the new selected index already exists in the support,
            % there is no need to implement the least-squares/lad step.
%             if ~ismember(search_set(sel_index), support)
                support = union(support, search_set(sel_index)); % update support

                if isequal(type, 'ladlassol0') || isequal(type, 'ladlassol1')
                    x_S = lad_solver(A(:, support), y, 'cvx');
                else
                    x_S = A(:,support) \ y;
                end
%             end

            status = 'Solved';

            x_vals(support, iter) = x_S;

            res(iter) = norm(y - A*x_vals(:, iter),2);
            norms(iter) = norm(x_vals(:, iter).*w, 1);
            stat{iter} = status;        
        end
end

