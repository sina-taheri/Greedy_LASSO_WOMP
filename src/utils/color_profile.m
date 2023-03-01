% This function creates a color pool that are discriminative enough for
% plots.
function colors = color_profile(n)
    p = n^(1/3);
    if mod(p, 1) == 0
        n = n + 1;
        p = n^(1/3);
    end
    lp = floor(p);
    up = ceil(p);
    
    a = up;
    b = lp;
    c = lp;
    if a*b*c < n
        b = up;
        if a*b*c < n
            c = up;
        end
    end
        
    colors = [];
    for i = 1:a
        for j = 1:b
            for k = 1:c
                colors = [colors; (i - 1)/a, (j - 1)/b, (k - 1)/c];
            end
        end
    end
    colors = colors(1:n, :);
end