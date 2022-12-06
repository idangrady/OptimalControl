function [u,J] = dptv( M, C, T)
h = length(M);
for l = 1:length(T)
    J_{h+1}{l} = T(l);
end
for k = h:-1:1
    ni = length(M{k}); % state dimension
    for i=1:ni
        nj = length(M{k}{i});
        caux = zeros(1,nj);
        for j = 1:nj
            caux(j) = C{k}{i}{j} + J_{k+1}{M{k}{i}{j}};
        end
        [a,b] = sort(caux);
        J_{k}{i} = a(1); J{k}{i} = J_{k}{i};
        u{k}{i}(1) = b(1);
        for ell = 2:length(a)
            if abs( a(ell) - a(1))<1e-8
                u{k}{i}(ell) = b(ell);
            else
                break;
            end
        end
    end
end
end