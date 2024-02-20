function [x] = upperTriangularSolve(a, b, c, v)
% Projekt 2, zadanie 13
% Piotr Jacak, 327354

% Funkcja pomocnicza, rozwiązująca układ równań Ax=v, gdzie A jest
% macierzą górną trójkątną z trzema przekątnymi, 
% jedną główną i dwiema nad nią.
% WEJŚCIE:
%   a - wektor elementów głównej przekątnej macierzy A
%   b - wektor elementów przekątnej nad główną przekątną macierzy A
%   c - wektor elementów przekątnej nad przekątną b macierzy A
%   v - wektor elementów prawej strony
% WYJŚCIE:
%   x - wektor rozwiązania

n = length(a);

x(n) = v(n) / a(n);
x(n-1) = (v(n-1) - b(n-1) * x(n)) / a(n-1);

for k = n-2: -1: 1
    x(k) = (v(k) - c(k) * x(k+2) - b(k) * x(k+1)) / a(k);
end
x = x';

end % function